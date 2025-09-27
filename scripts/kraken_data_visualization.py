import sys
import pandas as pd
import matplotlib.pyplot as plt
import textwrap
import os

# Usage:
# python kraken_data_visualization.py <path_to_kraken_report_file> <visualization_file> [comma_separated_taxids_to_exclude]
if len(sys.argv) < 3:
    print("Usage: python kraken_data_visualization.py <path_to_kraken_report_file> <visualization_file> [taxid1,taxid2,...]")
    sys.exit(1)

report_path = sys.argv[1]
fig_path = sys.argv[2]
exclude_ids = set()
if len(sys.argv) >= 4 and sys.argv[3].strip():
    exclude_ids = set(int(x) for x in sys.argv[3].split(",") if x.strip())

# Load the report file
df = pd.read_csv(
    report_path,
    sep='\t',
    header=None,
    names=['Percent', 'Reads', 'Direct_Assigned', 'Rank', 'TaxID', 'Name']
)

# Normalize types
df['Percent'] = pd.to_numeric(df['Percent'], errors='coerce')
df['TaxID'] = pd.to_numeric(df['TaxID'], errors='coerce').fillna(-1).astype(int)

# Keep Unclassified, Family, Genus, Species
df = df[df['Rank'].str.startswith(('U', 'F', 'G', 'S'))]

# Keep leaf nodes with actual classified reads
leaves = df[df['Direct_Assigned'] > 0].copy()

# Remove duplicates by TaxID
leaves = leaves.drop_duplicates(subset='TaxID')

# Recalculate percentages using Direct_Assigned
total_assigned = leaves['Direct_Assigned'].sum()
if total_assigned == 0:
    print("No directly assigned reads found.")
    sys.exit(0)

leaves['Percent'] = 100.0 * leaves['Direct_Assigned'] / total_assigned

# Drop leaf nodes with recalculated Percent < 1%
leaves = leaves[leaves['Percent'] >= 1]

# Calculate remainder and add 'Others' if needed
percent_sum = leaves['Percent'].sum()
remainder = 100.0 - percent_sum
if remainder > 0.01:
    other_row = {
        'Percent': remainder,
        'Reads': 0,
        'Direct_Assigned': 0,
        'Rank': 'O',
        'TaxID': -1,
        'Name': 'Others'
    }
    leaves = pd.concat([leaves, pd.DataFrame([other_row])], ignore_index=True)

# Clean names
leaves['Name'] = leaves['Name'].str.lstrip()
leaves['Name'] = leaves['Name'].replace('unclassified', 'Unclassified')

# Sort by percent
leaves.sort_values(by='Percent', ascending=False, inplace=True)

# Map each taxid under species to its species-level anchor taxid (the row with Rank == 'S')
species_anchor_for_taxid = {}
species_anchor = None
for _, row in df.iterrows():
    r = str(row['Rank'])
    tx = int(row['TaxID'])
    if r == 'S':  # species line (the anchor)
        species_anchor = tx
        species_anchor_for_taxid[tx] = tx
    elif r.startswith('S'):  # S1, S2, S3... belong to the current species anchor
        if species_anchor is not None:
            species_anchor_for_taxid[tx] = species_anchor
    else:
        # leaving the species block (e.g., back to G/F/U or other ranks)
        species_anchor = None


# print top 10 species >=1% not in exclude_ids or in the same species as any exclude_id
species_candidates = leaves[leaves['Rank'].str.startswith('S')].copy()

# figure out which species anchors should be excluded
excluded_species_anchors = {
    species_anchor_for_taxid.get(int(e), int(e)) for e in exclude_ids
}

# add each row's species anchor (default to itself if not found)
species_candidates['species_anchor'] = species_candidates['TaxID'].map(
    lambda t: species_anchor_for_taxid.get(int(t), int(t))
)

# drop rows whose species anchor is excluded
species_candidates = species_candidates[
    ~species_candidates['species_anchor'].isin(excluded_species_anchors)
]

species_top10 = species_candidates.sort_values(by='Percent', ascending=False).head(10)

if not species_top10.empty:
    header = (
        "\n\n\nNOTE:\n"
        "The following pathogens have >1% assigned reads but were NOT included for detailed analysis with WEPP:"
    )
    lines = [
        f"{i}) {name} (TaxID {tx}; {pct:.2f}%)"
        for i, (name, tx, pct) in enumerate(
            zip(species_top10['Name'], species_top10['TaxID'], species_top10['Percent']),
            start=1
        )
    ]
    print(header + "\n" + "\n".join(lines))
    print("\n\n")
else:
    print(
        "\n\nNOTE:\n"
        "No other pathogens have >=1% of directly assigned reads after exclusions."
    )


fig, ax = plt.subplots(figsize=(10, 8))
ax.set_position([0.05, 0.1, 0.5, 0.8])  # [left, bottom, width, height]

colors = plt.cm.tab20.colors

# Function to wrap long labels
def wrap_labels(labels, width=25):
    return ['\n'.join(textwrap.wrap(label, width)) for label in labels]

# Add percentages to the labels for legend
wrapped_labels_with_percent = [
    f"{name} ({percent:.2f}%)"
    for name, percent in zip(leaves['Name'], leaves['Percent'])
    if round(percent, 2) > 0
]
wrapped_labels = wrap_labels(wrapped_labels_with_percent)

# Pie chart
wedges, texts = ax.pie(
    leaves['Percent'],
    labels=None,
    colors=colors[:len(leaves)],
    startangle=140
)

# Legend
fig.legend(
    wedges,
    wrapped_labels,
    title="Pathogens",
    title_fontsize=16,
    loc="center right",
    bbox_to_anchor=(0.95, 0.5),
    labelspacing=1.2,
    fontsize=14
)

ax.axis('equal')
fig.suptitle('Classification Proportions by Pathogen', fontsize=20, fontweight='bold', ha='center', y=0.85)

os.makedirs(os.path.dirname(fig_path), exist_ok=True)
plt.savefig(fig_path, dpi=300, bbox_inches='tight')
