import sys
import pandas as pd
import matplotlib.pyplot as plt
import textwrap
import os
import math
import subprocess

# Usage:
# python kraken_data_visualization.py <path_to_kraken_report_file> <visualization_file> 
if len(sys.argv) < 4:
    print("Usage: python kraken_data_visualization.py <path_to_kraken_report_file> <visualization_file> <minimum_proportion_threshold_for_wepp>")
    sys.exit(1)

report_path = sys.argv[1]
fig_path = sys.argv[2]
min_proportion_threshold = float(sys.argv[3])
exclude_ids = set()

# Load the report file
df = pd.read_csv(
    report_path,
    sep='\t',
    header=None,
    names=['Percent', 'Reads', 'Direct_Assigned', 'Rank', 'TaxID', 'Name']
)
df['Name'] = df['Name'].str.strip()

# Get total reads from the 'root' entry in the report.
try:
    mapped_reads = df[df['Name'] == 'root']['Reads'].iloc[0]
    unclassified_reads = df[df['Name'] == 'unclassified']['Reads'].iloc[0]
    total_reads = mapped_reads + unclassified_reads
except IndexError:
    print("Error: 'root' or 'unclassified' entry not found in Kraken report. Cannot determine total reads.", file=sys.stderr)
    sys.exit(1)

if total_reads == 0:
    print("No reads found in the sample.", file=sys.stderr)
    sys.exit(1)

# Normalize types
df['Percent'] = pd.to_numeric(df['Percent'], errors='coerce')
df['TaxID'] = pd.to_numeric(df['TaxID'], errors='coerce').fillna(-1).astype(int)

# Keep Unclassified, Species, and Genus for pie chart leaves
df_filtered = df[df['Rank'].str.startswith(('U','S','G'))]

# Keep leaf nodes with actual classified reads
leaves = df_filtered[df_filtered['Direct_Assigned'] > 0].copy()

# Readjust the Genus reads to its species -> Account for Kraken2 errors
leaves = leaves.reset_index(drop=True)

rows_to_drop = []
adjustments = {}
i = 0
while i < len(leaves):
    row = leaves.iloc[i]
    rank = row['Rank']

    # Detect genus by prefix 'G'
    if isinstance(rank, str) and rank.startswith('G'):

        genus_reads = row['Direct_Assigned']
        j = i + 1
        species_indices = []

        # Find consecutive species under this genus
        while j < len(leaves):
            next_rank = leaves.iloc[j]['Rank']
            if isinstance(next_rank, str) and next_rank.startswith('S'):
                species_indices.append(j)
                j += 1
            else:
                break

        # CASE: 1 species under genus
        if len(species_indices) == 1:
            sp_idx = species_indices[0]
            adjustments[sp_idx] = adjustments.get(sp_idx, 0) + genus_reads

        # CASE: multiple species
        elif len(species_indices) > 1:
            share = math.floor(genus_reads / len(species_indices))
            for sp_idx in species_indices:
                adjustments[sp_idx] = adjustments.get(sp_idx, 0) + share

        # Remove genus row ALWAYS
        rows_to_drop.append(i)

    i += 1

# Apply adjustments
for idx, add_val in adjustments.items():
    leaves.at[idx, 'Direct_Assigned'] += add_val

# Drop all genus rows
leaves = leaves.drop(rows_to_drop).reset_index(drop=True)


# Recalculate percentages using the total number of reads in the sample
leaves['Percent'] = 100.0 * leaves['Direct_Assigned'] / total_reads

# Drop leaf nodes with recalculated Percent < min_proportion_threshold
leaves = leaves[leaves['Percent'] >= (min_proportion_threshold * 100)]

# Calculate remainder and add 'Others' if needed
percent_sum = leaves['Percent'].sum()
remainder = 100.0 - percent_sum
if remainder >= 1.0:
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
leaves['Name'] = leaves['Name'].replace('unclassified', 'Unclassified')

# Sort by percent
leaves.sort_values(by='Percent', ascending=False, inplace=True)

# print top 10 species >=1% 
species_candidates = leaves[leaves['Rank'].str.startswith('S')].copy()

# Read taxon IDs that have already been added to the pathogen list
added_taxons_path = "data/pathogens_for_wepp/added_taxons.tsv"
added_taxons = set()
if os.path.exists(added_taxons_path):
    with open(added_taxons_path, "r") as f:
        for line in f:
           line = line.strip()
           if not line:
               continue
           parts = line.split("\t")
           taxid = parts[0].strip()   # only column 1
           added_taxons.add(taxid)

if not species_candidates.empty:
    # Filter out species that are already in added_taxons.tsv
    candidates_to_print = species_candidates[~species_candidates['TaxID'].astype(str).isin(added_taxons)]

    if not candidates_to_print.empty:
        header = (
            "\n\n\nNOTE:\n"
            f"The following species have >{min_proportion_threshold * 100}% assigned reads but ABSENT from your haplotype analysis:"
        )
        lines = [
            f"{i}) {name} (TaxID {tx}; {pct:.2f}%)"
            for i, (name, tx, pct) in enumerate(
                zip(candidates_to_print['Name'], candidates_to_print['TaxID'], candidates_to_print['Percent']),
                start=1
            )
        ]
        print(header + "\n" + "\n".join(lines))
        print("\n\n")


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
fig.suptitle('Pathogen Proportions', fontsize=20, fontweight='bold', ha='center', y=0.85)

os.makedirs(os.path.dirname(fig_path), exist_ok=True)
plt.savefig(fig_path, dpi=300, bbox_inches='tight')

# Call add_ref_mat.py to update the database
subprocess.run(
            ["python", "scripts/add_ref_mat.py"],
            check=True
        )