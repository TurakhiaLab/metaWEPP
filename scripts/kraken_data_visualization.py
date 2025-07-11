import sys
import pandas as pd
import matplotlib.pyplot as plt
import textwrap

# Get file path from command-line argument
if len(sys.argv) < 2:
    print("Usage: python kraken_data_visualization.py <path_to_kraken_report_file>")
    sys.exit(1)

report_path = sys.argv[1]

# Load the report file
df = pd.read_csv(report_path, sep='\t', header=None,
                 names=['Percent', 'Reads', 'Direct_Assigned', 'Rank', 'TaxID', 'Name'])

# Keep Unclassified, Family, genus, and Species
df = df[df['Rank'].str.startswith(('U', 'F', 'G', 'S'))]

# Filter only leaf nodes with actual classified reads
leaves = df[df['Direct_Assigned'] > 0].copy()

# Remove duplicates
leaves = leaves.drop_duplicates(subset='TaxID')

# Sum of all retained percentages
retained_pct = leaves['Percent'].sum()

# Calculate the leftover percentage
other_pct = 100.0 - retained_pct

# Append 'Other' category only if it's non-zero
if other_pct >= 0.01:
    leaves = pd.concat([
        leaves,
        pd.DataFrame([{
            'Percent': other_pct,
            'Reads': 0,
            'Direct_Assigned': 0,
            'Rank': 'Other',
            'TaxID': -1,
            'Name': 'Other'
        }])
    ], ignore_index=True)

# Clean up whitespace from pathogen names
leaves['Name'] = leaves['Name'].str.lstrip()
leaves['Name'] = leaves['Name'].replace('unclassified', 'Unclassified')

# Sort by percent
leaves.sort_values(by='Percent', ascending=False, inplace=True)

fig, ax = plt.subplots(figsize=(10, 8))

# Shrink the pie chart's axes and move it to the left
ax.set_position([0.05, 0.1, 0.5, 0.8])  # [left, bottom, width, height]

colors = plt.cm.tab20.colors

# Function to wrap long labels to two lines
def wrap_labels(labels, width=25):
    return ['\n'.join(textwrap.wrap(label, width)) for label in labels]

# Add percentages to the labels for legend
wrapped_labels_with_percent = [
    f"{name} ({percent:.2f}%)"
    for name, percent in zip(leaves['Name'], leaves['Percent'])
    if round(percent, 2) > 0
]
wrapped_labels = wrap_labels(wrapped_labels_with_percent)

# Pie chart without percentages on slices
wedges, texts = ax.pie(
    leaves['Percent'],
    labels=None,
    colors=colors[:len(leaves)],
    startangle=140
)

# Legend showing labels with percentage
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

# Add centered title across the top
fig.suptitle('Classification Proportions by Pathogen', fontsize=20, fontweight='bold', ha='center', y=0.85)
plt.savefig("classification_proportions.png", dpi=300, bbox_inches='tight')