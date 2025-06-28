import sys
import pandas as pd
import matplotlib.pyplot as plt

# Get file path from command-line argument
if len(sys.argv) < 2:
    print("Usage: python kraken_data_visualization.py <path_to_kraken_report_file>")
    sys.exit(1)

report_path = sys.argv[1]

# Load the report file
df = pd.read_csv(report_path, sep='\t', header=None,
                 names=['Percent', 'Reads', 'Direct_Assigned', 'Rank', 'TaxID', 'Name'])

# Filter only leaf nodes with actual classified reads
leaves = df[df['Direct_Assigned'] > 0].copy()

# Remove duplicates
leaves = leaves.drop_duplicates(subset='TaxID')

# Sort by percent
leaves.sort_values(by='Percent', ascending=False, inplace=True)

# Clean up whitespace from pathogen names
leaves['Name'] = leaves['Name'].str.lstrip()

# Upper case first letter "unclassified"
leaves['Name'] = leaves['Name'].replace('unclassified', 'Unclassified')

fig, ax = plt.subplots(figsize=(10, 8))

# Shrink the pie chart's axes and move it to the left
ax.set_position([0.05, 0.1, 0.5, 0.8])  # [left, bottom, width, height]

colors = plt.cm.tab20.colors
wedges, texts, autotexts = ax.pie(
    leaves['Percent'],
    labels=None,
    autopct='%1.1f%%',
    colors=colors[:len(leaves)],
    startangle=140
)

for autotext in autotexts:
    autotext.set_fontweight('bold')

# Add legend to the right
fig.legend(
    wedges,
    leaves['Name'],
    title="Pathogens",
    loc="center right",
    bbox_to_anchor=(1, 0.5),
    labelspacing=1.2,
    fontsize=12
)

# Add centered title across the top
fig.suptitle('Classification Proportions by Pathogen', fontsize=20, fontweight='bold', ha='center')
plt.savefig("classification_proportions.png", dpi=300, bbox_inches='tight')
