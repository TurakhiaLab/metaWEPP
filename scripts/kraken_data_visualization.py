import sys
import pandas as pd
import matplotlib.pyplot as plt
import os

# Get arguments
report_path = sys.argv[1]

# Determine next filename
existing = [f for f in os.listdir(output_dir) if f.startswith("classification_proportions_") and f.endswith(".png")]
nums = [int(f.split("_")[-1].replace(".png", "")) for f in existing if f.split("_")[-1].replace(".png", "").isdigit()]
next_num = max(nums) + 1 if nums else 1

# Define output path
outfile = os.path.join(output_dir, f"classification_proportions_{next_num}.png")

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

# Plot pie chart
plt.figure(figsize=(10, 10))
colors = plt.cm.tab20.colors
wedges, texts, autotexts = plt.pie(
    leaves['Percent'],
    labels=None,
    autopct='%1.1f%%',
    colors=colors[:len(leaves)],
    startangle=140
)

# Make the percentage labels bold
for autotext in autotexts:
    autotext.set_fontweight('bold')

# Add cleaned labels to legend
plt.legend(
    wedges,
    leaves['Name'],
    title="Pathogens",
    loc="center left",
    bbox_to_anchor=(1, 0.5),
    labelspacing=1.2
)

plt.title('Classification Proportions by Pathogen')
plt.tight_layout()
plt.savefig("classification_proportions.png", dpi=300)