import pandas as pd # pyright: ignore[reportMissingModuleSource]
import matplotlib.pyplot as plt
import os

# Ensure results directory exists
os.makedirs("results", exist_ok=True)

# Load the metrics CSV
df = pd.read_csv("results/sequence_metrics.csv")

print("Loaded Data for Comparison:")
print(df[['organism', 'length', 'gc_percentage']])

# Create comparative GC percentage chart
plt.figure(figsize=(10, 6))
bars = plt.bar(df['organism'], df['gc_percentage'], color=['#4C72B0', '#55A868', '#C44E52', '#8172B3', '#CCB974'])

plt.xlabel('Bacterial Species', fontsize=12)
plt.ylabel('GC Content (%)', fontsize=12)
plt.title('Comparative GC Content of 16S rRNA Genes Across Bacterial Species', fontsize=14, fontweight='bold')
plt.xticks(rotation=15, ha='right')
plt.ylim(0, 100)

# Add value labels on top of bars
for bar in bars:
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2., height + 1, f'{height}%', ha='center', va='bottom', fontsize=10)

plt.tight_layout()

# Save chart
chart_path = "results/gc_comparison.png"
plt.savefig(chart_path)
print(f"\nSuccessfully generated and saved comparison chart to {chart_path}!") 