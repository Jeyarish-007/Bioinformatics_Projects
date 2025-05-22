import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_indicators
from collections import Counter

# Sample data
data = {
    'Gene Name': ['TP53', 'BRCA1', 'EGFR', 'KRAS', 'PTEN', 'MYC', 'RB1', 'PIK3CA', 'CDKN2A', 'BRAF', 'NF1', 'APC', 'AKT1', 'MET', 'ERBB2', 'NOTCH1', 'SMAD4', 'ARID1A', 'FBXW7', 'CTNNB1'],
    'Metastasis Frequency (%)': [25.3, 12.7, 8.5, 33.2, 5.8, 15.6, 9.2, 18.9, 7.4, 22.1, 11.5, 14.3, 6.7, 17.8, 13.2, 8.9, 10.4, 9.7, 12.6, 16.5],
    'Total Population': [100, 85, 200, 150, 120, 90, 110, 130, 95, 180, 105, 160, 75, 140, 115, 125, 135, 145, 155, 165],
    'Positive Patients': [25, 11, 17, 50, 7, 14, 10, 25, 7, 40, 12, 23, 5, 25, 15, 11, 14, 14, 20, 27],
    'Technique': ['NGS', 'WGS', 'NGS,WES', 'WES', 'WGS', 'NGS', 'WES', 'NGS,WGS', 'WGS,WES', 'NGS,WES', 'WGS', 'NGS,WGS,WES', 'WES', 'NGS', 'WGS,WES', 'NGS', 'WES', 'NGS,WGS', 'WGS', 'NGS,WES']
}

df = pd.DataFrame(data)

# Step 1: Split techniques into boolean columns
technique_df = df['Technique'].str.get_dummies(sep=',').astype(bool)

# Step 2: Combine metadata and technique flags
combined_df = pd.concat([
    df[['Metastasis Frequency (%)', 'Total Population', 'Positive Patients']],
    technique_df
], axis=1)

# Create UpSet-compatible data
upset_data = from_indicators(indicators=technique_df.columns.tolist(), data=combined_df)

# Create UpSet plot
upset = UpSet(upset_data, show_counts=True, sort_by='cardinality')

# Add strip plots for frequency and sample size
upset.add_catplot(value='Metastasis Frequency (%)', kind='strip', color='blue')
upset.add_catplot(value='Total Population', kind='strip', color='red')

# Plot
fig = plt.figure(figsize=(12, 8))
upset.plot(fig=fig)
plt.suptitle('UpSet Plot of Gene Techniques with Metastasis Frequency and Sample Size', y=1.02)

# Custom legend
plt.figtext(0.7, 0.9, "Blue dots: Metastasis Frequency (%)", fontsize=10, color='blue')
plt.figtext(0.7, 0.87, "Red dots: Total Population (sample size)", fontsize=10, color='red')

plt.tight_layout()
plt.show()

# Optional: Count how many genes used each technique
tech_count = {col: sum(technique_df[col]) for col in technique_df.columns}
print("Technique counts across genes:")
print(tech_count)

# Optional: Print which techniques were used per gene
print("\nGenes vs Techniques:")
for idx, row in technique_df.iterrows():
    print(f"{df.loc[idx, 'Gene Name']}: {', '.join(row[row].index)}")
