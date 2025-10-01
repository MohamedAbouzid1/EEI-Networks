#!/usr/bin/env python3
"""
EEI Prediction Results Visualization

This script creates visualizations of the EEI prediction results.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from collections import Counter
import os

# Set the style
plt.style.use('ggplot')
sns.set(font_scale=1.2)

# Create output directory for visualizations
os.makedirs("visualizations", exist_ok=True)

# Load the prediction results
results_file = "/cosybio/project/mabouzid/EEI_networks/EEI-Conservation-main/orthology_based_EEI_prediction/results/predicted_human_eeis.tsv"
df = pd.read_csv(results_file, sep='\t')

print(f"Loaded {len(df)} predicted EEIs")
print("\nBasic Statistics:")
print(f"Average confidence score: {df['confidence'].mean():.4f}")
print(f"Median confidence score: {df['confidence'].median():.4f}")
print(f"Min confidence score: {df['confidence'].min():.4f}")
print(f"Max confidence score: {df['confidence'].max():.4f}")

# 1. Confidence Score Distribution
plt.figure(figsize=(10, 6))
sns.histplot(df['confidence'], bins=20, kde=True)
plt.title('Distribution of Confidence Scores for Predicted EEIs')
plt.xlabel('Confidence Score')
plt.ylabel('Count')
plt.axvline(df['confidence'].mean(), color='red', linestyle='--', label=f'Mean: {df["confidence"].mean():.4f}')
plt.axvline(df['confidence'].median(), color='green', linestyle='-.', label=f'Median: {df["confidence"].median():.4f}')
plt.legend()
plt.tight_layout()
plt.savefig('visualizations/confidence_distribution.png', dpi=300)
plt.close()

# 2. Identity Scores Comparison
plt.figure(figsize=(10, 6))
sns.scatterplot(x='identity1', y='identity2', data=df, alpha=0.7, s=100, hue='confidence')
plt.title('Comparison of Identity Scores for Orthologous Exon Pairs')
plt.xlabel('Identity Score for Exon 1')
plt.ylabel('Identity Score for Exon 2')
plt.grid(True)
plt.tight_layout()
plt.savefig('visualizations/identity_comparison.png', dpi=300)
plt.close()

# 3. Orthology Types
orthology_types = []
for _, row in df.iterrows():
    orthology_types.append(f"{row['type1']}-{row['type2']}")

type_counts = Counter(orthology_types)
labels = list(type_counts.keys())
sizes = list(type_counts.values())

plt.figure(figsize=(10, 8))
plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90, shadow=True)
plt.axis('equal')
plt.title('Distribution of Orthology Types in Predicted EEIs')
plt.tight_layout()
plt.savefig('visualizations/orthology_types_pie.png', dpi=300)
plt.close()

# 4. Coverage Percentages
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
sns.boxplot(y=df['exon1_coverage_percent'])
plt.title('Exon 1 Coverage Percentage')
plt.ylabel('Coverage (%)')
plt.grid(True)

plt.subplot(1, 2, 2)
sns.boxplot(y=df['exon2_coverage_percent'])
plt.title('Exon 2 Coverage Percentage')
plt.ylabel('Coverage (%)')
plt.grid(True)

plt.tight_layout()
plt.savefig('visualizations/coverage_boxplots.png', dpi=300)
plt.close()

# 5. Coverage vs Jaccard Index
plt.figure(figsize=(10, 6))
avg_coverage = (df['exon1_coverage_percent'] + df['exon2_coverage_percent']) / 2
sns.scatterplot(x=avg_coverage, y=df['jaccard_percent'], alpha=0.7, s=100, hue='confidence')
plt.title('Average Coverage vs Jaccard Index')
plt.xlabel('Average Coverage (%)')
plt.ylabel('Jaccard Index (%)')
plt.grid(True)
plt.tight_layout()
plt.savefig('visualizations/coverage_vs_jaccard.png', dpi=300)
plt.close()

# 6. Top 20 Proteins
top_proteins1 = Counter(df['protein1'].str.split('_').str[0]).most_common(20)
proteins1_df = pd.DataFrame(top_proteins1, columns=['Protein', 'Count'])

plt.figure(figsize=(12, 8))
sns.barplot(x='Count', y='Protein', data=proteins1_df)
plt.title('Top 20 Proteins for Exon 1')
plt.xlabel('Count')
plt.ylabel('Protein ID')
plt.tight_layout()
plt.savefig('visualizations/top_proteins1.png', dpi=300)
plt.close()

# 7. Confidence vs Coverage
plt.figure(figsize=(10, 6))
sns.scatterplot(x='confidence', y=avg_coverage, data=df, alpha=0.7, s=100, hue='jaccard_percent')
plt.title('Confidence Score vs Average Coverage')
plt.xlabel('Confidence Score')
plt.ylabel('Average Coverage (%)')
plt.grid(True)
plt.tight_layout()
plt.savefig('visualizations/confidence_vs_coverage.png', dpi=300)
plt.close()

# 8. Summary Dashboard
plt.figure(figsize=(16, 12))

# Confidence Distribution
plt.subplot(2, 2, 1)
sns.histplot(df['confidence'], bins=20, kde=True)
plt.title('Confidence Score Distribution')
plt.xlabel('Confidence Score')
plt.ylabel('Count')

# Identity Comparison
plt.subplot(2, 2, 2)
sns.scatterplot(x='identity1', y='identity2', data=df, alpha=0.7, s=50, hue='confidence')
plt.title('Identity Score Comparison')
plt.xlabel('Identity Score for Exon 1')
plt.ylabel('Identity Score for Exon 2')

# Coverage vs Jaccard
plt.subplot(2, 2, 3)
sns.scatterplot(x=avg_coverage, y=df['jaccard_percent'], alpha=0.7, s=50, hue='confidence')
plt.title('Coverage vs Jaccard Index')
plt.xlabel('Average Coverage (%)')
plt.ylabel('Jaccard Index (%)')

# Orthology Types
plt.subplot(2, 2, 4)
sorted_types = sorted(type_counts.items(), key=lambda x: x[1], reverse=True)
types = [t[0] for t in sorted_types]
counts = [t[1] for t in sorted_types]
sns.barplot(x=counts, y=types)
plt.title('Orthology Types Distribution')
plt.xlabel('Count')
plt.ylabel('Orthology Type')

plt.tight_layout()
plt.savefig('visualizations/summary_dashboard.png', dpi=300)
plt.close()

print("\nVisualizations created in the 'visualizations' directory.")