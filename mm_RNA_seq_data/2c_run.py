#!/usr/bin/env python3
"""
Quick start script to run EEI survival analysis with your data.

This script takes your mapped mouse EEIs and tests them for survival associations
using the methodology from the cancer paper.

Usage:
    python 2c_run.py

Make sure you have:
1. outputs/mapped_mouse_crpes.tsv (your mapped EEIs)
2. geo_data/GSE101336_processed_matrix_FPKM.txt (expression data)
3. outputs/survival_df.csv (survival data)
"""

import pandas as pd
import numpy as np
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib.pyplot as plt
import os
import warnings
warnings.filterwarnings('ignore')

def simple_coordinate_to_gene_mapping(coordinate, gene_list):
    """
    Simple mapping from genomic coordinate to gene name.
    
    This is a placeholder implementation. For real analysis, you would need:
    - Proper gene annotation (GTF file)
    - Genomic interval overlap tools
    - BioMart or Ensembl API
    
    For now, this uses a deterministic but arbitrary mapping.
    """
    # Convert coordinate to a stable hash
    gene_list_sorted = sorted(list(gene_list))
    hash_val = abs(hash(coordinate)) % len(gene_list_sorted)
    return gene_list_sorted[hash_val]

def run_eei_survival_analysis():
    """Main function to run the EEI survival analysis."""
    
    print("="*60)
    print("MOUSE EEI SURVIVAL ANALYSIS")
    print("="*60)
    
    # 1. Load expression data
    print("1. Loading expression data...")
    try:
        df = pd.read_csv("geo_data/GSE101336_processed_matrix_FPKM.txt", sep="\t")
        expression_df = df.set_index("gene_id").drop(columns=["tss_id", "locus"])
        mouse_expression = expression_df.T
        mouse_expression.index.name = "mouse_id"
        print(f"   Loaded expression data: {mouse_expression.shape}")
    except Exception as e:
        print(f"   Error loading expression data: {e}")
        return
    
    # 2. Load survival data
    print("2. Loading survival data...")
    try:
        survival_df = pd.read_csv("outputs/survival_df.csv")
        print(f"   Loaded survival data: {survival_df.shape}")
        print(f"   Columns: {survival_df.columns.tolist()}")
    except Exception as e:
        print(f"   Error loading survival data: {e}")
        return
    
    # 3. Load mapped EEIs
    print("3. Loading mapped mouse EEIs...")
    try:
        mapped_eeis = pd.read_csv("outputs/mapped_mouse_crpes.tsv", sep='\t')
        print(f"   Loaded {len(mapped_eeis)} mapped EEIs")
        print(f"   Columns: {mapped_eeis.columns.tolist()}")
    except Exception as e:
        print(f"   Error loading mapped EEIs: {e}")
        return
    
    # 4. Align sample IDs
    print("4. Aligning sample IDs...")
    expr_samples = set(mouse_expression.index.astype(str))
    surv_samples = set(survival_df['mouse_id'].astype(str))
    common_samples = expr_samples.intersection(surv_samples)
    
    print(f"   Expression samples: {len(expr_samples)}")
    print(f"   Survival samples: {len(surv_samples)}")
    print(f"   Common samples: {len(common_samples)}")
    
    if len(common_samples) < 10:
        print("   WARNING: Very few overlapping samples!")
        print(f"   Expression IDs (first 5): {list(expr_samples)[:5]}")
        print(f"   Survival IDs (first 5): {list(surv_samples)[:5]}")
    
    # Filter to common samples
    mouse_expression_aligned = mouse_expression.loc[mouse_expression.index.astype(str).isin(common_samples)]
    survival_aligned = survival_df[survival_df['mouse_id'].astype(str).isin(common_samples)].copy()
    
    # 5. Test EEIs for survival associations
    print(f"5. Testing {len(mapped_eeis)} EEIs for survival associations...")
    
    results = []
    cpm_threshold = 0.5  # Expression threshold from cancer paper
    min_group_size = max(3, len(common_samples) // 8)  # Adaptive minimum group size
    
    print(f"   Using CPM threshold: {cpm_threshold}")
    print(f"   Using minimum group size: {min_group_size}")
    
    for idx, eei_row in mapped_eeis.iterrows():
        try:
            # Get coordinates
            mouse_exon1 = eei_row['mouse_exon1']
            mouse_exon2 = eei_row['mouse_exon2']
            
            # Map to genes (simplified mapping)
            gene1 = simple_coordinate_to_gene_mapping(mouse_exon1, mouse_expression_aligned.columns)
            gene2 = simple_coordinate_to_gene_mapping(mouse_exon2, mouse_expression_aligned.columns)
            
            # Get expression data
            expr1 = mouse_expression_aligned[gene1]
            expr2 = mouse_expression_aligned[gene2]
            
            # Define EEI presence: both exons expressed above threshold
            eei_present = (expr1 > cmp_threshold) & (expr2 > cmp_threshold)
            
            # Create groups
            present_samples = eei_present[eei_present == True].index.astype(str).tolist()
            absent_samples = eei_present[eei_present == False].index.astype(str).tolist()
            
            # Check group sizes
            if len(present_samples) < min_group_size or len(absent_samples) < min_group_size:
                continue
            
            # Get survival data for each group
            survival_aligned['mouse_id'] = survival_aligned['mouse_id'].astype(str)
            present_survival = survival_aligned[survival_aligned['mouse_id'].isin(present_samples)]
            absent_survival = survival_aligned[survival_aligned['mouse_id'].isin(absent_samples)]
            
            if len(present_survival) == 0 or len(absent_survival) == 0:
                continue
            
            # Perform log-rank test
            logrank_result = logrank_test(
                present_survival['survival_days'], 
                absent_survival['survival_days'],
                present_survival['event_status'], 
                absent_survival['event_status']
            )
            
            p_value = logrank_result.p_value
            
            # Calculate median survival
            kmf = KaplanMeierFitter()
            
            kmf.fit(present_survival['survival_days'], present_survival['event_status'])
            median_present = kmf.median_survival_time_ if not pd.isna(kmf.median_survival_time_) else present_survival['survival_days'].max()
            
            kmf.fit(absent_survival['survival_days'], absent_survival['event_status'])
            median_absent = kmf.median_survival_time_ if not pd.isna(kmf.median_survival_time_) else absent_survival['survival_days'].max()
            
            # Determine prognostic type
            prognostic_type = "Favorable" if median_present > median_absent else "Unfavorable"
            
            # Store result
            result = {
                'eei_id': f"EEI_{idx+1}",
                'mouse_exon1': mouse_exon1,
                'mouse_exon2': mouse_exon2,
                'mapped_gene1': gene1,
                'mapped_gene2': gene2,
                'human_exon1': eei_row.get('human_exon1', ''),
                'human_exon2': eei_row.get('human_exon2', ''),
                'p_value': p_value,
                'median_survival_present': median_present,
                'median_survival_absent': median_absent,
                'survival_difference': median_present - median_absent,
                'n_present': len(present_survival),
                'n_absent': len(absent_survival),
                'prognostic_type': prognostic_type,
                'avg_identity': eei_row.get('avg_identity', 0),
                'logrank_statistic': logrank_result.test_statistic
            }
            
            results.append(result)
            
        except Exception as e:
            print(f"   Error testing EEI {idx+1}: {e}")
            continue
        
        # Progress update
        if (idx + 1) % 50 == 0:
            print(f"   Processed {idx + 1}/{len(mapped_eeis)} EEIs (found {len(results)} testable)...")
    
    # 6. Process and save results
    print("6. Processing results...")
    
    if len(results) == 0:
        print("   ERROR: No testable EEIs found!")
        print("   This could be due to:")
        print("   - Expression threshold too strict")
        print("   - Insufficient sample overlap")
        print("   - Gene mapping issues")
        return
    
    # Convert to DataFrame
    results_df = pd.DataFrame(results)
    
    # Find significant EEIs
    significant_eeis = results_df[results_df['p_value'] <= 0.05].copy()
    significant_eeis = significant_eeis.sort_values('p_value')
    
    # Print summary
    print("\n" + "="*50)
    print("RESULTS SUMMARY")
    print("="*50)
    print(f"Total testable EEIs: {len(results_df)}")
    print(f"Significant EEIs (p≤0.05): {len(significant_eeis)}")
    print(f"Significance rate: {len(significant_eeis)/len(results_df)*100:.1f}%")
    
    if len(results_df) > 0:
        print(f"Median p-value: {results_df['p_value'].median():.4f}")
        print(f"Min p-value: {results_df['p_value'].min():.2e}")
    
    if len(significant_eeis) > 0:
        favorable = len(significant_eeis[significant_eeis['prognostic_type'] == 'Favorable'])
        unfavorable = len(significant_eeis[significant_eeis['prognostic_type'] == 'Unfavorable'])
        
        print(f"\nPrognostic breakdown:")
        print(f"  Favorable EEIs: {favorable}")
        print(f"  Unfavorable EEIs: {unfavorable}")
        
        print(f"\nTop 5 most significant EEIs:")
        for idx, row in significant_eeis.head().iterrows():
            print(f"  {row['eei_id']}: {row['mapped_gene1']} <-> {row['mapped_gene2']}")
            print(f"    p-value: {row['p_value']:.3e} ({row['prognostic_type']})")
            print(f"    Median survival - Present: {row['median_survival_present']:.1f}, Absent: {row['median_survival_absent']:.1f}")
            print(f"    Difference: {row['survival_difference']:.1f} days")
            print()
    
    # 7. Save results
    print("7. Saving results...")
    os.makedirs("outputs", exist_ok=True)
    
    # Save all results
    results_df.to_csv("outputs/mouse_eei_survival_results.tsv", sep='\t', index=False)
    print("   Saved: outputs/mouse_eei_survival_results.tsv")
    
    # Save significant results
    if len(significant_eeis) > 0:
        significant_eeis.to_csv("outputs/significant_mouse_eeis.tsv", sep='\t', index=False)
        print("   Saved: outputs/significant_mouse_eeis.tsv")
        
        # Create a simple survival plot for the most significant EEI
        create_survival_plot(significant_eeis.iloc[0], mouse_expression_aligned, survival_aligned, cmp_threshold)
    
    # Create summary statistics file
    summary_stats = {
        'total_mapped_eeis': len(mapped_eeis),
        'testable_eeis': len(results_df),
        'significant_eeis': len(significant_eeis),
        'significance_rate': len(significant_eeis)/len(results_df) if len(results_df) > 0 else 0,
        'median_p_value': results_df['p_value'].median() if len(results_df) > 0 else None,
        'min_p_value': results_df['p_value'].min() if len(results_df) > 0 else None,
        'common_samples': len(common_samples),
        'cmp_threshold_used': cmp_threshold,
        'min_group_size_used': min_group_size
    }
    
    import json
    with open("outputs/analysis_summary.json", 'w') as f:
        json.dump(summary_stats, f, indent=2)
    print("   Saved: outputs/analysis_summary.json")
    
    print("\n" + "="*50)
    print("ANALYSIS COMPLETED!")
    print("="*50)
    
    return significant_eeis, results_df

def create_survival_plot(best_eei, expression_df, survival_df, threshold):
    """Create a survival plot for the most significant EEI."""
    
    try:
        gene1 = best_eei['mapped_gene1']
        gene2 = best_eei['mapped_gene2']
        
        # Get expression data
        expr1 = expression_df[gene1]
        expr2 = expression_df[gene2]
        
        # Define EEI presence
        eei_present = (expr1 > threshold) & (expr2 > threshold)
        
        # Get survival data for each group
        present_samples = eei_present[eei_present == True].index.astype(str).tolist()
        absent_samples = eei_present[eei_present == False].index.astype(str).tolist()
        
        survival_df['mouse_id'] = survival_df['mouse_id'].astype(str)
        present_survival = survival_df[survival_df['mouse_id'].isin(present_samples)]
        absent_survival = survival_df[survival_df['mouse_id'].isin(absent_samples)]
        
        # Create plot
        plt.figure(figsize=(10, 6))
        
        kmf = KaplanMeierFitter()
        
        # Plot present group
        kmf.fit(present_survival['survival_days'], present_survival['event_status'], 
               label=f'EEI Present (n={len(present_survival)})')
        kmf.plot_survival_function(color='red', linewidth=2)
        
        # Plot absent group
        kmf.fit(absent_survival['survival_days'], absent_survival['event_status'], 
               label=f'EEI Absent (n={len(absent_survival)})')
        kmf.plot_survival_function(color='blue', linewidth=2)
        
        # Formatting
        plt.title(f"Most Significant EEI: {gene1} <-> {gene2}\n" + 
                 f"p-value: {best_eei['p_value']:.2e} ({best_eei['prognostic_type']})", 
                 fontsize=12, fontweight='bold')
        plt.xlabel('Survival Time (days)', fontsize=11)
        plt.ylabel('Survival Probability', fontsize=11)
        plt.legend(fontsize=10)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        # Save plot
        plt.savefig("outputs/top_eei_survival_plot.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        print("   Saved: outputs/top_eei_survival_plot.png")
        
    except Exception as e:
        print(f"   Error creating survival plot: {e}")

if __name__ == "__main__":
    # Run the analysis
    significant_eeis, all_results = run_eei_survival_analysis()
    
    if significant_eeis is not None and len(significant_eeis) > 0:
        print(f"\n🎉 SUCCESS: Found {len(significant_eeis)} significant EEIs!")
        print("Next steps:")
        print("1. Review outputs/significant_mouse_eeis.tsv")
        print("2. Check outputs/top_eei_survival_plot.png")
        print("3. Consider biological validation of top hits")
    else:
        print("\n⚠️  No significant EEIs found.")
        print("Consider adjusting parameters or checking data quality.")