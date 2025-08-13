import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from scipy.stats import chi2_contingency, spearmanr
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

class CrossSpeciesEEIAnalysis:
    """
    Improved pipeline for cross-species EEI analysis
    """
    
    def __init__(self):
        self.expression_df = None
        self.survival_df = None
        self.mouse_eei_df = None
        self.human_crpes = None
        self.ortholog_mapping = None
        self.results = {}
        
    def load_data(self, expression_df, survival_df, mouse_eei_df):
        """Load all required datasets"""
        self.expression_df = expression_df
        self.survival_df = survival_df
        self.mouse_eei_df = mouse_eei_df
        print(f"✅ Loaded expression data: {expression_df.shape}")
        print(f"✅ Loaded survival data: {len(survival_df)} samples")
        print(f"✅ Loaded mouse EEI network: {len(mouse_eei_df)} interactions")
        
    def map_samples_improved(self):
        """
        Improved sample mapping between expression and metadata
        """
        print("\n🔄 Mapping expression samples to metadata...")
        
        # Extract tissue and condition info from expression columns
        sample_mapping = {}
        tissue_patterns = []
        
        for col in self.expression_df.columns:
            tissue_condition, sample_id = col
            # Extract tissue type and cachexia status from name
            tissue = tissue_condition.split('_')[0]  # gastro, gWAT, Tumor
            cachexia_status = 'Positive' if '+CACS' in tissue_condition else 'Negative'
            rep_num = tissue_condition.split('_rep')[-1] if '_rep' in tissue_condition else '1'
            
            tissue_patterns.append({
                'expression_column': col,
                'tissue': tissue,
                'cachexia_status': cachexia_status,
                'replicate': rep_num,
                'full_name': tissue_condition
            })
        
        tissue_df = pd.DataFrame(tissue_patterns)
        
        # Match with survival data based on characteristics
        matched_samples = []
        for idx, row in self.survival_df.iterrows():
            characteristics = row['characteristics'] if 'characteristics' in row else []
            sample_cachexia = row['cachexia_status']
            
            # Find best match in expression data
            potential_matches = tissue_df[tissue_df['cachexia_status'] == sample_cachexia]
            if len(potential_matches) > 0:
                # Take first available match (could be improved with better logic)
                match_idx = len(matched_samples) % len(potential_matches)
                matched_sample = potential_matches.iloc[match_idx]
                matched_samples.append({
                    'gsm_id': row['sample_id'],
                    'expression_column': matched_sample['expression_column'],
                    'tissue': matched_sample['tissue'],
                    'cachexia_status': sample_cachexia,
                    'weight_loss_percent': row['weight_loss_percent']
                })
        
        self.sample_mapping_df = pd.DataFrame(matched_samples)
        print(f"✅ Mapped {len(self.sample_mapping_df)} samples")
        return self.sample_mapping_df
    
    def define_eei_presence_improved(self, expression_threshold_percentile=50):
        """
        Improved EEI presence definition
        """
        print(f"\n🔄 Defining EEI presence (threshold: {expression_threshold_percentile}th percentile)...")
        
        # Create sample mapping
        sample_mapping = self.map_samples_improved()
        
        # Filter to mapped samples only
        mapped_columns = [col for col in self.expression_df.columns 
                         if col in sample_mapping['expression_column'].values]
        
        eei_presence_data = []
        tested_eeis = 0
        
        for idx, eei in self.mouse_eei_df.iterrows():
            exon1, exon2 = eei['exon1'], eei['exon2']
            
            # Check if both exons exist in expression data
            if exon1 in self.expression_df.index and exon2 in self.expression_df.index:
                tested_eeis += 1
                
                # Get expression for both exons (mapped samples only)
                expr1 = self.expression_df.loc[exon1, mapped_columns]
                expr2 = self.expression_df.loc[exon2, mapped_columns]
                
                # Define thresholds
                threshold1 = np.percentile(expr1, expression_threshold_percentile)
                threshold2 = np.percentile(expr2, expression_threshold_percentile)
                
                # EEI present if both exons above threshold
                eei_present = (expr1 > threshold1) & (expr2 > threshold2)
                
                # Store results for each sample
                for col in mapped_columns:
                    sample_info = sample_mapping[
                        sample_mapping['expression_column'] == col
                    ].iloc[0]
                    
                    eei_presence_data.append({
                        'eei_id': f"{exon1}_{exon2}",
                        'exon1': exon1,
                        'exon2': exon2,
                        'sample_id': sample_info['gsm_id'],
                        'expression_column': col,
                        'tissue': sample_info['tissue'],
                        'cachexia_status': sample_info['cachexia_status'],
                        'weight_loss_percent': sample_info['weight_loss_percent'],
                        'eei_present': eei_present[col],
                        'exon1_expr': expr1[col],
                        'exon2_expr': expr2[col]
                    })
        
        self.eei_presence_df = pd.DataFrame(eei_presence_data)
        print(f"✅ Tested {tested_eeis} EEIs across {len(mapped_columns)} samples")
        print(f"✅ Created EEI presence matrix: {len(self.eei_presence_df)} entries")
        
        return self.eei_presence_df
    
    def test_eei_associations(self, min_samples=10):
        """
        Test EEI associations with cachexia using multiple approaches
        """
        print(f"\n🔄 Testing EEI associations (min_samples: {min_samples})...")
        
        if self.eei_presence_df is None:
            self.define_eei_presence_improved()
        
        results = []
        
        # Group by EEI for testing
        for eei_id, eei_data in self.eei_presence_df.groupby('eei_id'):
            if len(eei_data) < min_samples:
                continue
                
            exon1 = eei_data.iloc[0]['exon1']
            exon2 = eei_data.iloc[0]['exon2']
            
            # Test 1: Association with binary cachexia status
            try:
                contingency = pd.crosstab(
                    eei_data['eei_present'], 
                    eei_data['cachexia_status']
                )
                
                if contingency.shape == (2, 2):
                    chi2, p_binary, dof, expected = chi2_contingency(contingency)
                    
                    # Calculate effect size (Cramer's V)
                    n = contingency.sum().sum()
                    cramers_v = np.sqrt(chi2 / (n * (min(contingency.shape) - 1)))
                else:
                    p_binary, cramers_v = np.nan, np.nan
            except:
                p_binary, cramers_v = np.nan, np.nan
            
            # Test 2: Correlation with continuous weight loss
            try:
                valid_data = eei_data.dropna(subset=['weight_loss_percent'])
                if len(valid_data) >= min_samples:
                    corr_coef, p_continuous = spearmanr(
                        valid_data['eei_present'].astype(int),
                        valid_data['weight_loss_percent']
                    )
                else:
                    corr_coef, p_continuous = np.nan, np.nan
            except:
                corr_coef, p_continuous = np.nan, np.nan
            
            # Test 3: Tissue-specific analysis
            tissue_pvals = {}
            for tissue in eei_data['tissue'].unique():
                tissue_data = eei_data[eei_data['tissue'] == tissue]
                if len(tissue_data) >= 5:  # Lower threshold for tissue-specific
                    try:
                        tissue_contingency = pd.crosstab(
                            tissue_data['eei_present'],
                            tissue_data['cachexia_status']
                        )
                        if tissue_contingency.shape == (2, 2):
                            _, p_tissue, _, _ = chi2_contingency(tissue_contingency)
                            tissue_pvals[tissue] = p_tissue
                    except:
                        pass
            
            results.append({
                'eei_id': eei_id,
                'exon1': exon1,
                'exon2': exon2,
                'n_samples': len(eei_data),
                'n_tissues': eei_data['tissue'].nunique(),
                'p_binary_cachexia': p_binary,
                'cramers_v': cramers_v,
                'p_continuous_weight': p_continuous,
                'correlation_weight': corr_coef,
                'tissue_pvals': tissue_pvals,
                'min_tissue_pval': min(tissue_pvals.values()) if tissue_pvals else np.nan,
                'eei_prevalence': eei_data['eei_present'].mean(),
                'cachexic_with_eei': len(eei_data[
                    (eei_data['eei_present']) & 
                    (eei_data['cachexia_status'] == 'Positive')
                ]),
                'control_with_eei': len(eei_data[
                    (eei_data['eei_present']) & 
                    (eei_data['cachexia_status'] == 'Negative')
                ])
            })
        
        self.results_df = pd.DataFrame(results)
        
        # Apply multiple testing correction
        for col in ['p_binary_cachexia', 'p_continuous_weight', 'min_tissue_pval']:
            if col in self.results_df.columns:
                valid_pvals = self.results_df[col].dropna()
                if len(valid_pvals) > 0:
                    _, corrected_p, _, _ = multipletests(valid_pvals, method='fdr_bh')
                    corrected_dict = dict(zip(valid_pvals.index, corrected_p))
                    self.results_df[f'{col}_corrected'] = self.results_df.index.map(
                        lambda i: corrected_dict.get(i, np.nan)
                    )
        
        print(f"✅ Tested {len(self.results_df)} EEIs")
        
        # Summary statistics
        for test_type in ['binary_cachexia', 'continuous_weight', 'min_tissue_pval']:
            col = f'p_{test_type}_corrected'
            if col in self.results_df.columns:
                n_sig = (self.results_df[col] < 0.05).sum()
                print(f"   - {test_type}: {n_sig} significant (FDR < 0.05)")
        
        return self.results_df
    
    def create_visualizations(self):
        """
        Create comprehensive visualizations
        """
        print("\n🎨 Creating visualizations...")
        
        fig = plt.figure(figsize=(20, 15))
        
        # 1. Sample overview
        plt.subplot(3, 4, 1)
        if hasattr(self, 'sample_mapping_df'):
            sns.countplot(data=self.sample_mapping_df, x='tissue', hue='cachexia_status')
            plt.title('Sample Distribution by Tissue & Cachexia')
            plt.xticks(rotation=45)
        
        # 2. Weight loss distribution
        plt.subplot(3, 4, 2)
        if hasattr(self, 'sample_mapping_df'):
            sns.boxplot(data=self.sample_mapping_df, x='cachexia_status', y='weight_loss_percent')
            plt.title('Weight Loss by Cachexia Status')
        
        # 3. EEI prevalence distribution
        plt.subplot(3, 4, 3)
        if hasattr(self, 'results_df'):
            plt.hist(self.results_df['eei_prevalence'], bins=20, alpha=0.7)
            plt.xlabel('EEI Prevalence')
            plt.ylabel('Count')
            plt.title('Distribution of EEI Prevalence')
        
        # 4. P-value distributions
        plt.subplot(3, 4, 4)
        if hasattr(self, 'results_df'):
            valid_pvals = self.results_df['p_binary_cachexia'].dropna()
            if len(valid_pvals) > 0:
                plt.hist(valid_pvals, bins=20, alpha=0.7)
                plt.axhline(y=len(valid_pvals)/20, color='red', linestyle='--', 
                           label='Expected under null')
                plt.xlabel('P-value')
                plt.ylabel('Count')
                plt.title('P-value Distribution (Binary Cachexia)')
                plt.legend()
        
        # 5. Volcano plot
        plt.subplot(3, 4, 5)
        if hasattr(self, 'results_df'):
            valid_data = self.results_df.dropna(subset=['p_binary_cachexia', 'cramers_v'])
            if len(valid_data) > 0:
                x = valid_data['cramers_v']
                y = -np.log10(valid_data['p_binary_cachexia'])
                plt.scatter(x, y, alpha=0.6)
                plt.axhline(y=-np.log10(0.05), color='red', linestyle='--', label='p=0.05')
                plt.xlabel("Cramer's V (Effect Size)")
                plt.ylabel('-log10(P-value)')
                plt.title('Volcano Plot: EEI-Cachexia Associations')
                plt.legend()
        
        # 6. Correlation plot
        plt.subplot(3, 4, 6)
        if hasattr(self, 'results_df'):
            valid_data = self.results_df.dropna(subset=['p_continuous_weight', 'correlation_weight'])
            if len(valid_data) > 0:
                x = valid_data['correlation_weight']
                y = -np.log10(valid_data['p_continuous_weight'])
                plt.scatter(x, y, alpha=0.6)
                plt.axhline(y=-np.log10(0.05), color='red', linestyle='--', label='p=0.05')
                plt.xlabel('Correlation with Weight Loss')
                plt.ylabel('-log10(P-value)')
                plt.title('EEI-Weight Loss Correlations')
                plt.legend()
        
        # 7. Sample size vs significance
        plt.subplot(3, 4, 7)
        if hasattr(self, 'results_df'):
            valid_data = self.results_df.dropna(subset=['p_binary_cachexia'])
            if len(valid_data) > 0:
                colors = ['red' if p < 0.05 else 'blue' for p in valid_data['p_binary_cachexia']]
                plt.scatter(valid_data['n_samples'], 
                           -np.log10(valid_data['p_binary_cachexia']),
                           c=colors, alpha=0.6)
                plt.xlabel('Number of Samples')
                plt.ylabel('-log10(P-value)')
                plt.title('Sample Size vs Significance')
        
        # 8. Tissue-specific results
        plt.subplot(3, 4, 8)
        if hasattr(self, 'results_df'):
            valid_data = self.results_df.dropna(subset=['min_tissue_pval'])
            if len(valid_data) > 0:
                plt.hist(-np.log10(valid_data['min_tissue_pval']), bins=15, alpha=0.7)
                plt.axvline(x=-np.log10(0.05), color='red', linestyle='--', label='p=0.05')
                plt.xlabel('-log10(Min Tissue P-value)')
                plt.ylabel('Count')
                plt.title('Tissue-Specific Associations')
                plt.legend()
        
        # 9-12: Expression heatmaps for top EEIs (if any significant)
        if hasattr(self, 'results_df') and hasattr(self, 'eei_presence_df'):
            top_eeis = self.results_df.nsmallest(4, 'p_binary_cachexia')['eei_id'].tolist()
            
            for i, eei_id in enumerate(top_eeis):
                plt.subplot(3, 4, 9 + i)
                eei_data = self.eei_presence_df[self.eei_presence_df['eei_id'] == eei_id]
                
                # Create pivot table for heatmap
                pivot_data = eei_data.pivot_table(
                    values='eei_present', 
                    index='sample_id', 
                    columns='cachexia_status',
                    aggfunc='first'
                )
                
                if not pivot_data.empty:
                    sns.heatmap(pivot_data.T, cmap='RdYlBu_r', cbar_kws={'label': 'EEI Present'})
                    plt.title(f'Top EEI {i+1}: {eei_id}')
                    plt.ylabel('Cachexia Status')
        
        plt.tight_layout()
        plt.savefig('eei_analysis_comprehensive.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        return fig
    
    def generate_report(self):
        """
        Generate comprehensive analysis report
        """
        if not hasattr(self, 'results_df'):
            print("❌ No results available. Run test_eei_associations() first.")
            return
            
        print("\n📊 CROSS-SPECIES EEI ANALYSIS REPORT")
        print("=" * 50)
        
        # Data summary
        print(f"\n📋 DATA SUMMARY:")
        print(f"   • Expression samples: {self.expression_df.shape[1]}")
        print(f"   • Metadata samples: {len(self.survival_df)}")
        print(f"   • Mapped samples: {len(self.sample_mapping_df) if hasattr(self, 'sample_mapping_df') else 'N/A'}")
        print(f"   • Total EEIs tested: {len(self.results_df)}")
        print(f"   • Mouse EEI network size: {len(self.mouse_eei_df)}")
        
        # Results summary
        print(f"\n🎯 RESULTS SUMMARY:")
        for test_type in ['binary_cachexia', 'continuous_weight', 'min_tissue_pval']:
            col = f'p_{test_type}_corrected'
            if col in self.results_df.columns:
                n_total = self.results_df[col].notna().sum()
                n_sig_raw = (self.results_df[f'p_{test_type}'] < 0.05).sum()
                n_sig_corr = (self.results_df[col] < 0.05).sum()
                print(f"   • {test_type.replace('_', ' ').title()}:")
                print(f"     - Tested: {n_total}")
                print(f"     - Significant (raw p<0.05): {n_sig_raw}")
                print(f"     - Significant (FDR<0.05): {n_sig_corr}")
        
        # Top results
        if len(self.results_df) > 0:
            print(f"\n🏆 TOP 10 EEIs (by binary cachexia p-value):")
            top_results = self.results_df.nsmallest(10, 'p_binary_cachexia')
            for idx, row in top_results.iterrows():
                print(f"   {idx+1:2d}. {row['eei_id']}")
                print(f"       P-value: {row['p_binary_cachexia']:.2e}")
                print(f"       Effect size: {row['cramers_v']:.3f}")
                print(f"       Samples: {row['n_samples']}")
        
        # Recommendations
        print(f"\n💡 RECOMMENDATIONS:")
        if (self.results_df['p_binary_cachexia_corrected'] < 0.05).sum() == 0:
            print("   • No significant associations found after multiple testing correction")
            print("   • Consider: larger sample size, different thresholds, or pathway-level analysis")
            print("   • Focus on human-to-mouse ortholog mapping for targeted analysis")
        else:
            print("   • Significant associations detected - proceed with validation")
            print("   • Consider pathway enrichment analysis")
            print("   • Validate in additional mouse cancer datasets")
        
        return self.results_df

# Usage example:
def run_improved_analysis(expression_df, survival_df, mouse_eei_df):
    """
    Run the improved cross-species EEI analysis pipeline
    """
    print("🚀 Starting Cross-Species EEI Analysis Pipeline...")
    
    # Initialize analysis
    analyzer = CrossSpeciesEEIAnalysis()
    analyzer.load_data(expression_df, survival_df, mouse_eei_df)
    
    # Run analysis steps
    sample_mapping = analyzer.map_samples_improved()
    eei_presence = analyzer.define_eei_presence_improved()
    results = analyzer.test_eei_associations()
    
    # Create visualizations
    fig = analyzer.create_visualizations()
    
    # Generate report
    report = analyzer.generate_report()
    
    print("\n✅ Analysis complete!")
    
    return analyzer, results

# Example usage:
# analyzer, results = run_improved_analysis(expression_df, survival_df, mouse_eei_df)