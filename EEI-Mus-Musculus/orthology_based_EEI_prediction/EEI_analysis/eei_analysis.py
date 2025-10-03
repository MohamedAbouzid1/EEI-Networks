import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import networkx as nx
from sklearn.metrics import precision_recall_curve, auc, roc_curve, roc_auc_score
import os

class EEIAnalyzer:
    """Analyzer for Exon-Exon Interaction prediction results"""
    
    def __init__(self, 
                 human_eei_file=None, 
                 mouse_eei_file=None, 
                 egio_file=None, 
                 predicted_eei_file=None):
        """Initialize the analyzer with file paths"""
        self.human_eei_file = human_eei_file
        self.mouse_eei_file = mouse_eei_file
        self.egio_file = egio_file
        self.predicted_eei_file = predicted_eei_file
        
        # Data frames to be loaded
        self.human_eei_df = None
        self.mouse_eei_df = None
        self.egio_df = None
        self.predicted_eei_df = None
        
        # Statistics
        self.stats = {}
        
    def load_data(self):
        """Load all datasets"""
        print("Loading data files...")
        
        if self.human_eei_file:
            print(f"Loading human EEI network from {self.human_eei_file}")
            self.human_eei_df = pd.read_csv(self.human_eei_file, sep='\t')
            print(f"Loaded {len(self.human_eei_df)} human EEIs")
            
        if self.mouse_eei_file:
            print(f"Loading mouse EEI network from {self.mouse_eei_file}")
            self.mouse_eei_df = pd.read_csv(self.mouse_eei_file, sep='\t')
            print(f"Loaded {len(self.mouse_eei_df)} mouse EEIs")
            
        if self.egio_file:
            print(f"Loading EGIO orthology mappings from {self.egio_file}")
            self.egio_df = pd.read_csv(self.egio_file, sep='\t')
            print(f"Loaded {len(self.egio_df)} orthology mappings")
            
        if self.predicted_eei_file:
            print(f"Loading predicted EEIs from {self.predicted_eei_file}")
            self.predicted_eei_df = pd.read_csv(self.predicted_eei_file, sep='\t')
            print(f"Loaded {len(self.predicted_eei_df)} predicted EEIs")
    
    def compute_basic_stats(self):
        """Compute basic statistics for each dataset"""
        print("Computing basic statistics...")
        
        # Human EEI network stats
        if self.human_eei_df is not None:
            human_exons = set(self.human_eei_df['exon1'].unique()).union(
                set(self.human_eei_df['exon2'].unique()))
            
            self.stats['human_eei_count'] = len(self.human_eei_df)
            self.stats['human_exon_count'] = len(human_exons)
            
            # Handle different protein column naming conventions (CONTACT vs PISA)
            protein_cols = []
            if 'protein1' in self.human_eei_df.columns:
                protein_cols.append('protein1')
            if 'protein1.1' in self.human_eei_df.columns:
                protein_cols.append('protein1.1')
            elif 'protein2' in self.human_eei_df.columns:
                protein_cols.append('protein2')
            
            if protein_cols:
                unique_proteins = set()
                for col in protein_cols:
                    unique_proteins.update(self.human_eei_df[col].unique())
                self.stats['human_proteins'] = len(unique_proteins)
            else:
                self.stats['human_proteins'] = "N/A (protein columns not found)"
            
        # Mouse EEI network stats
        if self.mouse_eei_df is not None:
            mouse_exons = set(self.mouse_eei_df['exon1'].unique()).union(
                set(self.mouse_eei_df['exon2'].unique()))
            
            self.stats['mouse_eei_count'] = len(self.mouse_eei_df)
            self.stats['mouse_exon_count'] = len(mouse_exons)
            
            # Handle different protein column naming conventions (CONTACT vs PISA)
            protein_cols = []
            if 'protein1' in self.mouse_eei_df.columns:
                protein_cols.append('protein1')
            if 'protein1.1' in self.mouse_eei_df.columns:
                protein_cols.append('protein1.1')
            elif 'protein2' in self.mouse_eei_df.columns:
                protein_cols.append('protein2')
            
            if protein_cols:
                unique_proteins = set()
                for col in protein_cols:
                    unique_proteins.update(self.mouse_eei_df[col].unique())
                self.stats['mouse_proteins'] = len(unique_proteins)
            else:
                self.stats['mouse_proteins'] = "N/A (protein columns not found)"
            
        # EGIO orthology stats
        if self.egio_df is not None:
            self.stats['orthology_count'] = len(self.egio_df)
            self.stats['human_exons_with_orthologs'] = len(self.egio_df['hsaPos'].unique())
            self.stats['mouse_exons_with_orthologs'] = len(self.egio_df['musPos'].unique())
            
            # Count orthology mapping types
            type_counts = Counter(self.egio_df['Type'])
            self.stats['orthology_types'] = dict(type_counts)
            
            # Identity score distribution
            self.egio_df['Iden'] = pd.to_numeric(self.egio_df['Iden'], errors='coerce')
            self.stats['identity_mean'] = self.egio_df['Iden'].mean()
            self.stats['identity_median'] = self.egio_df['Iden'].median()
            self.stats['identity_std'] = self.egio_df['Iden'].std()
            
        # Predicted EEI stats
        if self.predicted_eei_df is not None:
            self.stats['predicted_eei_count'] = len(self.predicted_eei_df)
            self.stats['predicted_exons'] = len(set(self.predicted_eei_df['exon1'].unique()).union(
                set(self.predicted_eei_df['exon2'].unique())))
            
            # Confidence score stats
            if 'confidence' in self.predicted_eei_df.columns:
                self.stats['confidence_mean'] = self.predicted_eei_df['confidence'].mean()
                self.stats['confidence_median'] = self.predicted_eei_df['confidence'].median()
                self.stats['confidence_std'] = self.predicted_eei_df['confidence'].std()
        
        print("Basic statistics computed.")
        return self.stats
    
    def analyze_prediction_overlap(self):
        """Analyze overlap between predicted EEIs and existing human EEIs"""
        if self.human_eei_df is None or self.predicted_eei_df is None:
            print("Both human EEI and predicted EEI datasets must be loaded")
            return
            
        print("Analyzing prediction overlap with existing human EEIs...")
        
        # Create sets of EEI pairs for comparison
        human_eei_pairs = set()
        for _, row in self.human_eei_df.iterrows():
            pair = (row['exon1'], row['exon2'])
            human_eei_pairs.add(pair)
            # Add reverse pair too for undirected network
            human_eei_pairs.add((row['exon2'], row['exon1']))
            
        predicted_eei_pairs = set()
        for _, row in self.predicted_eei_df.iterrows():
            pair = (row['exon1'], row['exon2'])
            predicted_eei_pairs.add(pair)
            # Add reverse pair too for undirected network
            predicted_eei_pairs.add((row['exon2'], row['exon1']))
            
        # Calculate overlap
        overlapping_pairs = human_eei_pairs.intersection(predicted_eei_pairs)
        
        self.stats['overlapping_eei_count'] = len(overlapping_pairs) // 2  # Divide by 2 because we counted both directions
        self.stats['novel_eei_count'] = len(predicted_eei_pairs) // 2 - self.stats['overlapping_eei_count']
        
        if len(predicted_eei_pairs) > 0:
            self.stats['overlap_percentage'] = (self.stats['overlapping_eei_count'] / (len(predicted_eei_pairs) // 2)) * 100
        else:
            self.stats['overlap_percentage'] = 0
            
        print(f"Found {self.stats['overlapping_eei_count']} overlapping EEIs ({self.stats['overlap_percentage']:.2f}%)")
        print(f"Identified {self.stats['novel_eei_count']} novel predicted EEIs")
        
        return self.stats
    
    def analyze_conservation_patterns(self):
        """Analyze conservation patterns between mouse and human EEIs"""
        if self.human_eei_df is None or self.mouse_eei_df is None or self.egio_df is None:
            print("Human EEI, mouse EEI, and EGIO datasets must be loaded")
            return
            
        print("Analyzing conservation patterns...")
        
        # Create mouse to human exon mapping dictionary
        mouse_to_human = {}
        for _, row in self.egio_df.iterrows():
            if not pd.isna(row['musPos']) and not pd.isna(row['hsaPos']):
                mouse_to_human[row['musPos']] = row['hsaPos']
        
        # Count how many mouse EEIs have both exons mapped to human
        mappable_eeis = 0
        conserved_eeis = 0
        
        # Create sets of human EEI pairs
        human_eei_pairs = set()
        for _, row in self.human_eei_df.iterrows():
            pair = (row['exon1'], row['exon2'])
            human_eei_pairs.add(pair)
            # Add reverse pair too for undirected network
            human_eei_pairs.add((row['exon2'], row['exon1']))
        
        for _, row in self.mouse_eei_df.iterrows():
            mouse_exon1 = row['exon1']
            mouse_exon2 = row['exon2']
            
            if mouse_exon1 in mouse_to_human and mouse_exon2 in mouse_to_human:
                mappable_eeis += 1
                
                human_exon1 = mouse_to_human[mouse_exon1]
                human_exon2 = mouse_to_human[mouse_exon2]
                
                if (human_exon1, human_exon2) in human_eei_pairs:
                    conserved_eeis += 1
        
        self.stats['mappable_mouse_eeis'] = mappable_eeis
        self.stats['conserved_eeis'] = conserved_eeis
        
        if mappable_eeis > 0:
            self.stats['conservation_rate'] = (conserved_eeis / mappable_eeis) * 100
        else:
            self.stats['conservation_rate'] = 0
            
        print(f"Found {mappable_eeis} mouse EEIs with orthologous exons in human")
        print(f"Of these, {conserved_eeis} ({self.stats['conservation_rate']:.2f}%) are conserved in human")
        
        return self.stats
    
    def analyze_prediction_quality(self):
        """Analyze quality of predictions based on confidence scores"""
        if self.predicted_eei_df is None or 'confidence' not in self.predicted_eei_df.columns:
            print("Predicted EEIs with confidence scores must be loaded")
            return
            
        print("Analyzing prediction quality...")
        
        # Bin predictions by confidence score
        confidence_bins = [0.6, 0.7, 0.8, 0.9, 1.0]
        bin_counts = []
        
        for i in range(len(confidence_bins)):
            if i == 0:
                count = len(self.predicted_eei_df[self.predicted_eei_df['confidence'] < confidence_bins[i]])
            else:
                count = len(self.predicted_eei_df[(self.predicted_eei_df['confidence'] >= confidence_bins[i-1]) & 
                                                 (self.predicted_eei_df['confidence'] < confidence_bins[i])])
            bin_counts.append(count)
            
        # Add last bin
        bin_counts.append(len(self.predicted_eei_df[self.predicted_eei_df['confidence'] >= confidence_bins[-1]]))
        
        # Create bin labels
        bin_labels = [f'<{confidence_bins[0]}']
        for i in range(len(confidence_bins)-1):
            bin_labels.append(f'{confidence_bins[i]}-{confidence_bins[i+1]}')
        bin_labels.append(f'≥{confidence_bins[-1]}')
        
        self.stats['confidence_bins'] = bin_labels
        self.stats['bin_counts'] = bin_counts
        
        print("Confidence score distribution:")
        for label, count in zip(bin_labels, bin_counts):
            print(f"  {label}: {count} predictions")
        
        return self.stats
    
    def create_visualizations(self, output_dir='figures'):
        """Create visualizations of analysis results"""
        print(f"Creating visualizations in {output_dir}...")
        
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Set style
        plt.style.use('seaborn-v0_8-darkgrid')
        sns.set_context("talk")
        
        # 1. Basic count comparison
        if all(k in self.stats for k in ['human_eei_count', 'mouse_eei_count', 'predicted_eei_count']):
            plt.figure(figsize=(10, 6))
            counts = [self.stats['human_eei_count'], self.stats['mouse_eei_count'], self.stats['predicted_eei_count']]
            labels = ['Human EEIs', 'Mouse EEIs', 'Predicted EEIs']
            
            plt.bar(labels, counts)
            plt.title('EEI Network Sizes')
            plt.ylabel('Count')
            plt.tight_layout()
            plt.savefig(f"{output_dir}/network_sizes.png", dpi=300)
            plt.close()
        
        # 2. Orthology type distribution
        if 'orthology_types' in self.stats:
            plt.figure(figsize=(10, 6))
            types = list(self.stats['orthology_types'].keys())
            type_counts = list(self.stats['orthology_types'].values())
            
            plt.bar(types, type_counts)
            plt.title('Orthology Mapping Types')
            plt.ylabel('Count')
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.savefig(f"{output_dir}/orthology_types.png", dpi=300)
            plt.close()
        
        # 3. Identity score distribution
        if self.egio_df is not None and 'Iden' in self.egio_df.columns:
            plt.figure(figsize=(10, 6))
            sns.histplot(self.egio_df['Iden'].dropna(), bins=20, kde=True)
            plt.title('Identity Score Distribution')
            plt.xlabel('Identity Score')
            plt.ylabel('Count')
            plt.axvline(x=0.8, color='red', linestyle='--', label='Threshold (0.8)')
            plt.legend()
            plt.tight_layout()
            plt.savefig(f"{output_dir}/identity_distribution.png", dpi=300)
            plt.close()
            
        # 4. Confidence score distribution of predictions
        if self.predicted_eei_df is not None and 'confidence' in self.predicted_eei_df.columns:
            plt.figure(figsize=(10, 6))
            sns.histplot(self.predicted_eei_df['confidence'].dropna(), bins=20, kde=True)
            plt.title('Prediction Confidence Distribution')
            plt.xlabel('Confidence Score')
            plt.ylabel('Count')
            plt.tight_layout()
            plt.savefig(f"{output_dir}/confidence_distribution.png", dpi=300)
            plt.close()
            
        # 5. Confidence bin distribution
        if all(k in self.stats for k in ['confidence_bins', 'bin_counts']):
            plt.figure(figsize=(12, 6))
            plt.bar(self.stats['confidence_bins'], self.stats['bin_counts'])
            plt.title('Predictions by Confidence Score Range')
            plt.xlabel('Confidence Score Range')
            plt.ylabel('Number of Predictions')
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.savefig(f"{output_dir}/confidence_bins.png", dpi=300)
            plt.close()
        
        # 6. Conservation analysis
        if all(k in self.stats for k in ['mappable_mouse_eeis', 'conserved_eeis']):
            if self.stats['mappable_mouse_eeis'] > 0 or self.stats['conserved_eeis'] > 0:
                plt.figure(figsize=(8, 8))
                labels = ['Non-conserved EEIs', 'Conserved EEIs']
                sizes = [self.stats['mappable_mouse_eeis'] - self.stats['conserved_eeis'], 
                        self.stats['conserved_eeis']]
                
                plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90, colors=['#FF9999', '#66B2FF'])
                plt.axis('equal')
                plt.title('Conservation of Mappable Mouse EEIs in Human')
                plt.tight_layout()
                plt.savefig(f"{output_dir}/conservation_pie.png", dpi=300)
                plt.close()
            else:
                print("Skipping conservation pie chart - no valid data")
                
        # 7. Prediction overlap analysis
        if all(k in self.stats for k in ['novel_eei_count', 'overlapping_eei_count']):
        # Check if we have valid data for the pie chart
            if self.stats['novel_eei_count'] > 0 or self.stats['overlapping_eei_count'] > 0:
                plt.figure(figsize=(8, 8))
                labels = ['Novel Predictions', 'Already Known']
                sizes = [self.stats['novel_eei_count'], self.stats['overlapping_eei_count']]
                
                plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90, colors=['#66B2FF', '#FF9999'])
                plt.axis('equal')
                plt.title('Breakdown of Predicted EEIs')
                plt.tight_layout()
                plt.savefig(f"{output_dir}/prediction_overlap_pie.png", dpi=300)
                plt.close()
            else:
                print("Skipping prediction overlap pie chart - no valid data")
    
        
    def improved_analyze_prediction_quality(self):

        """Enhanced function to analyze prediction quality with confidence score breakdown by type"""
        
        if self.predicted_eei_df is None:
            print("Predicted EEI dataset must be loaded")
            return
                
        print("Analyzing prediction quality with relationship type breakdown...")
        
        # Get relationship types if available
        has_types = all(col in self.predicted_eei_df.columns for col in ['type1', 'type2'])
        
        if has_types:
            # Create relationship pattern
            self.predicted_eei_df['relationship_pattern'] = self.predicted_eei_df.apply(
                lambda x: f"{x['type1']}-{x['type2']}", axis=1
            )
            
            # Analyze confidence by relationship pattern
            pattern_stats = {}
            for pattern, group in self.predicted_eei_df.groupby('relationship_pattern'):
                if 'confidence' in group.columns:
                    pattern_stats[pattern] = {
                        'count': len(group),
                        'confidence_mean': group['confidence'].mean(),
                        'confidence_median': group['confidence'].median(),
                        'confidence_std': group['confidence'].std(),
                        'confidence_min': group['confidence'].min(),
                        'confidence_max': group['confidence'].max()
                    }
            
            self.stats['relationship_pattern_stats'] = pattern_stats
            
            # Create visualization of confidence distributions by pattern
            self.create_pattern_confidence_visualization()
        
        # Original confidence score analysis (keep this from the original function)
        if 'confidence' in self.predicted_eei_df.columns:
            # Bin predictions by confidence score
            confidence_bins = [0.6, 0.7, 0.8, 0.9, 1.0]
            bin_counts = []
            
            for i in range(len(confidence_bins)):
                if i == 0:
                    count = len(self.predicted_eei_df[self.predicted_eei_df['confidence'] < confidence_bins[i]])
                else:
                    count = len(self.predicted_eei_df[(self.predicted_eei_df['confidence'] >= confidence_bins[i-1]) & 
                                                    (self.predicted_eei_df['confidence'] < confidence_bins[i])])
                bin_counts.append(count)
                
            # Add last bin
            bin_counts.append(len(self.predicted_eei_df[self.predicted_eei_df['confidence'] >= confidence_bins[-1]]))
            
            # Create bin labels
            bin_labels = [f'<{confidence_bins[0]}']
            for i in range(len(confidence_bins)-1):
                bin_labels.append(f'{confidence_bins[i]}-{confidence_bins[i+1]}')
            bin_labels.append(f'≥{confidence_bins[-1]}')
            
            self.stats['confidence_bins'] = bin_labels
            self.stats['bin_counts'] = bin_counts
            
            print("Confidence score distribution:")
            for label, count in zip(bin_labels, bin_counts):
                print(f"  {label}: {count} predictions")
        
        return self.stats

    def create_pattern_confidence_visualization(self, output_dir='figures'):
        """Create visualizations of confidence scores by relationship pattern"""
        if 'relationship_pattern_stats' not in self.stats:
            return
            
        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
            
        patterns = list(self.stats['relationship_pattern_stats'].keys())
        counts = [self.stats['relationship_pattern_stats'][p]['count'] for p in patterns]
        
        # Check if confidence scores are available
        has_confidence = 'confidence_mean' in self.stats['relationship_pattern_stats'][patterns[0]]
        
        if has_confidence:
            means = [self.stats['relationship_pattern_stats'][p]['confidence_mean'] for p in patterns]
        
            # Sort by count
            sorted_indices = np.argsort(counts)[::-1]
            sorted_patterns = [patterns[i] for i in sorted_indices]
            sorted_counts = [counts[i] for i in sorted_indices]
            sorted_means = [means[i] for i in sorted_indices]
            
            # Top 10 patterns
            if len(sorted_patterns) > 10:
                top_patterns = sorted_patterns[:10]
                top_counts = sorted_counts[:10]
                top_means = sorted_means[:10]
            else:
                top_patterns = sorted_patterns
                top_counts = sorted_counts
                top_means = sorted_means
            
            # Create figure with two subplots
            plt.figure(figsize=(14, 10))
            
            # Count by pattern
            ax1 = plt.subplot(2, 1, 1)
            bars = ax1.bar(top_patterns, top_counts)
            ax1.set_title('Prediction Counts by Relationship Pattern')
            ax1.set_xlabel('Relationship Pattern')
            ax1.set_ylabel('Count')
            ax1.set_xticklabels(top_patterns, rotation=45, ha='right')
            
            # Mean confidence by pattern
            ax2 = plt.subplot(2, 1, 2)
            bars = ax2.bar(top_patterns, top_means)
            ax2.set_title('Mean Confidence by Relationship Pattern')
            ax2.set_xlabel('Relationship Pattern')
            ax2.set_ylabel('Mean Confidence')
            ax2.set_xticklabels(top_patterns, rotation=45, ha='right')
            ax2.set_ylim(0, 1)
            
            plt.tight_layout()
            plt.savefig(f"{output_dir}/pattern_confidence.png", dpi=300)
            plt.close()
        else:
            # Just show pattern counts if no confidence scores
            # Sort by count
            sorted_indices = np.argsort(counts)[::-1]
            sorted_patterns = [patterns[i] for i in sorted_indices]
            sorted_counts = [counts[i] for i in sorted_indices]
            
            # Top 10 patterns
            if len(sorted_patterns) > 10:
                top_patterns = sorted_patterns[:10]
                top_counts = sorted_counts[:10]
            else:
                top_patterns = sorted_patterns
                top_counts = sorted_counts
                
            plt.figure(figsize=(12, 6))
            bars = plt.bar(top_patterns, top_counts)
            plt.title('Prediction Counts by Relationship Pattern')
            plt.xlabel('Relationship Pattern')
            plt.ylabel('Count')
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            plt.savefig(f"{output_dir}/pattern_counts.png", dpi=300)
            plt.close()    

    def run_full_analysis(self):
        """Run the complete analysis pipeline"""
        self.load_data()
        self.compute_basic_stats()
        
        if self.human_eei_df is not None and self.predicted_eei_df is not None:
            self.analyze_prediction_overlap()
            
        if self.human_eei_df is not None and self.mouse_eei_df is not None and self.egio_df is not None:
            self.analyze_conservation_patterns()
            
        if self.predicted_eei_df is not None:
            # Use the improved version instead of the original
            self.improved_analyze_prediction_quality()
            
        self.create_visualizations()
        
        # Print summary statistics
        print("\nSummary Statistics:")
        for key, value in self.stats.items():
            if not isinstance(value, dict):
                print(f"  {key}: {value}")

if __name__ == "__main__":
    # Example usage
    analyzer = EEIAnalyzer(
        human_eei_file="path/to/human_eei_network.tsv",
        mouse_eei_file="path/to/mouse_eei_network.tsv",
        egio_file="path/to/egio_output.tsv",
        predicted_eei_file="path/to/predicted_eeis.tsv"
    )
    
    analyzer.run_full_analysis()