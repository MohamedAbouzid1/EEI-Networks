import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import os
import networkx as nx

def get_protein_pair(row):
    """Extract protein pair from a row, handling different data formats"""
    protein1 = None
    protein2 = None
    
    if 'protein1' in row and not pd.isna(row['protein1']):
        protein1 = row['protein1']
    
    # Try 'protein1.1' (CONTACT format) and 'protein2' (PISA format)
    if 'protein1.1' in row and not pd.isna(row['protein1.1']):
        protein2 = row['protein1.1']
    elif 'protein2' in row and not pd.isna(row['protein2']):
        protein2 = row['protein2']
        
    if protein1 and protein2:
        return tuple(sorted([protein1, protein2]))
    else:
        return None
    
def load_data(human_eei_file, mouse_eei_file, egio_file, predicted_file=None):
    """
    Load all necessary data files
    
    Args:
        human_eei_file: Path to human EEI network file
        mouse_eei_file: Path to mouse EEI network file
        egio_file: Path to EGIO orthology mapping file
        predicted_file: Optional path to predicted EEIs file
        
    Returns:
        Tuple of DataFrames
    """
    human_eei_df = pd.read_csv(human_eei_file, sep='\t')
    mouse_eei_df = pd.read_csv(mouse_eei_file, sep='\t')
    egio_df = pd.read_csv(egio_file, sep='\t')
    
    # Load predicted EEIs if provided
    predicted_df = None
    if predicted_file:
        predicted_df = pd.read_csv(predicted_file, sep='\t')
    
    return human_eei_df, mouse_eei_df, egio_df, predicted_df

def create_orthology_mappings(egio_df, min_identity=0.8):
    """
    Create bidirectional orthology mappings between mouse and human exons
    
    Args:
        egio_df: DataFrame with EGIO orthology mappings
        min_identity: Minimum identity score to consider orthology
        
    Returns:
        Dictionary with mappings in both directions
    """
    # Convert 'Iden' column to numeric
    egio_df['Iden'] = pd.to_numeric(egio_df['Iden'], errors='coerce')
    
    # Filter by minimum identity
    filtered_egio = egio_df[egio_df['Iden'] >= min_identity]
    
    # Create mapping dictionaries
    mouse_to_human = {}
    human_to_mouse = {}
    
    for _, row in filtered_egio.iterrows():
        if pd.notna(row['musPos']) and pd.notna(row['hsaPos']):
            mouse_to_human[row['musPos']] = {
                'hsaPos': row['hsaPos'],
                'identity': row['Iden'],
                'type': row['Type']
            }
            
            human_to_mouse[row['hsaPos']] = {
                'musPos': row['musPos'],
                'identity': row['Iden'],
                'type': row['Type']
            }
    
    return {
        'mouse_to_human': mouse_to_human,
        'human_to_mouse': human_to_mouse
    }

def build_eei_networks(human_eei_df, mouse_eei_df):
    """
    Build network representations of the EEI data
    
    Args:
        human_eei_df: DataFrame with human EEIs
        mouse_eei_df: DataFrame with mouse EEIs
        
    Returns:
        Dictionary with human and mouse networks
    """
    # Create human network
    human_network = nx.Graph()
    for _, row in human_eei_df.iterrows():
        exon1 = row['exon1']
        exon2 = row['exon2']
        
        # Add nodes with attributes
        human_network.add_node(exon1, type='exon')
        human_network.add_node(exon2, type='exon')
        
        # Add edge with attributes
        attrs = {}
        for col in human_eei_df.columns:
            if col not in ['exon1', 'exon2']:
                attrs[col] = row[col]
                
        # Add protein pair information
        protein_pair = get_protein_pair(row)
        if protein_pair:
            attrs['protein_pair'] = protein_pair
                
        human_network.add_edge(exon1, exon2, **attrs)
    
    # Create mouse network
    mouse_network = nx.Graph()
    for _, row in mouse_eei_df.iterrows():
        exon1 = row['exon1']
        exon2 = row['exon2']
        
        # Add nodes with attributes
        mouse_network.add_node(exon1, type='exon')
        mouse_network.add_node(exon2, type='exon')
        
        # Add edge with attributes
        attrs = {}
        for col in mouse_eei_df.columns:
            if col not in ['exon1', 'exon2']:
                attrs[col] = row[col]
        
        # Add protein pair information
        protein_pair = get_protein_pair(row)
        if protein_pair:
            attrs['protein_pair'] = protein_pair
                
        mouse_network.add_edge(exon1, exon2, **attrs)
    
    return {
        'human': human_network,
        'mouse': mouse_network
    }

def identify_conserved_eeis(networks, orthology_mappings):
    """
    Identify conserved EEIs between mouse and human
    
    Args:
        networks: Dictionary with human and mouse networks
        orthology_mappings: Dictionary with orthology mappings
        
    Returns:
        Dictionary with conserved EEI statistics
    """
    mouse_network = networks['mouse']
    human_network = networks['human']
    
    mouse_to_human = orthology_mappings['mouse_to_human']
    
    # Count different types of conservation
    total_mouse_eeis = mouse_network.number_of_edges()
    mappable_eeis = 0
    conserved_eeis = 0
    
    # Store conserved pairs for further analysis
    conserved_pairs = []
    
    # Check each mouse EEI
    for mouse_exon1, mouse_exon2 in mouse_network.edges():
        # Check if both exons have human orthologs
        if mouse_exon1 in mouse_to_human and mouse_exon2 in mouse_to_human:
            mappable_eeis += 1
            
            human_exon1 = mouse_to_human[mouse_exon1]['hsaPos']
            human_exon2 = mouse_to_human[mouse_exon2]['hsaPos']
            
            # Check if the corresponding human EEI exists
            if human_network.has_edge(human_exon1, human_exon2):
                conserved_eeis += 1
                
                # Store information about this conserved pair
                mouse_attrs = mouse_network.edges[mouse_exon1, mouse_exon2]
                human_attrs = human_network.edges[human_exon1, human_exon2]
                
                conserved_pairs.append({
                    'mouse_exon1': mouse_exon1,
                    'mouse_exon2': mouse_exon2,
                    'human_exon1': human_exon1,
                    'human_exon2': human_exon2,
                    'mouse_identity1': mouse_to_human[mouse_exon1]['identity'],
                    'mouse_identity2': mouse_to_human[mouse_exon2]['identity'],
                    'mouse_type1': mouse_to_human[mouse_exon1]['type'],
                    'mouse_type2': mouse_to_human[mouse_exon2]['type'],
                    'mouse_coverage1': mouse_attrs.get('exon1_coverage_percent', None),
                    'mouse_coverage2': mouse_attrs.get('exon2_coverage_percent', None),
                    'human_coverage1': human_attrs.get('exon1_coverage_percent', None),
                    'human_coverage2': human_attrs.get('exon2_coverage_percent', None)
                })
    
    conservation_stats = {
        'total_mouse_eeis': total_mouse_eeis,
        'mappable_eeis': mappable_eeis,
        'conserved_eeis': conserved_eeis,
        'mappable_rate': mappable_eeis / total_mouse_eeis * 100 if total_mouse_eeis > 0 else 0,
        'conservation_rate': conserved_eeis / mappable_eeis * 100 if mappable_eeis > 0 else 0,
        'overall_conservation_rate': conserved_eeis / total_mouse_eeis * 100 if total_mouse_eeis > 0 else 0,
        'conserved_pairs': pd.DataFrame(conserved_pairs) if conserved_pairs else None
    }
    
    return conservation_stats

def analyze_interface_coverage(conservation_stats):
    """
    Analyze coverage statistics for conserved EEI interfaces
    
    Args:
        conservation_stats: Dictionary with conservation statistics
        
    Returns:
        Dictionary with coverage analysis results
    """
    if not conservation_stats or 'conserved_pairs' not in conservation_stats or conservation_stats['conserved_pairs'] is None:
        return None
    
    conserved_df = conservation_stats['conserved_pairs']
    
    # Ensure numeric columns
    for col in ['mouse_coverage1', 'mouse_coverage2', 'human_coverage1', 'human_coverage2']:
        if col in conserved_df.columns:
            conserved_df[col] = pd.to_numeric(conserved_df[col], errors='coerce')
    
    # Calculate coverage differences
    if all(col in conserved_df.columns for col in ['mouse_coverage1', 'human_coverage1']):
        conserved_df['coverage_diff1'] = conserved_df['human_coverage1'] - conserved_df['mouse_coverage1']
    
    if all(col in conserved_df.columns for col in ['mouse_coverage2', 'human_coverage2']):
        conserved_df['coverage_diff2'] = conserved_df['human_coverage2'] - conserved_df['mouse_coverage2']
    
    # Calculate average coverage and differences
    coverage_stats = {}
    
    if 'mouse_coverage1' in conserved_df.columns and 'mouse_coverage2' in conserved_df.columns:
        coverage_stats['avg_mouse_coverage'] = (
            conserved_df['mouse_coverage1'].mean() + conserved_df['mouse_coverage2'].mean()
        ) / 2
    
    if 'human_coverage1' in conserved_df.columns and 'human_coverage2' in conserved_df.columns:
        coverage_stats['avg_human_coverage'] = (
            conserved_df['human_coverage1'].mean() + conserved_df['human_coverage2'].mean()
        ) / 2
    
    if 'coverage_diff1' in conserved_df.columns and 'coverage_diff2' in conserved_df.columns:
        coverage_stats['avg_coverage_diff'] = (
            conserved_df['coverage_diff1'].mean() + conserved_df['coverage_diff2'].mean()
        ) / 2
        
    # Count how many conserved EEIs have similar coverage (within 10%)
    if all(col in conserved_df.columns for col in ['coverage_diff1', 'coverage_diff2']):
        similar_coverage_count = sum(
            (abs(conserved_df['coverage_diff1']) <= 10) & 
            (abs(conserved_df['coverage_diff2']) <= 10)
        )
        
        coverage_stats['similar_coverage_count'] = similar_coverage_count
        coverage_stats['similar_coverage_percent'] = similar_coverage_count / len(conserved_df) * 100
    
    return coverage_stats

def analyze_orthology_types_in_conserved(conservation_stats):
    """
    Analyze the types of orthology relationships in conserved EEIs
    
    Args:
        conservation_stats: Dictionary with conservation statistics
        
    Returns:
        Dictionary with orthology type statistics
    """
    if not conservation_stats or 'conserved_pairs' not in conservation_stats or conservation_stats['conserved_pairs'] is None:
        return None
    
    conserved_df = conservation_stats['conserved_pairs']
    
    # Count orthology types
    type_counts = defaultdict(int)
    type_pair_counts = defaultdict(int)
    
    for _, row in conserved_df.iterrows():
        type1 = row.get('mouse_type1', 'unknown')
        type2 = row.get('mouse_type2', 'unknown')
        
        type_counts[type1] += 1
        type_counts[type2] += 1
        
        type_pair = tuple(sorted([type1, type2]))
        type_pair_counts[type_pair] += 1
    
    # Convert to regular dicts
    type_stats = {
        'exon_types': dict(type_counts),
        'type_pairs': dict(type_pair_counts)
    }
    
    return type_stats

def analyze_prediction_accuracy(conservation_stats, predicted_df):
    """
    Analyze how well the predictions captured conserved EEIs
    
    Args:
        conservation_stats: Dictionary with conservation statistics
        predicted_df: DataFrame with predicted EEIs
        
    Returns:
        Dictionary with prediction accuracy statistics
    """
    if not conservation_stats or 'conserved_pairs' not in conservation_stats or conservation_stats['conserved_pairs'] is None or predicted_df is None:
        return None
    
    conserved_df = conservation_stats['conserved_pairs']
    
    # Create sets of EEI pairs for comparison
    conserved_pairs = set()
    for _, row in conserved_df.iterrows():
        conserved_pairs.add((row['human_exon1'], row['human_exon2']))
        conserved_pairs.add((row['human_exon2'], row['human_exon1']))  # Add both directions
    
    predicted_pairs = set()
    for _, row in predicted_df.iterrows():
        predicted_pairs.add((row['exon1'], row['exon2']))
        predicted_pairs.add((row['exon2'], row['exon1']))  # Add both directions
    
    # Calculate overlap
    correctly_predicted = conserved_pairs.intersection(predicted_pairs)
    
    # Account for double-counting of undirected pairs
    accuracy_stats = {
        'conserved_pair_count': len(conserved_pairs) // 2,
        'predicted_pair_count': len(predicted_pairs) // 2,
        'correctly_predicted_count': len(correctly_predicted) // 2
    }
    
    # Calculate percentages
    if accuracy_stats['conserved_pair_count'] > 0:
        accuracy_stats['recall'] = accuracy_stats['correctly_predicted_count'] / accuracy_stats['conserved_pair_count'] * 100
    else:
        accuracy_stats['recall'] = 0
        
    if accuracy_stats['predicted_pair_count'] > 0:
        accuracy_stats['precision'] = accuracy_stats['correctly_predicted_count'] / accuracy_stats['predicted_pair_count'] * 100
    else:
        accuracy_stats['precision'] = 0
    
    # Calculate F1 score
    if accuracy_stats['precision'] + accuracy_stats['recall'] > 0:
        accuracy_stats['f1_score'] = 2 * (accuracy_stats['precision'] * accuracy_stats['recall']) / (accuracy_stats['precision'] + accuracy_stats['recall'])
    else:
        accuracy_stats['f1_score'] = 0
    
    return accuracy_stats

def create_visualizations(conservation_stats, coverage_stats, type_stats, accuracy_stats, output_dir='conservation_figures'):
    """
    Create visualizations of the conservation analysis results
    
    Args:
        conservation_stats: Dictionary with conservation statistics
        coverage_stats: Dictionary with coverage analysis results
        type_stats: Dictionary with orthology type statistics
        accuracy_stats: Dictionary with prediction accuracy statistics
        output_dir: Directory to save the figures
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Set style
    plt.style.use('seaborn-v0_8-darkgrid')
    sns.set_context("talk")
    
    # 1. Conservation rates
    plt.figure(figsize=(10, 6))
    rate_labels = ['Mappable Rate', 'Conservation Rate', 'Overall Conservation']
    rate_values = [
        conservation_stats['mappable_rate'],
        conservation_stats['conservation_rate'],
        conservation_stats['overall_conservation_rate']
    ]
    
    plt.bar(rate_labels, rate_values, color=['#66B2FF', '#FF9999', '#99CC99'])
    plt.title('EEI Conservation Rates')
    plt.ylabel('Percentage (%)')
    plt.ylim(0, 100)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/conservation_rates.png", dpi=300)
    plt.close()
    
    # 2. EEI counts
    plt.figure(figsize=(10, 6))
    count_labels = ['Total Mouse EEIs', 'Mappable EEIs', 'Conserved EEIs']
    count_values = [
        conservation_stats['total_mouse_eeis'],
        conservation_stats['mappable_eeis'],
        conservation_stats['conserved_eeis']
    ]
    
    plt.bar(count_labels, count_values, color=['#66B2FF', '#FF9999', '#99CC99'])
    plt.title('EEI Conservation Counts')
    plt.ylabel('Count')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/conservation_counts.png", dpi=300)
    plt.close()
    
    # 3. Conservation breakdown pie chart
    plt.figure(figsize=(8, 8))
    pie_labels = ['Non-mappable', 'Mappable but not conserved', 'Conserved']
    pie_values = [
        conservation_stats['total_mouse_eeis'] - conservation_stats['mappable_eeis'],
        conservation_stats['mappable_eeis'] - conservation_stats['conserved_eeis'],
        conservation_stats['conserved_eeis']
    ]
    
    plt.pie(pie_values, labels=pie_labels, autopct='%1.1f%%', startangle=90, 
            colors=['#FF9999', '#FFCC99', '#99CC99'])
    plt.axis('equal')
    plt.title('Mouse EEI Conservation Breakdown')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/conservation_pie.png", dpi=300)
    plt.close()
    
    # 4. Coverage comparison for conserved EEIs
    if coverage_stats and 'conserved_pairs' in conservation_stats and conservation_stats['conserved_pairs'] is not None:
        conserved_df = conservation_stats['conserved_pairs']
        
        if all(col in conserved_df.columns for col in ['mouse_coverage1', 'human_coverage1']):
            plt.figure(figsize=(12, 6))
            
            # Scatterplot of coverage comparison for exon1
            ax1 = plt.subplot(1, 2, 1)
            sns.scatterplot(x='mouse_coverage1', y='human_coverage1', data=conserved_df, ax=ax1)
            ax1.plot([0, 100], [0, 100], 'r--')  # Diagonal line for reference
            ax1.set_title('Exon 1 Coverage Comparison')
            ax1.set_xlabel('Mouse Exon Coverage (%)')
            ax1.set_ylabel('Human Exon Coverage (%)')
            ax1.set_xlim(0, 100)
            ax1.set_ylim(0, 100)
            
            # Scatterplot of coverage comparison for exon2
            ax2 = plt.subplot(1, 2, 2)
            sns.scatterplot(x='mouse_coverage2', y='human_coverage2', data=conserved_df, ax=ax2)
            ax2.plot([0, 100], [0, 100], 'r--')  # Diagonal line for reference
            ax2.set_title('Exon 2 Coverage Comparison')
            ax2.set_xlabel('Mouse Exon Coverage (%)')
            ax2.set_ylabel('Human Exon Coverage (%)')
            ax2.set_xlim(0, 100)
            ax2.set_ylim(0, 100)
            
            plt.tight_layout()
            plt.savefig(f"{output_dir}/coverage_comparison.png", dpi=300)
            plt.close()
        
        # Coverage difference distribution
        if 'coverage_diff1' in conserved_df.columns and 'coverage_diff2' in conserved_df.columns:
            plt.figure(figsize=(12, 6))
            
            # Histogram of coverage differences for exon1
            ax1 = plt.subplot(1, 2, 1)
            sns.histplot(conserved_df['coverage_diff1'], bins=20, kde=True, ax=ax1)
            ax1.axvline(x=0, color='r', linestyle='--')
            ax1.set_title('Exon 1 Coverage Difference')
            ax1.set_xlabel('Human - Mouse Coverage (%)')
            ax1.set_ylabel('Count')
            
            # Histogram of coverage differences for exon2
            ax2 = plt.subplot(1, 2, 2)
            sns.histplot(conserved_df['coverage_diff2'], bins=20, kde=True, ax=ax2)
            ax2.axvline(x=0, color='r', linestyle='--')
            ax2.set_title('Exon 2 Coverage Difference')
            ax2.set_xlabel('Human - Mouse Coverage (%)')
            ax2.set_ylabel('Count')
            
            plt.tight_layout()
            plt.savefig(f"{output_dir}/coverage_difference.png", dpi=300)
            plt.close()
    
    # 5. Orthology types in conserved EEIs
    if type_stats and 'exon_types' in type_stats:
        plt.figure(figsize=(10, 6))
        types = list(type_stats['exon_types'].keys())
        type_counts = list(type_stats['exon_types'].values())
        
        # Sort by count
        sorted_indices = np.argsort(type_counts)[::-1]
        sorted_types = [types[i] for i in sorted_indices]
        sorted_counts = [type_counts[i] for i in sorted_indices]
        
        plt.bar(sorted_types, sorted_counts)
        plt.title('Orthology Types in Conserved EEIs')
        plt.ylabel('Count')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(f"{output_dir}/conserved_orthology_types.png", dpi=300)
        plt.close()
        
        # Orthology type pairs
        if 'type_pairs' in type_stats:
            plt.figure(figsize=(12, 6))
            pair_types = list(type_stats['type_pairs'].keys())
            pair_counts = list(type_stats['type_pairs'].values())
            
            # Create readable labels
            pair_labels = [f"{p[0]}-{p[1]}" for p in pair_types]
            
            # Sort by count
            sorted_indices = np.argsort(pair_counts)[::-1]
            sorted_pairs = [pair_labels[i] for i in sorted_indices]
            sorted_pair_counts = [pair_counts[i] for i in sorted_indices]
            
            # Limit to top 10 pairs
            if len(sorted_pairs) > 10:
                sorted_pairs = sorted_pairs[:10]
                sorted_pair_counts = sorted_pair_counts[:10]
            
            plt.bar(sorted_pairs, sorted_pair_counts)
            plt.title('Top Orthology Type Pairs in Conserved EEIs')
            plt.ylabel('Count')
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.savefig(f"{output_dir}/conserved_type_pairs.png", dpi=300)
            plt.close()
    
    # 6. Prediction accuracy metrics
    if accuracy_stats:
        plt.figure(figsize=(10, 6))
        metrics = ['Recall', 'Precision', 'F1 Score']
        values = [
            accuracy_stats['recall'],
            accuracy_stats['precision'],
            accuracy_stats['f1_score']
        ]
        
        plt.bar(metrics, values, color=['#66B2FF', '#FF9999', '#99CC99'])
        plt.title('Prediction Accuracy Metrics')
        plt.ylabel('Percentage (%)')
        plt.ylim(0, 100)
        plt.tight_layout()
        plt.savefig(f"{output_dir}/prediction_accuracy.png", dpi=300)
        plt.close()
        
        # Prediction counts
        plt.figure(figsize=(10, 6))
        count_labels = ['Conserved EEIs', 'Predicted EEIs', 'Correctly Predicted']
        count_values = [
            accuracy_stats['conserved_pair_count'],
            accuracy_stats['predicted_pair_count'],
            accuracy_stats['correctly_predicted_count']
        ]
        
        plt.bar(count_labels, count_values, color=['#66B2FF', '#FF9999', '#99CC99'])
        plt.title('Prediction Counts')
        plt.ylabel('Count')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/prediction_counts.png", dpi=300)
        plt.close()

def analyze_conservation_by_type(conservation_stats, type_stats):
    """Analyze conservation patterns for different relationship types"""
    if not conservation_stats or 'conserved_pairs' not in conservation_stats:
        return None
        
    conserved_df = conservation_stats['conserved_pairs']
    
    # Group by type pairs
    if 'mouse_type1' in conserved_df.columns and 'mouse_type2' in conserved_df.columns:
        # Create relationship pattern column
        conserved_df['relationship_pattern'] = conserved_df.apply(
            lambda x: f"{x['mouse_type1']}-{x['mouse_type2']}", axis=1
        )
        
        # Calculate conservation stats by pattern
        pattern_stats = {}
        for pattern, group in conserved_df.groupby('relationship_pattern'):
            avg_coverage = None
            
            # Check if coverage columns exist
            if all(col in group.columns for col in ['mouse_coverage1', 'mouse_coverage2', 'human_coverage1', 'human_coverage2']):
                # Convert to numeric if needed
                for col in ['mouse_coverage1', 'mouse_coverage2', 'human_coverage1', 'human_coverage2']:
                    group[col] = pd.to_numeric(group[col], errors='coerce')
                
                # Calculate average coverage
                avg_coverage = (
                    group['mouse_coverage1'].mean() + 
                    group['mouse_coverage2'].mean() +
                    group['human_coverage1'].mean() + 
                    group['human_coverage2'].mean()
                ) / 4
            
            pattern_stats[pattern] = {
                'count': len(group),
                'avg_coverage': avg_coverage
            }
            
        return {'pattern_conservation': pattern_stats}
    
    return None

def run_conservation_analysis(human_eei_file, mouse_eei_file, egio_file, predicted_file=None, output_dir='conservation_results'):
    """
    Run a complete conservation analysis between mouse and human EEIs
    
    Args:
        human_eei_file: Path to human EEI network file
        mouse_eei_file: Path to mouse EEI network file
        egio_file: Path to EGIO orthology mapping file
        predicted_file: Optional path to predicted EEIs file
        output_dir: Directory to save the results
    """
    print(f"Starting conservation analysis between mouse and human EEIs")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Load data
    human_eei_df, mouse_eei_df, egio_df, predicted_df = load_data(
        human_eei_file, mouse_eei_file, egio_file, predicted_file
    )
    
    # Create orthology mappings
    print("Creating orthology mappings...")
    orthology_mappings = create_orthology_mappings(egio_df)
    
    # Build networks
    print("Building EEI networks...")
    networks = build_eei_networks(human_eei_df, mouse_eei_df)
    
    # Identify conserved EEIs
    print("Identifying conserved EEIs...")
    conservation_stats = identify_conserved_eeis(networks, orthology_mappings)
    
    # Analyze interface coverage
    print("Analyzing interface coverage...")
    coverage_stats = analyze_interface_coverage(conservation_stats)
    
    # Analyze orthology types
    print("Analyzing orthology types...")
    type_stats = analyze_orthology_types_in_conserved(conservation_stats)
    
    # Analyze conservation by type
    print("Analyzing conservation by relationship type...")
    type_conservation_stats = analyze_conservation_by_type(conservation_stats, type_stats)
    
    # Analyze prediction accuracy if predictions are provided
    accuracy_stats = None
    if predicted_df is not None:
        print("Analyzing prediction accuracy...")
        accuracy_stats = analyze_prediction_accuracy(conservation_stats, predicted_df)
    
    # Create visualizations
    print("Creating visualizations...")
    create_visualizations(conservation_stats, coverage_stats, type_stats, accuracy_stats, output_dir)
    
    # Create additional visualization for conservation by type
    if type_conservation_stats and 'pattern_conservation' in type_conservation_stats:
        patterns = list(type_conservation_stats['pattern_conservation'].keys())
        counts = [type_conservation_stats['pattern_conservation'][p]['count'] for p in patterns]
        
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
        plt.title('Conserved EEIs by Relationship Pattern')
        plt.xlabel('Relationship Pattern')
        plt.ylabel('Count')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/conservation_by_pattern.png", dpi=300)
        plt.close()
        
        # Create visualization of coverage by pattern if available
        has_coverage = all('avg_coverage' in type_conservation_stats['pattern_conservation'][p] and 
                         type_conservation_stats['pattern_conservation'][p]['avg_coverage'] is not None 
                         for p in top_patterns)
        
        if has_coverage:
            avg_coverages = [type_conservation_stats['pattern_conservation'][p]['avg_coverage'] for p in top_patterns]
            
            plt.figure(figsize=(12, 6))
            bars = plt.bar(top_patterns, avg_coverages)
            plt.title('Average Interface Coverage by Relationship Pattern')
            plt.xlabel('Relationship Pattern')
            plt.ylabel('Average Coverage (%)')
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            plt.savefig(f"{output_dir}/coverage_by_pattern.png", dpi=300)
            plt.close()
    
    # Print summary
    print("\nConservation Analysis Summary:")
    print(f"  Total mouse EEIs: {conservation_stats['total_mouse_eeis']}")
    print(f"  Mappable EEIs: {conservation_stats['mappable_eeis']} ({conservation_stats['mappable_rate']:.2f}%)")
    print(f"  Conserved EEIs: {conservation_stats['conserved_eeis']} ({conservation_stats['conservation_rate']:.2f}% of mappable)")
    print(f"  Overall conservation rate: {conservation_stats['overall_conservation_rate']:.2f}%")
    
    if coverage_stats:
        print("\nCoverage Analysis:")
        if 'avg_mouse_coverage' in coverage_stats:
            print(f"  Average mouse interface coverage: {coverage_stats['avg_mouse_coverage']:.2f}%")
        if 'avg_human_coverage' in coverage_stats:
            print(f"  Average human interface coverage: {coverage_stats['avg_human_coverage']:.2f}%")
        if 'similar_coverage_percent' in coverage_stats:
            print(f"  Interfaces with similar coverage: {coverage_stats['similar_coverage_percent']:.2f}%")
    
    if type_conservation_stats and 'pattern_conservation' in type_conservation_stats:
        print("\nConservation by Relationship Pattern:")
        # Print top 5 patterns by count
        top_5 = sorted(type_conservation_stats['pattern_conservation'].items(), 
                       key=lambda x: x[1]['count'], reverse=True)[:5]
        for pattern, stats in top_5:
            print(f"  {pattern}: {stats['count']} conserved EEIs")
    
    if accuracy_stats:
        print("\nPrediction Accuracy:")
        print(f"  Recall: {accuracy_stats['recall']:.2f}%")
        print(f"  Precision: {accuracy_stats['precision']:.2f}%")
        print(f"  F1 Score: {accuracy_stats['f1_score']:.2f}%")
    
    print("\nFull conservation analysis results saved to:", output_dir)
    
    # Return statistics for further analysis
    return {
        'conservation_stats': conservation_stats,
        'coverage_stats': coverage_stats,
        'type_stats': type_stats,
        'type_conservation_stats': type_conservation_stats,
        'accuracy_stats': accuracy_stats
    }

if __name__ == "__main__":
    # Example usage
    run_conservation_analysis(
        human_eei_file="path/to/human_eei_network.tsv",
        mouse_eei_file="path/to/mouse_eei_network.tsv",
        egio_file="path/to/egio_output.tsv",
        predicted_file="path/to/predicted_eeis.tsv"
    )