import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from collections import defaultdict
import os

def load_prediction_data(predicted_file, human_eei_file, mouse_eei_file, egio_file):
    """
    Load all datasets needed for evaluation
    
    Args:
        predicted_file: Path to the predicted EEIs file
        human_eei_file: Path to the human EEI network file
        mouse_eei_file: Path to the mouse EEI network file
        egio_file: Path to the EGIO orthology mapping file
        
    Returns:
        Tuple of DataFrames (predicted_df, human_eei_df, mouse_eei_df, egio_df)
    """
    predicted_df = pd.read_csv(predicted_file, sep='\t')
    human_eei_df = pd.read_csv(human_eei_file, sep='\t')
    mouse_eei_df = pd.read_csv(mouse_eei_file, sep='\t')
    egio_df = pd.read_csv(egio_file, sep='\t')
    
    print(f"Loaded {len(predicted_df)} predicted EEIs")
    print(f"Loaded {len(human_eei_df)} human EEIs")
    print(f"Loaded {len(mouse_eei_df)} mouse EEIs")
    print(f"Loaded {len(egio_df)} orthology mappings")
    
    return predicted_df, human_eei_df, mouse_eei_df, egio_df

def evaluate_coverage(predicted_df, human_eei_df, mouse_eei_df, egio_df):
    """
    Evaluate how well the predictions cover the potential exon-exon interactions
    
    Args:
        predicted_df: DataFrame with predicted EEIs
        human_eei_df: DataFrame with human EEIs
        mouse_eei_df: DataFrame with mouse EEIs
        egio_df: DataFrame with orthology mappings
        
    Returns:
        Dictionary with coverage statistics
    """
    # Create mappings
    mouse_to_human = {}
    human_to_mouse = {}
    
    # Only consider high-quality orthology mappings (identity >= 0.8)
    egio_df['Iden'] = pd.to_numeric(egio_df['Iden'], errors='coerce')
    high_quality = egio_df[egio_df['Iden'] >= 0.8]
    
    for _, row in high_quality.iterrows():
        if pd.notna(row['musPos']) and pd.notna(row['hsaPos']):
            mouse_to_human[row['musPos']] = row['hsaPos']
            human_to_mouse[row['hsaPos']] = row['musPos']
    
    # Count potential EEIs (mouse EEIs that could be transferred to human)
    potential_eeis = set()
    for _, row in mouse_eei_df.iterrows():
        if row['exon1'] in mouse_to_human and row['exon2'] in mouse_to_human:
            human_exon1 = mouse_to_human[row['exon1']]
            human_exon2 = mouse_to_human[row['exon2']]
            potential_eeis.add((human_exon1, human_exon2))
            potential_eeis.add((human_exon2, human_exon1))  # Add both directions
    
    # Convert to undirected pairs (count each pair only once)
    potential_eeis_undirected = set()
    for e1, e2 in potential_eeis:
        if (e1, e2) not in potential_eeis_undirected and (e2, e1) not in potential_eeis_undirected:
            potential_eeis_undirected.add((e1, e2))
    
    # Extract predicted EEI pairs
    predicted_pairs = set()
    for _, row in predicted_df.iterrows():
        predicted_pairs.add((row['exon1'], row['exon2']))
        predicted_pairs.add((row['exon2'], row['exon1']))  # Add both directions
    
    # Convert to undirected pairs
    predicted_pairs_undirected = set()
    for e1, e2 in predicted_pairs:
        if (e1, e2) not in predicted_pairs_undirected and (e2, e1) not in predicted_pairs_undirected:
            predicted_pairs_undirected.add((e1, e2))
    
    # Extract existing human EEI pairs
    human_pairs = set()
    for _, row in human_eei_df.iterrows():
        human_pairs.add((row['exon1'], row['exon2']))
        human_pairs.add((row['exon2'], row['exon1']))  # Add both directions
    
    # Convert to undirected pairs
    human_pairs_undirected = set()
    for e1, e2 in human_pairs:
        if (e1, e2) not in human_pairs_undirected and (e2, e1) not in human_pairs_undirected:
            human_pairs_undirected.add((e1, e2))
    
    # Calculate statistics
    coverage_stats = {
        'potential_eeis': len(potential_eeis_undirected),
        'predicted_eeis': len(predicted_pairs_undirected),
        'existing_human_eeis': len(human_pairs_undirected),
        'coverage_rate': len(predicted_pairs_undirected) / len(potential_eeis_undirected) * 100 if potential_eeis_undirected else 0
    }
    
    # Calculate overlap
    predicted_in_existing = predicted_pairs_undirected.intersection(human_pairs_undirected)
    potential_in_existing = potential_eeis_undirected.intersection(human_pairs_undirected)
    
    coverage_stats['predicted_in_existing'] = len(predicted_in_existing)
    coverage_stats['potential_in_existing'] = len(potential_in_existing)
    
    if human_pairs_undirected:
        coverage_stats['predicted_overlap_rate'] = len(predicted_in_existing) / len(predicted_pairs_undirected) * 100 if predicted_pairs_undirected else 0
        coverage_stats['potential_overlap_rate'] = len(potential_in_existing) / len(potential_eeis_undirected) * 100 if potential_eeis_undirected else 0
    
    return coverage_stats

def analyze_protein_interactions(predicted_df, human_eei_df, mouse_eei_df):
    """
    Analyze the protein interactions involved in the predicted EEIs
    
    Args:
        predicted_df: DataFrame with predicted EEIs
        human_eei_df: DataFrame with human EEIs
        mouse_eei_df: DataFrame with mouse EEIs
        
    Returns:
        Dictionary with protein-level statistics
    """
    # Extract protein pairs from mouse EEIs
    mouse_protein_pairs = set()
    for _, row in mouse_eei_df.iterrows():
        # Check if the required protein columns exist
        protein1 = None
        protein2 = None
        
        if 'protein1' in row:
            protein1 = row['protein1']
        
        # Try 'protein1.1' (CONTACT format) and 'protein2' (PISA format)
        if 'protein1.1' in row:
            protein2 = row['protein1.1']
        elif 'protein2' in row:
            protein2 = row['protein2']
            
        if protein1 and protein2:
            protein_pair = tuple(sorted([protein1, protein2]))
            mouse_protein_pairs.add(protein_pair)
    
    # Extract protein pairs from human EEIs
    human_protein_pairs = set()
    for _, row in human_eei_df.iterrows():
        # Check if the required protein columns exist
        protein1 = None
        protein2 = None
        
        if 'protein1' in row:
            protein1 = row['protein1']
        
        # Try 'protein1.1' (CONTACT format) and 'protein2' (PISA format)
        if 'protein1.1' in row:
            protein2 = row['protein1.1']
        elif 'protein2' in row:
            protein2 = row['protein2']
            
        if protein1 and protein2:
            protein_pair = tuple(sorted([protein1, protein2]))
            human_protein_pairs.add(protein_pair)
    
    # Extract protein pairs from predicted EEIs
    predicted_protein_pairs = set()
    for _, row in predicted_df.iterrows():
        # Check if the required protein columns exist
        protein1 = None
        protein2 = None
        
        if 'protein1' in row:
            protein1 = row['protein1']
        
        # Try 'protein1.1' (CONTACT format) and 'protein2' (PISA format)
        if 'protein1.1' in row:
            protein2 = row['protein1.1']
        elif 'protein2' in row:
            protein2 = row['protein2']
            
        if protein1 and protein2:
            protein_pair = tuple(sorted([protein1, protein2]))
            predicted_protein_pairs.add(protein_pair)
    
    # Calculate statistics
    protein_stats = {
        'mouse_protein_pairs': len(mouse_protein_pairs),
        'human_protein_pairs': len(human_protein_pairs),
        'predicted_protein_pairs': len(predicted_protein_pairs) if predicted_protein_pairs else 0
    }
    
    # Calculate overlaps
    if predicted_protein_pairs:
        predicted_in_human = predicted_protein_pairs.intersection(human_protein_pairs)
        protein_stats['predicted_protein_overlap'] = len(predicted_in_human)
        protein_stats['predicted_protein_overlap_rate'] = len(predicted_in_human) / len(predicted_protein_pairs) * 100
    
    return protein_stats

def analyze_confidence_vs_accuracy(predicted_df, human_eei_df):
    """
    Analyze how confidence scores relate to prediction accuracy
    
    Args:
        predicted_df: DataFrame with predicted EEIs
        human_eei_df: DataFrame with human EEIs
        
    Returns:
        DataFrame with confidence bins and accuracy metrics
    """
    # Extract human EEI pairs
    human_pairs = set()
    for _, row in human_eei_df.iterrows():
        human_pairs.add((row['exon1'], row['exon2']))
        human_pairs.add((row['exon2'], row['exon1']))  # Add both directions
    
    # Group predictions by confidence score
    if 'confidence' not in predicted_df.columns:
        print("Confidence scores not found in predictions")
        return None
    
    # Create bins for confidence scores
    bins = [0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0]
    bin_labels = [f"{bins[i]:.2f}-{bins[i+1]:.2f}" for i in range(len(bins)-1)]
    bin_labels.insert(0, f"<{bins[0]:.2f}")
    bin_labels.append(f"≥{bins[-1]:.2f}")
    
    # Initialize results
    results = {
        'bin': bin_labels,
        'count': [0] * len(bin_labels),
        'validated': [0] * len(bin_labels),
        'accuracy': [0.0] * len(bin_labels)
    }
    
    # Count predictions and validated predictions in each bin
    for _, row in predicted_df.iterrows():
        pair = (row['exon1'], row['exon2'])
        confidence = row['confidence']
        
        # Determine bin index
        bin_index = 0
        if confidence >= bins[-1]:
            bin_index = len(bins)
        else:
            for i, threshold in enumerate(bins):
                if confidence < threshold:
                    bin_index = i
                    break
        
        # Update counts
        results['count'][bin_index] += 1
        
        # Check if prediction is validated (exists in human EEIs)
        if pair in human_pairs or (row['exon2'], row['exon1']) in human_pairs:
            results['validated'][bin_index] += 1
    
    # Calculate accuracy for each bin
    for i in range(len(bin_labels)):
        if results['count'][i] > 0:
            results['accuracy'][i] = results['validated'][i] / results['count'][i] * 100
    
    return pd.DataFrame(results)

def analyze_orthology_types(predicted_df, egio_df):
    """
    Analyze the types of orthology relationships in the predictions
    
    Args:
        predicted_df: DataFrame with predicted EEIs
        egio_df: DataFrame with orthology mappings
        
    Returns:
        Dictionary with orthology type statistics
    """
    # Create mapping from exon to orthology type
    exon_to_type = {}
    for _, row in egio_df.iterrows():
        if pd.notna(row['hsaPos']) and pd.notna(row['Type']):
            exon_to_type[row['hsaPos']] = row['Type']
    
    # Count orthology types in predictions
    type_counts = defaultdict(int)
    type_pair_counts = defaultdict(int)
    
    for _, row in predicted_df.iterrows():
        exon1_type = exon_to_type.get(row['exon1'], 'unknown')
        exon2_type = exon_to_type.get(row['exon2'], 'unknown')
        
        type_counts[exon1_type] += 1
        type_counts[exon2_type] += 1
        
        type_pair = tuple(sorted([exon1_type, exon2_type]))
        type_pair_counts[type_pair] += 1
    
    # Convert to regular dicts
    type_stats = {
        'exon_types': dict(type_counts),
        'type_pairs': dict(type_pair_counts)
    }
    
    return type_stats

def create_visualizations(predicted_df, coverage_stats, protein_stats, confidence_df, type_stats, output_dir='eval_figures'):
    """
    Create visualizations of the evaluation results
    
    Args:
        predicted_df: DataFrame with predicted EEIs
        coverage_stats: Dictionary with coverage statistics
        protein_stats: Dictionary with protein-level statistics
        confidence_df: DataFrame with confidence vs accuracy results
        type_stats: Dictionary with orthology type statistics
        output_dir: Directory to save the figures
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Set style
    plt.style.use('seaborn-v0_8-darkgrid')
    sns.set_context("talk")
    
    # 1. Coverage analysis
    plt.figure(figsize=(10, 6))
    coverage_labels = ['Potential EEIs', 'Predicted EEIs', 'Existing Human EEIs']
    coverage_values = [
        coverage_stats['potential_eeis'],
        coverage_stats['predicted_eeis'],
        coverage_stats['existing_human_eeis']
    ]
    
    plt.bar(coverage_labels, coverage_values, color=['#66B2FF', '#FF9999', '#99CC99'])
    plt.title('EEI Network Coverage')
    plt.ylabel('Count')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/eei_coverage.png", dpi=300)
    plt.close()
    
    # 2. Coverage percentages
    plt.figure(figsize=(8, 6))
    percent_labels = ['Coverage Rate', 'Predicted in Existing', 'Potential in Existing']
    percent_values = [
        coverage_stats['coverage_rate'],
        coverage_stats.get('predicted_overlap_rate', 0),
        coverage_stats.get('potential_overlap_rate', 0)
    ]
    
    plt.bar(percent_labels, percent_values, color=['#66B2FF', '#FF9999', '#99CC99'])
    plt.title('EEI Coverage Percentages')
    plt.ylabel('Percentage')
    plt.ylim(0, 100)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/coverage_percentages.png", dpi=300)
    plt.close()
    
    # 3. Protein interaction coverage
    if protein_stats:
        plt.figure(figsize=(10, 6))
        protein_labels = ['Mouse Protein Pairs', 'Human Protein Pairs', 'Predicted Protein Pairs']
        protein_values = [
            protein_stats['mouse_protein_pairs'],
            protein_stats['human_protein_pairs'],
            protein_stats['predicted_protein_pairs']
        ]
        
        plt.bar(protein_labels, protein_values, color=['#66B2FF', '#FF9999', '#99CC99'])
        plt.title('Protein Interaction Coverage')
        plt.ylabel('Count')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/protein_coverage.png", dpi=300)
        plt.close()
        
        # Protein overlap percentage
        if 'predicted_protein_overlap_rate' in protein_stats:
            plt.figure(figsize=(6, 6))
            labels = ['Novel', 'Known']
            sizes = [
                100 - protein_stats['predicted_protein_overlap_rate'],
                protein_stats['predicted_protein_overlap_rate']
            ]
            
            plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90, colors=['#66B2FF', '#FF9999'])
            plt.axis('equal')
            plt.title('Predicted Protein Interactions Overlap')
            plt.tight_layout()
            plt.savefig(f"{output_dir}/protein_overlap_pie.png", dpi=300)
            plt.close()
    
    # 4. Confidence vs. Accuracy
    if confidence_df is not None:
        plt.figure(figsize=(12, 6))
        
        # Bar chart for counts
        ax1 = plt.subplot(1, 2, 1)
        sns.barplot(x='bin', y='count', data=confidence_df, ax=ax1)
        ax1.set_title('Predictions by Confidence Range')
        ax1.set_xlabel('Confidence Range')
        ax1.set_ylabel('Number of Predictions')
        ax1.tick_params(axis='x', rotation=45)
        
        # Line chart for accuracy
        ax2 = plt.subplot(1, 2, 2)
        sns.lineplot(x='bin', y='accuracy', data=confidence_df, marker='o', ax=ax2)
        ax2.set_title('Accuracy by Confidence Range')
        ax2.set_xlabel('Confidence Range')
        ax2.set_ylabel('Accuracy (%)')
        ax2.set_ylim(0, 100)
        ax2.tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        plt.savefig(f"{output_dir}/confidence_vs_accuracy.png", dpi=300)
        plt.close()
    
    # 5. Orthology types
    if type_stats and 'exon_types' in type_stats:
        plt.figure(figsize=(10, 6))
        types = list(type_stats['exon_types'].keys())
        type_counts = list(type_stats['exon_types'].values())
        
        # Sort by count
        sorted_indices = np.argsort(type_counts)[::-1]
        sorted_types = [types[i] for i in sorted_indices]
        sorted_counts = [type_counts[i] for i in sorted_indices]
        
        plt.bar(sorted_types, sorted_counts)
        plt.title('Orthology Types in Predictions')
        plt.ylabel('Count')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(f"{output_dir}/orthology_types.png", dpi=300)
        plt.close()
        
        # Type pairs chart
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
            plt.title('Top Orthology Type Pairs in Predictions')
            plt.ylabel('Count')
            plt.xticks(rotation=45)
            plt.tight_layout()
            plt.savefig(f"{output_dir}/orthology_type_pairs.png", dpi=300)
            plt.close()
    
    # 6. Confidence score distribution
    if 'confidence' in predicted_df.columns:
        plt.figure(figsize=(10, 6))
        sns.histplot(predicted_df['confidence'], bins=20, kde=True)
        plt.title('Prediction Confidence Distribution')
        plt.xlabel('Confidence Score')
        plt.ylabel('Count')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/confidence_distribution.png", dpi=300)
        plt.close()
        
        # 7. Identity score distributions
        if 'identity1' in predicted_df.columns and 'identity2' in predicted_df.columns:
            plt.figure(figsize=(12, 6))
            
            # Identity1 distribution
            ax1 = plt.subplot(1, 2, 1)
            sns.histplot(predicted_df['identity1'], bins=20, kde=True, ax=ax1)
            ax1.set_title('Exon 1 Identity Score Distribution')
            ax1.set_xlabel('Identity Score')
            ax1.set_ylabel('Count')
            
            # Identity2 distribution
            ax2 = plt.subplot(1, 2, 2)
            sns.histplot(predicted_df['identity2'], bins=20, kde=True, ax=ax2)
            ax2.set_title('Exon 2 Identity Score Distribution')
            ax2.set_xlabel('Identity Score')
            ax2.set_ylabel('Count')
            
            plt.tight_layout()
            plt.savefig(f"{output_dir}/identity_distributions.png", dpi=300)
            plt.close()

def analyze_orthology_types(predicted_df, egio_df):
    """Improved function to better analyze relationship types"""
    # Create mapping from exon to orthology type
    exon_to_type = {}
    for _, row in egio_df.iterrows():
        if pd.notna(row['hsaPos']) and pd.notna(row['Type']):
            exon_to_type[row['hsaPos']] = row['Type']
    
    # Check if type columns already exist
    if 'type1' in predicted_df.columns and 'type2' in predicted_df.columns:
        # Use existing types
        relationship_patterns = []
        for _, row in predicted_df.iterrows():
            type1 = row['type1']
            type2 = row['type2']
            relationship_patterns.append(f"{type1}-{type2}")
    else:
        # Add types from EGIO mapping
        relationship_patterns = []
        for _, row in predicted_df.iterrows():
            type1 = exon_to_type.get(row['exon1'], 'unknown')
            type2 = exon_to_type.get(row['exon2'], 'unknown')
            relationship_patterns.append(f"{type1}-{type2}")
    
    predicted_df_with_patterns = predicted_df.copy()
    predicted_df_with_patterns['relationship_pattern'] = relationship_patterns
    
    # Count by pattern
    pattern_counts = predicted_df_with_patterns['relationship_pattern'].value_counts().to_dict()
    
    # Enhanced analysis for different patterns
    pattern_confidence = {}
    if 'confidence' in predicted_df_with_patterns.columns:
        for pattern in pattern_counts.keys():
            subset = predicted_df_with_patterns[predicted_df_with_patterns['relationship_pattern'] == pattern]
            pattern_confidence[pattern] = {
                'mean': subset['confidence'].mean(),
                'median': subset['confidence'].median(),
                'std': subset['confidence'].std(),
                'count': len(subset)
            }
            
    return {
        'pattern_counts': pattern_counts,
        'pattern_confidence': pattern_confidence,
        'predicted_df_with_patterns': predicted_df_with_patterns
    }

def compare_prediction_sets(predicted_df1, predicted_df2, labels=["With Threshold", "Without Threshold"]):
    """Compare two sets of predicted EEIs"""
    print(f"Comparing prediction sets: {labels[0]} vs {labels[1]}")
    
    # Create sets for easy comparison
    set1 = set()
    for _, row in predicted_df1.iterrows():
        pair = (row['exon1'], row['exon2'])
        set1.add(pair)
        set1.add((row['exon2'], row['exon1']))  # Add both directions for undirected network
    
    set2 = set()
    for _, row in predicted_df2.iterrows():
        pair = (row['exon1'], row['exon2'])
        set2.add(pair)
        set2.add((row['exon2'], row['exon1']))  # Add both directions
    
    # Find common and unique predictions
    common = set1.intersection(set2)
    unique_to_set1 = set1.difference(set2)
    unique_to_set2 = set2.difference(set1)
    
    # Divide by 2 to account for both directions
    results = {
        f"{labels[0]}_total": len(set1) // 2,
        f"{labels[1]}_total": len(set2) // 2,
        "common": len(common) // 2,
        f"unique_to_{labels[0]}": len(unique_to_set1) // 2,
        f"unique_to_{labels[1]}": len(unique_to_set2) // 2
    }
    
    # Add percentages
    results[f"{labels[0]}_overlap_percent"] = (len(common) / len(set1) * 100) if set1 else 0
    results[f"{labels[1]}_overlap_percent"] = (len(common) / len(set2) * 100) if set2 else 0
    
    return results

def analyze_clusters_in_predictions(predicted_df, human_eei_df, min_cluster_size=3):
    """Find clusters of predicted EEIs that might represent conserved protein complexes"""
    import networkx as nx
    from networkx.algorithms import community
    
    print(f"Performing cluster analysis (min cluster size: {min_cluster_size})...")
    
    # Create network from predictions
    G = nx.Graph()
    for _, row in predicted_df.iterrows():
        G.add_edge(row['exon1'], row['exon2'])
    
    # Find communities/clusters
    print("Detecting communities in the prediction network...")
    communities = list(community.greedy_modularity_communities(G))
    print(f"Found {len(communities)} communities/clusters")
    
    # Filter by minimum size
    large_communities = [c for c in communities if len(c) >= min_cluster_size]
    print(f"Found {len(large_communities)} clusters with size >= {min_cluster_size}")
    
    # Analyze top clusters
    cluster_stats = []
    for i, cluster in enumerate(large_communities[:10]):  # Analyze top 10 large clusters
        cluster_size = len(cluster)
        
        # Get confidence scores for edges in this cluster
        cluster_eeis = []
        avg_confidence = 0
        confidence_count = 0
        
        for _, row in predicted_df.iterrows():
            if row['exon1'] in cluster and row['exon2'] in cluster:
                cluster_eeis.append((row['exon1'], row['exon2']))
                if 'confidence' in row and pd.notna(row['confidence']):
                    avg_confidence += row['confidence']
                    confidence_count += 1
        
        if confidence_count > 0:
            avg_confidence /= confidence_count
            
        # Check how many are already known in human EEIs
        known_count = 0
        for e1, e2 in cluster_eeis:
            for _, row in human_eei_df.iterrows():
                if (row['exon1'] == e1 and row['exon2'] == e2) or (row['exon1'] == e2 and row['exon2'] == e1):
                    known_count += 1
                    break
        
        cluster_stats.append({
            'cluster_id': i+1,
            'size': cluster_size,
            'eei_count': len(cluster_eeis),
            'avg_confidence': avg_confidence if confidence_count > 0 else "N/A",
            'known_count': known_count,
            'known_ratio': known_count / len(cluster_eeis) if cluster_eeis else 0
        })
    
    return {
        'total_clusters': len(communities),
        'large_clusters': len(large_communities),
        'cluster_stats': cluster_stats
    }

def visualize_clusters(cluster_results, output_dir):
    """Create visualizations for cluster analysis results"""
    if not cluster_results or 'cluster_stats' not in cluster_results:
        return
        
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Set style
    plt.style.use('seaborn-v0_8-darkgrid')
    sns.set_context("talk")
    
    # Extract data
    cluster_ids = [str(c['cluster_id']) for c in cluster_results['cluster_stats']]
    cluster_sizes = [c['size'] for c in cluster_results['cluster_stats']]
    eei_counts = [c['eei_count'] for c in cluster_results['cluster_stats']]
    known_ratios = [c['known_ratio'] * 100 for c in cluster_results['cluster_stats']]
    
    # Create figure with multiple subplots
    plt.figure(figsize=(15, 10))
    
    # Cluster size
    ax1 = plt.subplot(2, 1, 1)
    bars = ax1.bar(cluster_ids, cluster_sizes)
    ax1.set_title('Cluster Sizes')
    ax1.set_xlabel('Cluster ID')
    ax1.set_ylabel('Number of Exons')
    
    # Known ratio
    ax2 = plt.subplot(2, 1, 2)
    bars = ax2.bar(cluster_ids, known_ratios)
    ax2.set_title('Known EEIs in Clusters (%)')
    ax2.set_xlabel('Cluster ID')
    ax2.set_ylabel('Known Ratio (%)')
    ax2.set_ylim(0, 100)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/cluster_analysis.png", dpi=300)
    plt.close()
    
    # Create a summary table visualization
    plt.figure(figsize=(12, 6))
    ax = plt.subplot(1, 1, 1)
    ax.axis('tight')
    ax.axis('off')
    
    table_data = []
    for c in cluster_results['cluster_stats']:
        if isinstance(c['avg_confidence'], float):
            avg_conf = f"{c['avg_confidence']:.3f}"
        else:
            avg_conf = c['avg_confidence']
            
        table_data.append([
            c['cluster_id'],
            c['size'],
            c['eei_count'],
            avg_conf,
            f"{c['known_count']}/{c['eei_count']}",
            f"{c['known_ratio']*100:.1f}%"
        ])
    
    table = ax.table(
        cellText=table_data,
        colLabels=['Cluster ID', 'Exon Count', 'EEI Count', 'Avg Confidence', 'Known EEIs', 'Known Ratio'],
        loc='center',
        cellLoc='center'
    )
    
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.5)
    
    plt.title('Cluster Analysis Summary')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/cluster_summary_table.png", dpi=300)
    plt.close()

def run_evaluation(predicted_file, human_eei_file, mouse_eei_file, egio_file, 
                  predicted_with_threshold_file=None, analyze_clusters=False, 
                  min_cluster_size=3, output_dir='eval_results'):
    """
    Run a complete evaluation of the EEI predictions
    
    Args:
        predicted_file: Path to the predicted EEIs file
        human_eei_file: Path to the human EEI network file
        mouse_eei_file: Path to the mouse EEI network file
        egio_file: Path to the EGIO orthology mapping file
        predicted_with_threshold_file: Optional path to predicted EEIs with threshold
        analyze_clusters: Whether to perform cluster analysis
        min_cluster_size: Minimum size for clusters to analyze
        output_dir: Directory to save the results
    """
    print(f"Starting evaluation of predictions in {predicted_file}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Load data
    predicted_df, human_eei_df, mouse_eei_df, egio_df = load_prediction_data(
        predicted_file, human_eei_file, mouse_eei_file, egio_file
    )
    
    # Load comparison dataset if provided
    predicted_with_threshold_df = None
    if predicted_with_threshold_file:
        print(f"Loading comparison predictions from {predicted_with_threshold_file}")
        predicted_with_threshold_df = pd.read_csv(predicted_with_threshold_file, sep='\t')
        print(f"Loaded {len(predicted_with_threshold_df)} predictions with threshold")
    
    # Evaluate coverage
    print("Evaluating coverage...")
    coverage_stats = evaluate_coverage(predicted_df, human_eei_df, mouse_eei_df, egio_df)
    
    # Analyze protein interactions
    print("Analyzing protein interactions...")
    protein_stats = analyze_protein_interactions(predicted_df, human_eei_df, mouse_eei_df)
    
    # Analyze confidence vs accuracy
    print("Analyzing confidence vs accuracy...")
    confidence_df = analyze_confidence_vs_accuracy(predicted_df, human_eei_df)
    
    # Analyze orthology types
    print("Analyzing orthology types...")
    type_stats = analyze_orthology_types(predicted_df, egio_df)
    
    # Compare prediction sets if provided
    comparison_stats = None
    if predicted_with_threshold_df is not None:
        print("Comparing prediction sets...")
        comparison_stats = compare_prediction_sets(
            predicted_with_threshold_df, predicted_df, 
            ["With Threshold", "Without Threshold"]
        )
    
    # Analyze clusters if requested
    cluster_stats = None
    if analyze_clusters:
        print("Analyzing clusters in predictions...")
        cluster_stats = analyze_clusters_in_predictions(
            predicted_df, human_eei_df, min_cluster_size
        )
        visualize_clusters(cluster_stats, output_dir)
    
    # Create visualizations
    print("Creating visualizations...")
    create_visualizations(
        predicted_df, coverage_stats, protein_stats, confidence_df, type_stats, output_dir
    )
    
    # Create additional visualizations for orthology types
    if type_stats and 'predicted_df_with_patterns' in type_stats:
        # Create visualization of pattern distributions
        patterns = list(type_stats['pattern_counts'].keys())
        counts = list(type_stats['pattern_counts'].values())
        
        # Sort by count
        sorted_indices = np.argsort(counts)[::-1]
        sorted_patterns = [patterns[i] for i in sorted_indices]
        sorted_counts = [counts[i] for i in sorted_indices]
        
        # Top 15 patterns
        if len(sorted_patterns) > 15:
            top_patterns = sorted_patterns[:15]
            top_counts = sorted_counts[:15]
        else:
            top_patterns = sorted_patterns
            top_counts = sorted_counts
        
        plt.figure(figsize=(14, 8))
        bars = plt.bar(top_patterns, top_counts)
        plt.title('Top Relationship Patterns in Predictions')
        plt.xlabel('Relationship Pattern')
        plt.ylabel('Count')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/relationship_patterns.png", dpi=300)
        plt.close()
    
    # Print summary
    print("\nEvaluation Summary:")
    print(f"  Total predictions: {len(predicted_df)}")
    print(f"  Coverage rate: {coverage_stats['coverage_rate']:.2f}%")
    if 'predicted_overlap_rate' in coverage_stats:
        print(f"  Overlap with existing human EEIs: {coverage_stats['predicted_overlap_rate']:.2f}%")
    if 'confidence' in predicted_df.columns:
        print(f"  Average confidence score: {predicted_df['confidence'].mean():.4f}")
    
    if comparison_stats:
        print("\nPrediction Set Comparison:")
        print(f"  With threshold: {comparison_stats['With Threshold_total']} predictions")
        print(f"  Without threshold: {comparison_stats['Without Threshold_total']} predictions")
        print(f"  Common predictions: {comparison_stats['common']}")
        print(f"  Unique to with threshold: {comparison_stats['unique_to_With Threshold']}")
        print(f"  Unique to without threshold: {comparison_stats['unique_to_Without Threshold']}")
    
    if cluster_stats:
        print("\nCluster Analysis:")
        print(f"  Total clusters found: {cluster_stats['total_clusters']}")
        print(f"  Clusters with size >= {min_cluster_size}: {cluster_stats['large_clusters']}")
        print("  Top clusters:")
        for c in cluster_stats['cluster_stats']:
            print(f"    Cluster {c['cluster_id']}: {c['size']} exons, {c['eei_count']} EEIs, {c['known_ratio']*100:.1f}% known")
    
    print("\nFull evaluation results saved to:", output_dir)
    
    # Save statistics to file
    stats_file = os.path.join(output_dir, 'evaluation_stats.txt')
    with open(stats_file, 'w') as f:
        f.write("=== EEI Prediction Evaluation Statistics ===\n\n")
        
        f.write("Coverage Statistics:\n")
        for key, value in coverage_stats.items():
            f.write(f"  {key}: {value}\n")
        f.write("\n")
        
        f.write("Protein Interaction Statistics:\n")
        for key, value in protein_stats.items():
            f.write(f"  {key}: {value}\n")
        f.write("\n")
        
        if type_stats:
            f.write("Orthology Type Statistics:\n")
            f.write("  Pattern Counts:\n")
            for pattern, count in sorted(type_stats['pattern_counts'].items(), key=lambda x: x[1], reverse=True):
                f.write(f"    {pattern}: {count}\n")
            f.write("\n")
        
        if comparison_stats:
            f.write("Prediction Set Comparison:\n")
            for key, value in comparison_stats.items():
                f.write(f"  {key}: {value}\n")
            f.write("\n")
        
        if cluster_stats:
            f.write("Cluster Analysis:\n")
            f.write(f"  Total clusters: {cluster_stats['total_clusters']}\n")
            f.write(f"  Large clusters (>= {min_cluster_size}): {cluster_stats['large_clusters']}\n")
            f.write("  Top clusters:\n")
            for c in cluster_stats['cluster_stats']:
                f.write(f"    Cluster {c['cluster_id']}: {c['size']} exons, {c['eei_count']} EEIs, {c['known_ratio']*100:.1f}% known\n")
    
    # Return statistics for further analysis
    return {
        'coverage_stats': coverage_stats,
        'protein_stats': protein_stats,
        'type_stats': type_stats,
        'comparison_stats': comparison_stats,
        'cluster_stats': cluster_stats
    }


if __name__ == "__main__":
    # Example usage
    run_evaluation(
        predicted_file="path/to/predicted_eeis.tsv",
        human_eei_file="path/to/human_eei_network.tsv",
        mouse_eei_file="path/to/mouse_eei_network.tsv",
        egio_file="path/to/egio_output.tsv"
    )