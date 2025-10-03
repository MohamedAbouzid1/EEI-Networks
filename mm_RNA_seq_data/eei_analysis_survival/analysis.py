import pandas as pd
import numpy as np
from scipy import stats
from find_expression import find_expression_with_mapping

def create_response_based_survival_proxy(metadata_df, expression_df, eei_df, coord_to_exon):
    """
    Create a survival-like analysis using treatment response patterns from the paper.
    
    Based on the paper:
    - E0771 shows robust regression (good prognosis)
    - 4T1 shows growth stasis (poor prognosis)
    - Early IFN response predicts later regression
    """
    
    # 1. Define response groups based on model and treatment
    def assign_response_group(row):
        """Assign response groups based on paper findings"""
        if row['model'] == 'E0771':
            if row['treatment'] == 'Vehicle':
                return 'poor_control'
            elif 'Cyclophosphamide' in row['treatment']:
                # E0771 + CPA = regression (good outcome)
                return 'good_responder'
        elif row['model'] == '4T1':
            if row['treatment'] == 'Vehicle':
                return 'poor_control'
            elif 'Cyclophosphamide' in row['treatment']:
                # 4T1 + CPA = stasis (moderate outcome)
                return 'moderate_responder'
        return 'unknown'
    
    metadata_df['response_group'] = metadata_df.apply(assign_response_group, axis=1)
    
    # 2. Create pseudo-survival scores based on expected outcomes
    survival_scores = {
        'good_responder': 1.0,      # E0771 + CPA (regression)
        'moderate_responder': 0.5,   # 4T1 + CPA (stasis)
        'poor_control': 0.0          # Vehicle controls
    }
    
    metadata_df['survival_proxy'] = metadata_df['response_group'].map(survival_scores)
    
    return metadata_df

def analyze_eei_ifn_correlation(eei_df, expression_df, metadata_df, coord_to_exon):
    """
    Analyze EEI correlation with IFN response genes from the paper.
    The paper shows early IFN response predicts treatment efficacy.
    """
    
    # Key ISGs from the paper (Figure 6B)
    key_isgs = ['Mx1', 'Cxcl11', 'Rsad2', 'Ifit1', 'Oasl1', 'Cxcl10']
    
    # Additional immune markers from paper
    immune_markers = {
        'cd8_t_cells': ['Cd8a', 'Cd8b1'],
        'nk_cells': ['Nkg46', 'Klrb1c'],
        'ifn_signaling': ['Ifng', 'Ifna', 'Ifnb1'],
        'checkpoint': ['Pdcd1', 'Cd274']  # PD-1, PD-L1
    }
    
    results = []
    
    for idx, eei in eei_df.iterrows():
        mouse_exon1 = eei['mouse_exon1']
        mouse_exon2 = eei['mouse_exon2']
        
        # Get expression for EEI exons
        exon1_expr = find_expression_with_mapping(mouse_exon1, expression_df, coord_to_exon)
        exon2_expr = find_expression_with_mapping(mouse_exon2, expression_df, coord_to_exon)
        
        if exon1_expr is None or exon2_expr is None:
            continue
        
        # Calculate EEI expression pattern
        eei_expr = (exon1_expr + exon2_expr) / 2  # Simple average
        
        # Test correlation with response groups
        result = {
            'mouse_exon1': mouse_exon1,
            'mouse_exon2': mouse_exon2,
            'human_exon1': eei.get('human_exon1', ''),
            'human_exon2': eei.get('human_exon2', '')
        }
        
        # Analyze by response group
        for group in ['good_responder', 'moderate_responder', 'poor_control']:
            group_samples = metadata_df[metadata_df['response_group'] == group]['sample_id']
            if len(group_samples) > 0:
                # Use pandas indexing to get common samples
                common_samples = group_samples[group_samples.isin(eei_expr.index)]
                group_expr = eei_expr[common_samples]
                result[f'{group}_mean'] = group_expr.mean()
                result[f'{group}_std'] = group_expr.std()
        
        # Calculate response score
        if 'good_responder_mean' in result and 'poor_control_mean' in result:
            result['response_score'] = result['good_responder_mean'] - result['poor_control_mean']
            
            # Statistical test
            good_samples = metadata_df[metadata_df['response_group'] == 'good_responder']['sample_id']
            poor_samples = metadata_df[metadata_df['response_group'] == 'poor_control']['sample_id']
            
            good_common = good_samples[good_samples.isin(eei_expr.index)]
            poor_common = poor_samples[poor_samples.isin(eei_expr.index)]
            good_expr = eei_expr[good_common]
            poor_expr = eei_expr[poor_common]
            
            if len(good_expr) > 1 and len(poor_expr) > 1:
                t_stat, p_val = stats.ttest_ind(good_expr, poor_expr)
                result['p_value'] = p_val
                result['t_statistic'] = t_stat
        
        results.append(result)
    
    results_df = pd.DataFrame(results)

    # Apply Benjamini-Hochberg FDR correction if any p-values are present
    if 'p_value' in results_df.columns and results_df['p_value'].notna().any():
        pvals = results_df['p_value'].values
        # Benjamini-Hochberg procedure
        n = np.sum(~np.isnan(pvals))
        ranks = np.argsort(np.argsort(np.where(np.isnan(pvals), np.inf, pvals))) + 1
        with np.errstate(invalid='ignore'):
            bh = pvals * n / ranks
        # Ensure monotonicity
        order = np.argsort(np.where(np.isnan(pvals), np.inf, pvals))
        bh_ordered = bh[order]
        for i in range(len(bh_ordered) - 2, -1, -1):
            if not np.isnan(bh_ordered[i]):
                bh_ordered[i] = np.nanmin([bh_ordered[i], bh_ordered[i + 1]])
        bh_adj = np.empty_like(pvals, dtype=float)
        bh_adj[:] = np.nan
        bh_adj[order] = bh_ordered
        # Cap at 1.0
        with np.errstate(invalid='ignore'):
            bh_adj = np.minimum(bh_adj, 1.0)
        results_df['p_adj_fdr_bh'] = bh_adj

    return results_df

def identify_time_dependent_eeis(eei_df, expression_df, metadata_df, coord_to_exon):
    """
    Identify EEIs that change over treatment time course.
    The paper shows early changes predict later response.
    """
    
    # Focus on samples with timepoint data
    time_samples = metadata_df[metadata_df['timepoint'].notna()]
    
    # Define early vs late timepoints
    early_timepoints = ['day1', 'day2', 'day3']
    late_timepoints = ['day6', 'day12', '2_cycles', '4_cycles', '7_cycles']
    
    time_dependent_eeis = []
    
    for idx, eei in eei_df.iterrows():
        mouse_exon1 = eei['mouse_exon1']
        mouse_exon2 = eei['mouse_exon2']
        
        # Get expression
        exon1_expr = find_expression_with_mapping(mouse_exon1, expression_df, coord_to_exon)
        exon2_expr = find_expression_with_mapping(mouse_exon2, expression_df, coord_to_exon)
        
        if exon1_expr is None or exon2_expr is None:
            continue
        
        # Calculate EEI presence
        eei_present = (exon1_expr > 1.0) & (exon2_expr > 1.0)  # CPM threshold
        
        # Analyze temporal patterns
        early_samples = time_samples[time_samples['timepoint'].isin(early_timepoints)]['sample_id']
        late_samples = time_samples[time_samples['timepoint'].isin(late_timepoints)]['sample_id']
        
        if len(early_samples) > 0 and len(late_samples) > 0:
            early_common = early_samples[early_samples.isin(eei_present.index)]
            late_common = late_samples[late_samples.isin(eei_present.index)]
            early_presence = eei_present[early_common].mean()
            late_presence = eei_present[late_common].mean()
            
            temporal_change = late_presence - early_presence
            
            if abs(temporal_change) > 0.3:  # Significant change threshold
                time_dependent_eeis.append({
                    'mouse_exon1': mouse_exon1,
                    'mouse_exon2': mouse_exon2,
                    'early_presence': early_presence,
                    'late_presence': late_presence,
                    'temporal_change': temporal_change,
                    'pattern': 'gained' if temporal_change > 0 else 'lost'
                })
    
    return pd.DataFrame(time_dependent_eeis)

def create_pseudo_survival_analysis(eei_response_df, metadata_df):
    """
    Create survival-like curves using response groups as proxy.
    """
    
    # Group EEIs by their association with good response
    p_col = 'p_adj_fdr_bh' if 'p_adj_fdr_bh' in eei_response_df.columns else 'p_value'
    if p_col not in eei_response_df.columns:
        return None
    significant_eeis = eei_response_df[eei_response_df[p_col] < 0.05].copy()
    
    if len(significant_eeis) == 0:
        print("No significant EEIs found for survival analysis")
        return None
    
    # Sort by response score
    significant_eeis = significant_eeis.sort_values('response_score', ascending=False)
    
    # Top positive (associated with good response) and negative (associated with poor response)
    top_positive = significant_eeis.head(10)
    top_negative = significant_eeis.tail(10)
    
    # Create risk groups based on EEI expression patterns
    def calculate_risk_score(sample_expr, positive_eeis, negative_eeis):
        """Calculate risk based on EEI expression patterns"""
        pos_score = 0
        neg_score = 0
        
        # Implementation would calculate aggregate score
        # based on expression of positive and negative EEIs
        
        return pos_score - neg_score
    
    return {
        'top_positive_eeis': top_positive,
        'top_negative_eeis': top_negative,
        'significant_count': len(significant_eeis)
    }   