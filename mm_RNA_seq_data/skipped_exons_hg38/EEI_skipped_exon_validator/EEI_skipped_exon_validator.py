#!/usr/bin/env python3

import pandas as pd
import numpy as np
import os
from pathlib import Path
import argparse
from typing import Dict, List, Tuple
import warnings
warnings.filterwarnings('ignore')

class EEISkippedExonValidator:
    """
    Class to handle validation of predicted EEI networks against skipped exons data
    """
    
    def __init__(self, tolerance: int = 10):
        """
        Initialize validator
        
        Parameters:
        - tolerance: Allowed coordinate mismatch in base pairs (default 10)
        """
        self.tolerance = tolerance
        self.stats = {}
        
    def parse_coordinate_string(self, coord_string: str) -> Tuple[str, int, int, str]:
        """
        Parse coordinate format: 'chr1:3069168:3069296:1'
        
        Returns: (chromosome, start, end, strand)
        """
        try:
            parts = coord_string.split(':')
            if len(parts) == 4:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                strand = parts[3]
                return chrom, start, end, strand
            else:
                return None, None, None, None
        except:
            return None, None, None, None
    
    def load_coordinate_mapping(self, coord_file: str) -> Dict:
        """
        Load the exon coordinate mapping file
        Format: exon_id, coord, reverse_lookup
        """
        print(f"Loading coordinate mapping from {coord_file}...")
        coord_df = pd.read_csv(coord_file, sep='\t')
        
        # Create mapping dictionary
        coord_map = {}
        for _, row in coord_df.iterrows():
            coord_map[row['exon_id']] = row['coord']
            # Also map reverse_lookup if different
            if row['reverse_lookup'] != row['exon_id']:
                coord_map[row['reverse_lookup']] = row['coord']
        
        print(f"Loaded {len(coord_map)} exon coordinate mappings")
        return coord_map
    
    def check_coordinate_overlap(self, coord1: str, skip_chrom: str, 
                                skip_start: int, skip_end: int) -> bool:
        """
        Check if EEI exon overlaps with skipped exon
        """
        # Parse the coordinate string
        chrom, start, end, strand = self.parse_coordinate_string(coord1)
        
        if chrom is None:
            return False
        
        # Normalize chromosome names (remove 'chr' if needed for comparison)
        if skip_chrom.startswith('chr') and not chrom.startswith('chr'):
            chrom = 'chr' + chrom
        elif not skip_chrom.startswith('chr') and chrom.startswith('chr'):
            chrom = chrom.replace('chr', '')
        
        if chrom != skip_chrom:
            return False
        
        # Check for overlap with tolerance
        return (start - self.tolerance <= skip_end and 
                end + self.tolerance >= skip_start)
    
    def find_matches(self, eei_df: pd.DataFrame, skipped_df: pd.DataFrame,
                    coord_map: Dict) -> pd.DataFrame:
        """
        Find matches between EEI exons and skipped exons using coordinates
        """
        matches = []
        total_eeis = len(eei_df)
        print(f"Searching for matches in {total_eeis} EEIs...")
        
        # Progress reporting
        report_interval = max(1, total_eeis // 10)
        
        for idx, eei_row in eei_df.iterrows():
            if idx % report_interval == 0:
                print(f"  Processed {idx}/{total_eeis} EEIs...")
            
            # Get coordinates for both exons
            exon1_coord = coord_map.get(eei_row['exon1'])
            exon2_coord = coord_map.get(eei_row['exon2'])
            
            # Skip if we don't have coordinates
            if not exon1_coord and not exon2_coord:
                continue
            
            # Check against all skipped exons
            for _, skip_row in skipped_df.iterrows():
                # Check exon1
                if exon1_coord and self.check_coordinate_overlap(
                    exon1_coord, skip_row['chrom'], 
                    skip_row['A_start'], skip_row['A_end']
                ):
                    match_info = {
                        'eei_exon1': eei_row['exon1'],
                        'eei_exon2': eei_row['exon2'],
                        'eei_exon1_coord': exon1_coord,
                        'eei_exon2_coord': exon2_coord,
                        'matched_exon': 'exon1',
                        'skipped_event_id': skip_row['event_id'],
                        'skipped_coords': f"{skip_row['chrom']}:{skip_row['A_start']}-{skip_row['A_end']}",
                        'skipped_length': skip_row['A_len'],
                        'confidence': eei_row.get('confidence', None),
                        'identity1': eei_row.get('identity1', None),
                        'identity2': eei_row.get('identity2', None)
                    }
                    matches.append(match_info)
                
                # Check exon2
                if exon2_coord and self.check_coordinate_overlap(
                    exon2_coord, skip_row['chrom'], 
                    skip_row['A_start'], skip_row['A_end']
                ):
                    match_info = {
                        'eei_exon1': eei_row['exon1'],
                        'eei_exon2': eei_row['exon2'],
                        'eei_exon1_coord': exon1_coord,
                        'eei_exon2_coord': exon2_coord,
                        'matched_exon': 'exon2',
                        'skipped_event_id': skip_row['event_id'],
                        'skipped_coords': f"{skip_row['chrom']}:{skip_row['A_start']}-{skip_row['A_end']}",
                        'skipped_length': skip_row['A_len'],
                        'confidence': eei_row.get('confidence', None),
                        'identity1': eei_row.get('identity1', None),
                        'identity2': eei_row.get('identity2', None)
                    }
                    matches.append(match_info)
        
        print(f"  Found {len(matches)} total matches")
        return pd.DataFrame(matches)
    
    def analyze_psi_values(self, matches_df: pd.DataFrame, 
                          psi_dir: str) -> pd.DataFrame:
        """
        Load and analyze PSI values for matched exons
        """
        psi_summaries = []
        unique_events = matches_df['skipped_event_id'].unique()
        print(f"Analyzing PSI values for {len(unique_events)} unique events...")
        
        for event_id in unique_events:
            psi_file = Path(psi_dir) / f"counts_psi__{event_id}.tsv"
            
            if psi_file.exists():
                try:
                    psi_df = pd.read_csv(psi_file, sep='\t')
                    
                    # Calculate statistics
                    psi_stats = {
                        'event_id': event_id,
                        'mean_psi': psi_df['PSI'].mean(),
                        'std_psi': psi_df['PSI'].std(),
                        'min_psi': psi_df['PSI'].min(),
                        'max_psi': psi_df['PSI'].max(),
                        'median_psi': psi_df['PSI'].median(),
                        'n_samples': len(psi_df),
                        'n_tissues': psi_df['tissue'].nunique() if 'tissue' in psi_df.columns else 0,
                        'high_variance': psi_df['PSI'].std() > 0.2,
                        'constitutive': psi_df['PSI'].mean() > 0.9 and psi_df['PSI'].std() < 0.1,
                        'cassette': psi_df['PSI'].std() > 0.3
                    }
                    
                    # Get tissue-specific stats if available
                    if 'tissue' in psi_df.columns:
                        tissue_stats = psi_df.groupby('tissue')['PSI'].agg(['mean', 'std'])
                        psi_stats['max_tissue_diff'] = tissue_stats['mean'].max() - tissue_stats['mean'].min()
                    
                    psi_summaries.append(psi_stats)
                    
                except Exception as e:
                    print(f"Error reading PSI file for {event_id}: {e}")
        
        return pd.DataFrame(psi_summaries)
    
    def calculate_statistics(self, eei_df: pd.DataFrame, 
                            matches_df: pd.DataFrame,
                            method_name: str = "Unknown") -> Dict:
        """
        Calculate comprehensive validation statistics
        """
        # Remove duplicate EEI pairs
        unique_eeis_matched = matches_df[['eei_exon1', 'eei_exon2']].drop_duplicates()
        
        # Count unique exons
        all_exons = set(eei_df['exon1']) | set(eei_df['exon2'])
        matched_exons = set(matches_df['eei_exon1']) | set(matches_df['eei_exon2'])
        
        stats = {
            'method': method_name,
            'total_eeis': len(eei_df),
            'total_unique_exons': len(all_exons),
            'eeis_with_skipped_exons': len(unique_eeis_matched),
            'unique_exons_matched': len(matched_exons),
            'percentage_eeis_validated': (len(unique_eeis_matched) / len(eei_df) * 100) if len(eei_df) > 0 else 0,
            'percentage_exons_matched': (len(matched_exons) / len(all_exons) * 100) if len(all_exons) > 0 else 0,
            'total_matches': len(matches_df),
            'unique_skipped_events': matches_df['skipped_event_id'].nunique() if len(matches_df) > 0 else 0,
            'exon1_matches': len(matches_df[matches_df['matched_exon'] == 'exon1']) if len(matches_df) > 0 else 0,
            'exon2_matches': len(matches_df[matches_df['matched_exon'] == 'exon2']) if len(matches_df) > 0 else 0,
        }
        
        # Add confidence score statistics if available
        if 'confidence' in eei_df.columns:
            matched_eeis_df = eei_df[
                (eei_df[['exon1', 'exon2']].apply(tuple, axis=1).isin(
                    unique_eeis_matched[['eei_exon1', 'eei_exon2']].apply(tuple, axis=1)))
            ]
            if len(matched_eeis_df) > 0:
                stats['mean_confidence_matched'] = matched_eeis_df['confidence'].mean()
                stats['mean_confidence_all'] = eei_df['confidence'].mean()
                stats['confidence_diff'] = stats['mean_confidence_matched'] - stats['mean_confidence_all']
        
        return stats

def run_single_method_validation(eei_file: str, skipped_file: str, 
                                coord_file: str, psi_dir: str,
                                method_name: str, output_dir: Path,
                                tolerance: int = 10):
    """
    Run validation for a single EEI prediction method
    """
    print(f"\n{'='*60}")
    print(f"Processing {method_name} method")
    print('='*60)
    
    # Initialize validator
    validator = EEISkippedExonValidator(tolerance=tolerance)
    
    # Load coordinate mapping
    coord_map = validator.load_coordinate_mapping(coord_file)
    
    # Load data
    print(f"Loading EEI data...")
    eei_df = pd.read_csv(eei_file, sep='\t')
    print(f"  Loaded {len(eei_df)} EEIs")
    
    print(f"Loading skipped exons data...")
    skipped_df = pd.read_csv(skipped_file, sep='\t')
    print(f"  Loaded {len(skipped_df)} skipped exon events")
    
    # Find matches
    matches_df = validator.find_matches(eei_df, skipped_df, coord_map)
    
    # Analyze PSI if available
    if psi_dir and len(matches_df) > 0:
        print("\nAnalyzing PSI values...")
        psi_summary = validator.analyze_psi_values(matches_df, psi_dir)
        
        if len(psi_summary) > 0:
            matches_df = matches_df.merge(
                psi_summary,
                left_on='skipped_event_id',
                right_on='event_id',
                how='left'
            )
    
    # Calculate statistics
    stats = validator.calculate_statistics(eei_df, matches_df, method_name)
    
    # Save results
    output_file = output_dir / f"{method_name.lower()}_validation_results.tsv"
    matches_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nSaved results to {output_file}")
    
    # Print summary
    print(f"\n{method_name} Validation Summary:")
    print(f"  Total EEIs: {stats['total_eeis']}")
    print(f"  EEIs with skipped exons: {stats['eeis_with_skipped_exons']} ({stats['percentage_eeis_validated']:.1f}%)")
    print(f"  Unique exons matched: {stats['unique_exons_matched']} ({stats['percentage_exons_matched']:.1f}%)")
    print(f"  Unique skipped events: {stats['unique_skipped_events']}")
    
    return matches_df, stats

def main():
    """
    Main function to run the validation pipeline
    """
    # Set up paths for your specific environment
    base_dir = Path("/cosybio_project/mabouzid/EEI_networks/EEI-Conservation-main")
    
    # Input files
    eei_files = {
        'PISA': base_dir / "orthology_based_EEI_prediction/results_PISA_Based/predicted_human_eeis.tsv",
        'EPPIC': base_dir / "orthology_based_EEI_prediction/results_EPPIC_Based/predicted_human_eeis_fixed_no_conf_threshold_all_iden_0.tsv",
        'Contact': base_dir / "orthology_based_EEI_prediction/results_CONTACT_BASED/predicted_human_eeis.tsv"
    }
    
    skipped_file = base_dir / "mm_RNA_seq_data/skipped_exons_hg38/skipped_exon_coords.tsv"
    coord_file = base_dir / "orthology_based_EEI_prediction/results_CONTACT_BASED/exon_coord_mapping_human.tsv"
    psi_dir = base_dir / "mm_RNA_seq_data/skipped_exons_hg38/ex_psi"
    output_dir = base_dir / "mm_RNA_seq_data/skipped_exons_hg38/validation_results"
    
    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("="*60)
    print("EEI-SKIPPED EXON VALIDATION PIPELINE")
    print("="*60)
    print(f"\nBase directory: {base_dir}")
    print(f"Output directory: {output_dir}")
    
    # Check if files exist
    if not coord_file.exists():
        print(f"ERROR: Coordinate mapping file not found: {coord_file}")
        return
    
    if not skipped_file.exists():
        print(f"ERROR: Skipped exons file not found: {skipped_file}")
        return
    
    # Process each method
    all_results = {}
    all_stats = []
    
    for method, eei_file in eei_files.items():
        if eei_file.exists():
            matches_df, stats = run_single_method_validation(
                str(eei_file), str(skipped_file), str(coord_file),
                str(psi_dir), method, output_dir, tolerance=10
            )
            all_results[method] = matches_df
            all_stats.append(stats)
        else:
            print(f"\nWarning: {method} file not found: {eei_file}")
    
    # Save combined statistics
    if all_stats:
        stats_df = pd.DataFrame(all_stats)
        stats_file = output_dir / "validation_statistics_summary.tsv"
        stats_df.to_csv(stats_file, sep='\t', index=False)
        print(f"\n{'='*60}")
        print(f"Saved summary statistics to {stats_file}")
        
        # Print comparison
        print("\n" + "="*60)
        print("METHOD COMPARISON")
        print("="*60)
        print(stats_df[['method', 'total_eeis', 'eeis_with_skipped_exons', 
                       'percentage_eeis_validated', 'unique_skipped_events']].to_string(index=False))
    
    # Analyze overlap between methods
    if len(all_results) > 1:
        print("\n" + "="*60)
        print("CROSS-METHOD ANALYSIS")
        print("="*60)
        
        for method1 in all_results:
            for method2 in all_results:
                if method1 < method2:  # Only compare each pair once
                    # Get unique matched EEI pairs
                    eeis1 = set(all_results[method1][['eei_exon1', 'eei_exon2']].apply(tuple, axis=1))
                    eeis2 = set(all_results[method2][['eei_exon1', 'eei_exon2']].apply(tuple, axis=1))
                    
                    overlap = len(eeis1 & eeis2)
                    only1 = len(eeis1 - eeis2)
                    only2 = len(eeis2 - eeis1)
                    
                    print(f"\n{method1} vs {method2}:")
                    print(f"  Shared validated EEIs: {overlap}")
                    print(f"  {method1}-only: {only1}")
                    print(f"  {method2}-only: {only2}")
    
    print("\n" + "="*60)
    print("VALIDATION COMPLETED SUCCESSFULLY!")
    print("="*60)

if __name__ == "__main__":
    main()