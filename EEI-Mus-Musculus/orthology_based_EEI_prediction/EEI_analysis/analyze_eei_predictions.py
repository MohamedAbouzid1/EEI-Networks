#!/usr/bin/env python3
"""
EEI Prediction Analysis Pipeline

This script runs a complete analysis of Exon-Exon Interaction (EEI) predictions
based on homology data between mouse and human.
"""

import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import time
from datetime import datetime

# Import our modules - adjust imports based on your file organization
import sys
sys.path.append('.')  # Add current directory to path

# Try importing our modules
try:
    from preprocess_data import *
    from eei_analysis import EEIAnalyzer
    from evaluate_predictions import run_evaluation
    from conservation_analysis import run_conservation_analysis
except ImportError as e:
    print(f"Error importing modules: {e}")
    print("Make sure all required script files are in the correct location.")
    sys.exit(1)

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='EEI Prediction Analysis Pipeline')
    
    # Input files
    parser.add_argument('--human-eei', required=True, help='Path to human EEI network file')
    parser.add_argument('--mouse-eei', required=True, help='Path to mouse EEI network file')
    parser.add_argument('--egio', required=True, help='Path to EGIO orthology mapping file')
    parser.add_argument('--predicted', required=True, help='Path to predicted EEIs file')
    parser.add_argument('--predicted-with-threshold', help='Path to predicted EEIs with confidence threshold')
    
    # Mapping files
    parser.add_argument('--human-exon-map', required=True, help='Path to human exon ID to coordinate mapping file')
    parser.add_argument('--mouse-exon-map', required=True, help='Path to mouse exon ID to coordinate mapping file')
    
    # Parameters
    parser.add_argument('--min-identity', type=float, default=0.8, 
                        help='Minimum identity score for orthology (default: 0.8)')
    parser.add_argument('--min-coverage', type=float, default=5.0, 
                        help='Minimum interface coverage percentage (default: 5.0)')
    parser.add_argument('--analyze-clusters', action='store_true',
                        help='Enable cluster analysis of predictions')
    parser.add_argument('--min-cluster-size', type=int, default=3,
                        help='Minimum size for clusters to analyze')
    
    # Output directories
    parser.add_argument('--output-dir', default='results', 
                        help='Main output directory (default: results)')
    parser.add_argument('--preprocessed-dir', default=None, 
                        help='Directory for preprocessed data (default: {output_dir}/preprocessed)')
    parser.add_argument('--analysis-dir', default=None, 
                        help='Directory for analysis results (default: {output_dir}/analysis)')
    parser.add_argument('--evaluation-dir', default=None, 
                        help='Directory for evaluation results (default: {output_dir}/evaluation)')
    parser.add_argument('--conservation-dir', default=None, 
                        help='Directory for conservation results (default: {output_dir}/conservation)')
    parser.add_argument('--figures-dir', default=None, 
                        help='Directory for figures (default: {output_dir}/figures)')
    
    # Analysis options
    parser.add_argument('--skip-preprocessing', action='store_true', 
                        help='Skip preprocessing step')
    parser.add_argument('--skip-analysis', action='store_true', 
                        help='Skip general analysis step')
    parser.add_argument('--skip-evaluation', action='store_true', 
                        help='Skip prediction evaluation step')
    parser.add_argument('--skip-conservation', action='store_true', 
                        help='Skip conservation analysis step')
    
    return parser.parse_args()

def setup_directories(args):
    """Setup output directories"""
    # Create main output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Set default subdirectories if not specified
    if args.preprocessed_dir is None:
        args.preprocessed_dir = os.path.join(args.output_dir, 'preprocessed')
    
    if args.analysis_dir is None:
        args.analysis_dir = os.path.join(args.output_dir, 'analysis')
    
    if args.evaluation_dir is None:
        args.evaluation_dir = os.path.join(args.output_dir, 'evaluation')
    
    if args.conservation_dir is None:
        args.conservation_dir = os.path.join(args.output_dir, 'conservation')
    
    if args.figures_dir is None:
        args.figures_dir = os.path.join(args.output_dir, 'figures')
    
    # Create subdirectories
    os.makedirs(args.preprocessed_dir, exist_ok=True)
    os.makedirs(args.analysis_dir, exist_ok=True)
    os.makedirs(args.evaluation_dir, exist_ok=True)
    os.makedirs(args.conservation_dir, exist_ok=True)
    os.makedirs(args.figures_dir, exist_ok=True)
    
    return args

def run_pipeline(args):
    """Run the complete analysis pipeline"""
    start_time = time.time()
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    print(f"=== EEI Prediction Analysis Pipeline ===")
    print(f"Started at: {timestamp}")
    print(f"Output directory: {args.output_dir}")
    print(f"Parameters: min_identity={args.min_identity}, min_coverage={args.min_coverage}")
    if args.analyze_clusters:
        print(f"Cluster analysis enabled with min_cluster_size={args.min_cluster_size}")
    if args.predicted_with_threshold:
        print(f"Comparison with thresholded predictions enabled")
    print("=" * 40)
    
    # Setup directories
    args = setup_directories(args)
    
    # Step 1: Load exon ID mappings
    print("\n--- Loading Exon ID Mappings ---")
    mappings = load_exon_mappings(args.human_exon_map, args.mouse_exon_map)
    
    # Step 2: Preprocess data
    if not args.skip_preprocessing:
        print("\n--- Step 2: Preprocessing Data ---")
        # Load raw data
        human_eei_df = pd.read_csv(args.human_eei, sep='\t')
        mouse_eei_df = pd.read_csv(args.mouse_eei, sep='\t')
        egio_df = pd.read_csv(args.egio, sep='\t')
        predicted_df = pd.read_csv(args.predicted, sep='\t') if args.predicted else None
        
        # Load comparison data if provided
        predicted_with_threshold_df = None
        if args.predicted_with_threshold:
            predicted_with_threshold_df = pd.read_csv(args.predicted_with_threshold, sep='\t')
            print(f"Loaded {len(predicted_with_threshold_df)} predictions with threshold")
        
        print(f"Loaded {len(human_eei_df)} human EEIs, {len(mouse_eei_df)} mouse EEIs, "
              f"and {len(egio_df)} orthology mappings")
        if predicted_df is not None:
            print(f"Loaded {len(predicted_df)} predicted EEIs")
        
        # Convert EEIs to coordinate format
        print("\n--- Converting EEIs to Coordinate Format ---")
        human_eei_coord_df = convert_eei_to_coordinates(
            human_eei_df, mappings['human_id_to_coord'])
        mouse_eei_coord_df = convert_eei_to_coordinates(
            mouse_eei_df, mappings['mouse_id_to_coord'])
        
        # No need to convert EGIO, as it already uses coordinates
        # But we should check if the coordinates match
        print("\n--- Checking Coordinate Compatibility ---")
        human_eei_coords = set(human_eei_coord_df['exon1'].unique()).union(
            set(human_eei_coord_df['exon2'].unique()))
        mouse_eei_coords = set(mouse_eei_coord_df['exon1'].unique()).union(
            set(mouse_eei_coord_df['exon2'].unique()))
        
        egio_human_coords = set(egio_df['hsaPos'].dropna().unique())
        egio_mouse_coords = set(egio_df['musPos'].dropna().unique())
        
        human_overlap = human_eei_coords.intersection(egio_human_coords)
        mouse_overlap = mouse_eei_coords.intersection(egio_mouse_coords)
        
        print(f"Human EEI coordinates overlapping with EGIO: {len(human_overlap)} "
              f"({len(human_overlap)/len(human_eei_coords)*100:.1f}% of human EEI coords)")
        print(f"Mouse EEI coordinates overlapping with EGIO: {len(mouse_overlap)} "
              f"({len(mouse_overlap)/len(mouse_eei_coords)*100:.1f}% of mouse EEI coords)")
        
        if len(human_overlap) == 0 or len(mouse_overlap) == 0:
            print("\nWARNING: No coordinate overlap found. Check if the coordinate formats match exactly.")
            print("Sample human EEI coords:", list(human_eei_coords)[:3])
            print("Sample EGIO human coords:", list(egio_human_coords)[:3])
            print("Sample mouse EEI coords:", list(mouse_eei_coords)[:3])
            print("Sample EGIO mouse coords:", list(egio_mouse_coords)[:3])
        
        # Convert predicted EEIs if needed
        if predicted_df is not None:
            predicted_coord_df = convert_eei_to_coordinates(
                predicted_df, mappings['human_id_to_coord'])
        else:
            predicted_coord_df = None
            
        # Convert comparison predictions if provided
        if predicted_with_threshold_df is not None:
            predicted_with_threshold_coord_df = convert_eei_to_coordinates(
                predicted_with_threshold_df, mappings['human_id_to_coord'])
        else:
            predicted_with_threshold_coord_df = None
        
        # Save coordinate-converted files
        os.makedirs(args.preprocessed_dir, exist_ok=True)
        human_eei_coord_file = os.path.join(args.preprocessed_dir, 'human_eei_coord.tsv')
        mouse_eei_coord_file = os.path.join(args.preprocessed_dir, 'mouse_eei_coord.tsv')
        egio_file = args.egio  # EGIO already uses coordinates
        predicted_coord_file = os.path.join(args.preprocessed_dir, 'predicted_coord.tsv') if predicted_df is not None else None
        
        human_eei_coord_df.to_csv(human_eei_coord_file, sep='\t', index=False)
        mouse_eei_coord_df.to_csv(mouse_eei_coord_file, sep='\t', index=False)
        if predicted_coord_df is not None:
            predicted_coord_df.to_csv(predicted_coord_file, sep='\t', index=False)
            
        # Save comparison predictions if provided
        predicted_with_threshold_coord_file = None
        if predicted_with_threshold_coord_df is not None:
            predicted_with_threshold_coord_file = os.path.join(args.preprocessed_dir, 'predicted_with_threshold_coord.tsv')
            predicted_with_threshold_coord_df.to_csv(predicted_with_threshold_coord_file, sep='\t', index=False)
        
        print(f"Saved coordinate-converted files to {args.preprocessed_dir}")
        
        # Use the coordinate-converted files for analysis
        human_eei_file = human_eei_coord_file
        mouse_eei_file = mouse_eei_coord_file
        predicted_file = predicted_coord_file
        predicted_with_threshold_file = predicted_with_threshold_coord_file
    else:
        print("\n--- Skipping Preprocessing ---")
        # Use original input files
        human_eei_file = args.human_eei
        mouse_eei_file = args.mouse_eei
        egio_file = args.egio
        predicted_file = args.predicted
        predicted_with_threshold_file = args.predicted_with_threshold
    
    # Step 2: General analysis
    if not args.skip_analysis:
        print("\n--- Step 2: General EEI Analysis ---")
        analyzer = EEIAnalyzer(
            human_eei_file=human_eei_file,
            mouse_eei_file=mouse_eei_file,
            egio_file=egio_file,
            predicted_eei_file=predicted_file
        )
        
        analyzer.load_data()
        analyzer.compute_basic_stats()
        analyzer.analyze_prediction_overlap()
        analyzer.analyze_conservation_patterns()
        
        # Use improved prediction quality analysis
        if hasattr(analyzer, 'improved_analyze_prediction_quality'):
            analyzer.improved_analyze_prediction_quality()
        else:
            # Fall back to regular analysis if improved version not available
            analyzer.analyze_prediction_quality()
            
        analyzer.create_visualizations(output_dir=args.figures_dir)
        
        # Save statistics to file
        stats_file = os.path.join(args.analysis_dir, 'analysis_stats.txt')
        with open(stats_file, 'w') as f:
            f.write("=== EEI Analysis Statistics ===\n\n")
            for key, value in analyzer.stats.items():
                if isinstance(value, dict):
                    f.write(f"{key}:\n")
                    for subkey, subvalue in value.items():
                        f.write(f"  {subkey}: {subvalue}\n")
                else:
                    f.write(f"{key}: {value}\n")
        
        print(f"Analysis statistics saved to {stats_file}")
    else:
        print("\n--- Skipping General Analysis ---")
    
    # Step 3: Evaluation
    if not args.skip_evaluation:
        print("\n--- Step 3: Prediction Evaluation ---")
        
        # Check if we need to pass predicted_with_threshold_file
        eval_kwargs = {
            'predicted_file': predicted_file,
            'human_eei_file': human_eei_file,
            'mouse_eei_file': mouse_eei_file,
            'egio_file': egio_file,
            'output_dir': args.evaluation_dir
        }
        
        # Add optional parameters if available
        if predicted_with_threshold_file:
            eval_kwargs['predicted_with_threshold_file'] = predicted_with_threshold_file
            
        if args.analyze_clusters:
            eval_kwargs['analyze_clusters'] = True
            eval_kwargs['min_cluster_size'] = args.min_cluster_size
        
        # Run evaluation with all applicable parameters
        evaluation_result = run_evaluation(**eval_kwargs)
        
        print(f"Evaluation completed and saved to {args.evaluation_dir}")
    else:
        print("\n--- Skipping Evaluation ---")
    
    # Step 4: Conservation analysis
    if not args.skip_conservation:
        print("\n--- Step 4: Conservation Analysis ---")
        conservation_result = run_conservation_analysis(
            human_eei_file=human_eei_file,
            mouse_eei_file=mouse_eei_file,
            egio_file=egio_file,
            predicted_file=predicted_file,
            output_dir=args.conservation_dir
        )
        
        print(f"Conservation analysis completed and saved to {args.conservation_dir}")
    else:
        print("\n--- Skipping Conservation Analysis ---")
    
    # Final timing
    end_time = time.time()
    elapsed_time = end_time - start_time
    print("\n" + "=" * 40)
    print(f"Pipeline completed in {elapsed_time:.2f} seconds")
    print(f"Results saved to {args.output_dir}")
    print("=" * 40)

if __name__ == "__main__":
    args = parse_args()
    run_pipeline(args)