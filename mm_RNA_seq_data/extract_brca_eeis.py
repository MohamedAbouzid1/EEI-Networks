#!/usr/bin/env python3
"""
Script to extract BRCA-specific exon-exon interactions from 3_NETHIGH_CRPES dataset
Using only standard Python libraries
"""

import os
import glob
from collections import defaultdict

def analyze_brca_files():
    """Analyze BRCA-specific files in the 3_NETHIGH_CRPES dataset"""
    
    base_path = "data/Final_survival_filt/3_NETHIGH_CRPES/threshold_0.5"
    
    # Find all BRCA files
    brca_files = glob.glob(os.path.join(base_path, "BRCA_*.txt"))
    
    print("=== BRCA-Specific Exon-Exon Interactions Analysis ===\n")
    print(f"Found {len(brca_files)} BRCA-specific files:\n")
    
    for file_path in sorted(brca_files):
        filename = os.path.basename(file_path)
        file_size = os.path.getsize(file_path)
        
        # Read first few lines to understand structure
        try:
            with open(file_path, 'r') as f:
                lines = f.readlines()
                header = lines[0].strip()
                num_lines = len(lines)
                
                print(f"File: {filename}")
                print(f"Size: {file_size:,} bytes")
                print(f"Lines: {num_lines:,}")
                print(f"Header: {header}")
                
                # Show first few data rows
                if num_lines > 1:
                    print("First few data rows:")
                    for i in range(1, min(6, num_lines)):
                        print(f"  {lines[i].strip()}")
                print("-" * 50)
                
        except Exception as e:
            print(f"Error reading {filename}: {e}")
            print("-" * 50)

def extract_brca_summary():
    """Create a summary of BRCA exon-exon interactions"""
    
    print("\n=== BRCA Exon-Exon Interactions Summary ===\n")
    
    # Define file types and their descriptions
    file_types = {
        "BRCA_Gained_.txt": "BRCA Gained Interactions (all patients)",
        "BRCA_Lost_.txt": "BRCA Lost Interactions (all patients)", 
        "BRCA_Gained_Surv_.txt": "BRCA Gained Interactions with Survival Data",
        "BRCA_Lost_Surv_.txt": "BRCA Lost Interactions with Survival Data"
    }
    
    base_path = "data/Final_survival_filt/3_NETHIGH_CRPES/threshold_0.5"
    
    for filename, description in file_types.items():
        file_path = os.path.join(base_path, filename)
        
        if os.path.exists(file_path):
            try:
                # Read the file
                with open(file_path, 'r') as f:
                    lines = f.readlines()
                
                if len(lines) > 1:
                    header = lines[0].strip().split('\t')
                    data_lines = lines[1:]
                    
                    print(f"📁 {filename}")
                    print(f"   Description: {description}")
                    print(f"   Total interactions: {len(data_lines):,}")
                    print(f"   Columns: {header}")
                    
                    # Analyze weights if present
                    if 'weight' in header:
                        weight_idx = header.index('weight')
                        weights = []
                        for line in data_lines:
                            parts = line.strip().split('\t')
                            if len(parts) > weight_idx:
                                try:
                                    weights.append(float(parts[weight_idx]))
                                except:
                                    pass
                        
                        if weights:
                            print(f"   Weight range: {min(weights):.3f} - {max(weights):.3f}")
                    
                    # Analyze patients if present
                    if 'patient' in header:
                        patient_idx = header.index('patient')
                        patients = set()
                        for line in data_lines:
                            parts = line.strip().split('\t')
                            if len(parts) > patient_idx:
                                patients.add(parts[patient_idx])
                        
                        print(f"   Unique patients: {len(patients)}")
                    
                    # Analyze p-values if present
                    if 'pval' in header:
                        pval_idx = header.index('pval')
                        significant = 0
                        for line in data_lines:
                            parts = line.strip().split('\t')
                            if len(parts) > pval_idx:
                                try:
                                    pval = float(parts[pval_idx])
                                    if pval < 0.05:
                                        significant += 1
                                except:
                                    pass
                        
                        print(f"   Significant interactions (p<0.05): {significant:,}")
                    
                    print()
                
            except Exception as e:
                print(f"Error processing {filename}: {e}")
                print()

def create_brca_combined_dataset():
    """Create a combined BRCA dataset with all interactions"""
    
    print("=== Creating Combined BRCA Dataset ===\n")
    
    base_path = "data/Final_survival_filt/3_NETHIGH_CRPES/threshold_0.5"
    
    # Read all BRCA files
    brca_data = []
    total_interactions = 0
    
    for file_type in ["Gained", "Lost"]:
        # Regular file
        file_path = os.path.join(base_path, f"BRCA_{file_type}_.txt")
        if os.path.exists(file_path):
            with open(file_path, 'r') as f:
                lines = f.readlines()
            
            if len(lines) > 1:
                header = lines[0].strip().split('\t')
                data_lines = lines[1:]
                
                for line in data_lines:
                    parts = line.strip().split('\t')
                    if len(parts) >= len(header):
                        row_data = dict(zip(header, parts))
                        row_data['interaction_type'] = file_type
                        row_data['has_survival'] = 'False'
                        brca_data.append(row_data)
                
                print(f"Loaded {len(data_lines):,} {file_type} interactions")
                total_interactions += len(data_lines)
        
        # Survival file
        surv_file_path = os.path.join(base_path, f"BRCA_{file_type}_Surv_.txt")
        if os.path.exists(surv_file_path):
            with open(surv_file_path, 'r') as f:
                lines = f.readlines()
            
            if len(lines) > 1:
                header = lines[0].strip().split('\t')
                data_lines = lines[1:]
                
                for line in data_lines:
                    parts = line.strip().split('\t')
                    if len(parts) >= len(header):
                        row_data = dict(zip(header, parts))
                        row_data['interaction_type'] = file_type
                        row_data['has_survival'] = 'True'
                        brca_data.append(row_data)
                
                print(f"Loaded {len(data_lines):,} {file_type} interactions with survival data")
                total_interactions += len(data_lines)
    
    # Count unique exon pairs
    unique_pairs = set()
    for row in brca_data:
        if 'V1' in row and 'V2' in row:
            unique_pairs.add((row['V1'], row['V2']))
    
    print(f"\nCombined dataset:")
    print(f"Total BRCA interactions: {total_interactions:,}")
    print(f"Unique exon pairs: {len(unique_pairs):,}")
    
    # Save combined dataset
    if brca_data:
        output_file = "BRCA_combined_EEIs.txt"
        with open(output_file, 'w') as f:
            # Write header
            if brca_data:
                header = list(brca_data[0].keys())
                f.write('\t'.join(header) + '\n')
                
                # Write data
                for row in brca_data:
                    f.write('\t'.join(str(row.get(col, '')) for col in header) + '\n')
        
        print(f"Saved combined dataset to: {output_file}")
        
        return brca_data
    
    return None

if __name__ == "__main__":
    # Analyze BRCA files
    analyze_brca_files()
    
    # Create summary
    extract_brca_summary()
    
    # Create combined dataset
    combined_data = create_brca_combined_dataset()
    
    print("\n=== Extraction Complete ===")
    print("The BRCA-specific exon-exon interactions have been extracted and analyzed.")
    print("You can find the combined dataset in 'BRCA_combined_EEIs.txt'")
