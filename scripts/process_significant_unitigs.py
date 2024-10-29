import pandas as pd
import subprocess
import os
import argparse
from tqdm import tqdm
import warnings

def warn_if_file_exists(filepath):
    if os.path.exists(filepath):
        warnings.warn(f"Warning: File {filepath} already exists and will be overwritten.", UserWarning)

def process_significant_unitigs(results_dir, output_folder, input_matrix):
    comparisons = ['AvsI', 'AvsO', 'OvsI']
    ref_sequence = '/pasteur/appa/scratch/cduitama/Damselflies/unitig_matrix_pa/unitigs'
    seqkit_path = '/pasteur/appa/homes/cduitama/venvs/Damselflies/bin/seqkit'

    print(f"Ensuring output folder exists: {output_folder}")
    os.makedirs(output_folder, exist_ok=True)

    print(f"Reading input matrix: {input_matrix}")
    unitigs_df = pd.read_csv(input_matrix, usecols=['Unitigs_id'])
    print(f"Read {len(unitigs_df)} unitig IDs from input matrix")

    for comparison in tqdm(comparisons, desc="Processing comparisons"):
        print(f"\nProcessing comparison: {comparison}")
        
        results_file = f'{results_dir}/{comparison}_chi2_bin_thresh_0.00.csv'
        print(f"Reading results file: {results_file}")
        results_df = pd.read_csv(results_file)
        print(f"Read {len(results_df)} results")
        
        print("Adding Unitigs_id to results")
        results_df['Unitigs_id'] = unitigs_df['Unitigs_id']
        
        print("Filtering significant unitigs")
        significant_df = results_df[results_df['adjusted_p_value'] < 0.05]
        print(f"Found {len(significant_df)} significant unitigs")
        
        significant_ids_file = f'{output_folder}/significant_unitig_ids_{comparison}.txt'
        warn_if_file_exists(significant_ids_file)
        print(f"Saving significant Unitigs_id to: {significant_ids_file}")
        significant_df['Unitigs_id'].to_csv(significant_ids_file, index=False, header=False)
        
        output_fasta = f'{output_folder}/significant_unitigs_{comparison}.fasta.gz'
        warn_if_file_exists(output_fasta)
        print(f"Extracting sequences using seqkit to: {output_fasta}")
        seqkit_command = f'{seqkit_path} grep -f {significant_ids_file} {ref_sequence} -o {output_fasta}'
        subprocess.run(seqkit_command, shell=True, check=True)
        
        stats_file = f'{output_folder}/significant_unitigs_stats_{comparison}.txt'
        warn_if_file_exists(stats_file)
        print(f"Running seqkit stats and saving to: {stats_file}")
        stats_command = f'{seqkit_path} stats {output_fasta} > {stats_file}'
        subprocess.run(stats_command, shell=True, check=True)

    print("\nAll comparisons processed successfully.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process significant unitigs from chi-square test results.")
    parser.add_argument("--results_dir", required=True, help="Directory containing the chi-square test results")
    parser.add_argument("--output_folder", required=True, help="Folder to store the output files")
    parser.add_argument("--input_matrix", required=True, help="Path to the input matrix CSV file")
    
    args = parser.parse_args()
    
    print("Starting significant unitigs processing")
    
    process_significant_unitigs(args.results_dir, args.output_folder, args.input_matrix)
    
    print("Script execution completed.")