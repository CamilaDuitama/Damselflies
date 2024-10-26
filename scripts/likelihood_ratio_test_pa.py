import os
import argparse
import pandas as pd
import dask.dataframe as dd
from dask.distributed import Client, LocalCluster
from scipy.stats import bernoulli, binomtest, false_discovery_control
import numpy as np
from scipy import stats
from collections import defaultdict

# Constants for metadata and column names files
METADATA_FILE = "/pasteur/appa/scratch/cduitama/Damselflies/metadata/Damselflies_metadata_with_morph.csv"
COLUMN_NAMES_FILE = "/pasteur/appa/scratch/cduitama/Damselflies/metadata/fof_filtered_pa.txt"

def setup_dask_client(n_workers, memory_limit):
    """
    Initialize a Dask Cluster and client for distributed computation.
    
    Args:
    n_workers (int): Number of Dask worker processes.
    memory_limit (str): Memory limit per worker (e.g., '50GB').

    Returns:
    client (Client): Initialized Dask distributed client object.
    """
    cluster = LocalCluster(n_workers=n_workers, threads_per_worker=1, memory_limit=memory_limit)
    client = Client(cluster)
    print("Dask client setup complete.")
    print(f"Dask dashboard available at: {client.dashboard_link}")
    return client

def load_metadata(filepath):
    """
    Load metadata from a CSV file into a Pandas DataFrame.
    
    Args:
    filepath (str): Path to CSV file containing metadata.

    Returns:
    DataFrame: Pandas DataFrame containing loaded metadata.
    """
    metadata = pd.read_csv(filepath)
    print("Metadata loaded.")
    return metadata


def map_runs_to_morph(metadata):
    """
    Map experimental runs to corresponding morph groups.
    
    Args:
    metadata (DataFrame): DataFrame containing 'Run' and 'Morph' columns.

    Returns:
    dict: Dictionary mapping morph groups to lists of runs.
    """
    run_to_morph = dict(zip(metadata['Run'], metadata['Morph']))
    morph_groups = defaultdict(list)
    for run, morph in run_to_morph.items():
        morph_groups[morph].append(run)
    print("Runs mapped to morph groups.")
    return morph_groups

def compute_likelihood_ratio_test_bernoulli(row, groups):
    """
    Compute Bernoulli Likelihood Ratio Test for a given row of binary data.

    Args:
        row (pandas.Series): A row of binary data (values: 0 or 1) from a DataFrame.
        groups (dict): A dictionary mapping group names to the corresponding column indices.

    Returns:
        pandas.Series: A series containing the LRT statistic, p-value and effect size.
        
    This function calculates Bernoulli-based probabilities and likelihoods for binary input data,
    computes the likelihood ratio statistic based on these Bernoulli likelihoods, and the effect size, 
    specifically the absolute difference between success probabilities of two groups.
    """
    # Group data and validation
    group_data = {group: row[samples].dropna().values for group, samples in groups.items() if all(sample in row.index for sample in samples)}
  
    if not group_data or not all(len(data) > 0 for data in group_data.values()) or len(group_data) < 2:
        return pd.Series({'lr_statistic': np.nan, 'p_value': np.nan, 'effect_size': np.nan})

    # Data collection
    y, group_labels = [], []
    for label, values in enumerate(group_data.values()):
        y.extend(values)
        group_labels.extend([label]*len(values))
    y = np.array(y)
    group_labels = np.array(group_labels)
    # Unique groups and their probabilities
    unique_groups = np.unique(group_labels)
    prob_hat = [np.mean(y[group_labels == group]) for group in unique_groups]
    
    # Null probability and log-likelihood
    prob_null = np.mean(y)
    llf_null = np.sum(bernoulli.logpmf(y, prob_null))

    # Full Log-likelihood
    llf_full = sum([np.sum(bernoulli.logpmf(y[group_labels == group], prob)) for group, prob in zip(unique_groups, prob_hat)])

    # LR statistic and p-value
    lr_statistic = -2 * (llf_null - llf_full)
    p_value = binomtest(sum(y), len(y), prob_null, alternative='two-sided').pvalue

    # Effect size as absolute difference of success probabilities
    effect_size = abs(prob_hat[0] - prob_hat[1])

    return pd.Series({'lr_statistic': lr_statistic, 'p_value': p_value, 'effect_size': effect_size})

def perform_analysis(df, morph_groups, comparison):

    if comparison != 'overall':
        desired_groups = [group.strip() for group in comparison.split('vs')]
        groups = {k: morph_groups[k] for k in morph_groups if k in desired_groups}
        if not groups:
            print(f"Error: Specified groups in comparison ({desired_groups}) not found.")
            return None
    else:
        groups = morph_groups

    relevant_columns = [column for group in groups.values() for column in group]
    if 'Unitigs_id' in relevant_columns: relevant_columns.remove('Unitigs_id')  # Ensure unitigs_id is not included

    df_slice = df[relevant_columns]
    df_binarized = df_slice > 0.8
    df_binarized = df_binarized.astype(int)
    #df_binarized=df_slice
    print("Dataframe was sliced and binarized.")
    print(df_binarized.head())
    try:
        meta = pd.Series({
            'lr_statistic': np.float64,
            'p_value': np.float64,
            'effect_size': np.float64
        }, index=['lr_statistic', 'p_value', 'effect_size'])        
        results = df_binarized.apply(lambda row: compute_likelihood_ratio_test_bernoulli(row, groups), axis=1,meta=meta)
    except Exception as e:
        print(f"Error during analysis: {e}")
        return None

    # Compute the adjusted p-values and add them to the DataFrame 
    results_df=results.compute()
    p_values = results_df['p_value'] # compute to create a numpy array
    adjusted_p_values = false_discovery_control(p_values, method='bh')
    results_df['adjusted_p_value'] = adjusted_p_values
    return results_df


def check_columns(file_path):
    # Read only a small sample of the CSV to get the columns and dtypes
    small_df = dd.read_csv(file_path, sample=1000000)  # only sample 1 MB of data
    actual_columns = small_df.columns
    num_columns = len(actual_columns)
    print(f"Sample data check: The CSV file appears to have {num_columns} columns.")
    return num_columns, actual_columns

def main(args):
    """
    Main function to set up the environment, process input, and save outputs.

    Args:
    args (Namespace): Argument namespace containing command line arguments.
    """
    client = setup_dask_client(args.threads, args.memory)
    metadata = load_metadata(METADATA_FILE)
    morph_groups = map_runs_to_morph(metadata)

    # Read column names from COLUMN_NAMES_FILE
    with open(COLUMN_NAMES_FILE, 'r') as f:
        column_names = [line.strip() for line in f]
    column_names = [name.replace("data/unitigs/", "").replace(".unitigs.fa.gz", "") for name in column_names]
    
    # Insert 'Unitigs_id' at the beginning of your list of column names if it's not already included
    if column_names[0] != 'Unitigs_id':
        column_names.insert(0, 'Unitigs_id')
    print("Column names prepared.")
    
    # Efficient column number and names check with Dask
    num_columns, actual_columns = check_columns(args.unitig_matrix)

    print(f"File {args.unitig_matrix} detected to have {num_columns} columns.")
    print(f"Prepared {len(column_names)} column names from the provided names file.")

    if len(column_names) != num_columns:
        raise ValueError("Mismatch between the number of prepared column names and the number of columns detected in the file.")
    
    dtype_dict = {col: float for col in column_names if col != "Unitigs_id"}
    dtype_dict['Unitigs_id'] = str 
    df = dd.read_csv(args.unitig_matrix, header=None, skiprows=1, dtype=dtype_dict,names=column_names)
    # Correctly specifying `names=` parameter according to actual columns in file
    print("Unitig data loaded into Dask DataFrame.")

    results = perform_analysis(df, morph_groups, args.comparison)

    if results is not None:
        # Ensure the output directory exists
        os.makedirs(args.output_dir, exist_ok=True)

        # Create the final output path
        formatted_comparison = args.comparison.strip().replace(' ', '_').replace('/', '').replace('\\', '')
        filename = f"{formatted_comparison}.csv"
        output_file_path = os.path.join(args.output_dir, filename)

        print(f"Attempting to save results to: {output_file_path}")

        # Compute results and save to CSV
        try:
            results.to_csv(output_file_path, index=False)
            print(f"Results have been saved to {output_file_path}")
        except Exception as e:
            print(f"Failed to save results: {e}")
    else:
        print("No results were generated.")
    client.close()
    print("Dask client closed.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform likelihood ratio tests on unitig data.")
    parser.add_argument("-t", "--threads", type=int, required=True, help="Number of threads to use")
    parser.add_argument("-m", "--memory", required=True, help="Amount of memory to use (e.g., '32GB')")
    parser.add_argument("-u", "--unitig_matrix", required=True, help="Path to the unitig matrix file")
    parser.add_argument("-c", "--comparison", required=True,
                    help="Type of comparison to perform, formatted as 'Group1vsGroup2' (e.g., 'AvsB') or 'overall'")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory to save the results")
    args = parser.parse_args()
    main(args)