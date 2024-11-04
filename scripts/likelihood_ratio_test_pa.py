import os
import sys
import argparse
import pandas as pd
import logging
from datetime import datetime
import dask.dataframe as dd
from dask.distributed import Client, LocalCluster
from scipy.stats import bernoulli, binomtest, kruskal, false_discovery_control, chi2, chi2_contingency, norm, poisson
import numpy as np
from scipy import stats
from collections import defaultdict

# Constants for metadata and column names files
METADATA_FILE = "/pasteur/appa/scratch/cduitama/Damselflies/metadata/Damselflies_metadata_with_morph.csv"
COLUMN_NAMES_FILE = "/pasteur/appa/scratch/cduitama/Damselflies/metadata/fof_filtered_pa.txt"

def setup_logging(output_dir):
    """
    Set up logging configuration.

    Args:
    output_dir (str): Directory to save the log file.

    Returns:
    str: Path to the created log file.
    """
    log_filename = f"run_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    log_path = os.path.join(output_dir, log_filename)
    
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        handlers=[
                            logging.FileHandler(log_path),
                            logging.StreamHandler()
                        ])
    return log_path

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
    logging.info("Dask client setup complete.")
    logging.info(f"Dask dashboard available at: {client.dashboard_link}")
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
    logging.info("Metadata loaded.")
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
    logging.info("Runs mapped to morph groups.")
    return morph_groups

def compute_chi_square_test(row, groups):
    """
    Compute Chi-square test for a given row of binary data with two groups.

    Args:
        row (pandas.Series): A row of binary data (values: 0 or 1) from a DataFrame.
        groups (dict): A dictionary mapping group names to the corresponding column indices.

    Returns:
        pandas.Series: A series containing the chi-square statistic, p-value, and effect size (Cramer's V).
    """
    if len(groups) != 2:
        raise ValueError("This function is designed for exactly two groups.")

    group_data = {group: row[samples].dropna().values for group, samples in groups.items() if all(sample in row.index for sample in samples)}
  
    # Check if we have data for both groups
    if len(group_data) < 2 or not all(len(data) > 0 for data in group_data.values()):
        return pd.Series({'chi2_statistic': np.nan, 'p_value': np.nan, 'effect_size': np.nan})

    # Prepare contingency table with Yates' correction
    contingency_table = np.array([
        [np.sum(group_data[group] == 0) + 0.5 for group in groups],
        [np.sum(group_data[group] == 1) + 0.5 for group in groups]
    ])

    # Check if all numbers are identical across both groups
    if np.all(contingency_table == contingency_table[0][0]):
        return pd.Series({'chi2_statistic': np.nan, 'p_value': np.nan, 'effect_size': np.nan})

    try:
        # Compute Chi-square test
        chi2_statistic, p_value, dof, expected = chi2_contingency(contingency_table)

        # Compute effect size (Cramer's V)
        n = np.sum(contingency_table)
        min_dim = min(contingency_table.shape) - 1
        effect_size = np.sqrt(chi2_statistic / (n * min_dim))

    except ValueError as e:
        logging.warning(f"ValueError in Chi-square test: {e}")
        return pd.Series({'chi2_statistic': np.nan, 'p_value': np.nan, 'effect_size': np.nan})

    return pd.Series({'chi2_statistic': chi2_statistic, 'p_value': p_value, 'effect_size': effect_size})

def compute_likelihood_ratio_test_bernoulli(row, groups):
    """
    Compute Bernoulli Likelihood Ratio Test for a given row of binary data.

    Args:
        row (pandas.Series): A row of binary data (values: 0 or 1) from a DataFrame.
        groups (dict): A dictionary mapping group names to the corresponding column indices.

    Returns:
        pandas.Series: A series containing the LRT statistic, p-value and effect size.
    """
    # Group data preparation: Create a dictionary of group data, dropping NaN values
    group_data = {group: row[samples].dropna().values for group, samples in groups.items() if all(sample in row.index for sample in samples)}
  
    # Input validation: Check if there's enough valid data for at least two groups
    if not group_data or not all(len(data) > 0 for data in group_data.values()) or len(group_data) < 2:
        return pd.Series({'lr_statistic': np.nan, 'p_value': np.nan, 'effect_size': np.nan})

    # Data collection: Combine all group data into single arrays
    y, group_labels = [], []
    for label, values in enumerate(group_data.values()):
        y.extend(values)
        group_labels.extend([label]*len(values))
    y = np.array(y)  # Array of all binary outcomes
    group_labels = np.array(group_labels)  # Corresponding group labels

    # Calculate probabilities for each unique group
    unique_groups = np.unique(group_labels)
    prob_hat = [np.mean(y[group_labels == group]) for group in unique_groups]
    
    # Adjust probabilities to avoid 0 and 1 (which can cause issues in log calculations)
    epsilon = 1e-10
    prob_hat = [max(min(p, 1-epsilon), epsilon) for p in prob_hat]
    
    # Calculate overall probability (null hypothesis)
    prob_null = np.mean(y)
    prob_null = max(min(prob_null, 1-epsilon), epsilon)  # Adjust to avoid 0 and 1

    # Calculate log-likelihood for null hypothesis (all groups have same probability)
    llf_null = np.sum(bernoulli.logpmf(y, prob_null))

    # Calculate log-likelihood for full model (each group has its own probability)
    llf_full = sum([np.sum(bernoulli.logpmf(y[group_labels == group], prob)) 
                    for group, prob in zip(unique_groups, prob_hat)])

    # Compute likelihood ratio test statistic
    lr_statistic = -2 * (llf_null - llf_full)
    
    # Calculate p-value using binomial test
    p_value = binomtest(sum(y), len(y), prob_null, alternative='two-sided').pvalue

    # Calculate effect size (absolute difference between group probabilities)
    effect_size = abs(prob_hat[0] - prob_hat[1])

    # Return results as a pandas Series
    return pd.Series({'lr_statistic': lr_statistic, 'p_value': p_value, 'effect_size': effect_size})

def compute_kruskal_wallis_test(row, groups):
    """
    Compute Kruskal-Wallis H-test for a given row of binary data.

    Args:
        row (pandas.Series): A row of binary data (values: 0 or 1) from a DataFrame.
        groups (dict): A dictionary mapping group names to the corresponding column indices.

    Returns:
        pandas.Series: A series containing the H-statistic, p-value and effect size.
    """
    group_data = {group: row[samples].dropna().values for group, samples in groups.items() if all(sample in row.index for sample in samples)}
  
    # Check if we have at least two groups with data
    if len(group_data) < 2 or not all(len(data) > 0 for data in group_data.values()):
        return pd.Series({'h_statistic': np.nan, 'p_value': np.nan, 'effect_size': np.nan})

    # Prepare data for Kruskal-Wallis test
    samples = [data for data in group_data.values()]
    
    # Check if all numbers are identical across ALL samples
    all_data = np.concatenate(samples)
    if np.all(all_data == all_data[0]):
        return pd.Series({'h_statistic': np.nan, 'p_value': np.nan, 'effect_size': np.nan})

    try:
        # Compute Kruskal-Wallis H-test
        h_statistic, p_value = kruskal(*samples)

        # Compute effect size (eta-squared)
        n = sum(len(s) for s in samples)
        effect_size = (h_statistic - len(samples) + 1) / (n - len(samples))
    except ValueError as e:
        logging.warning(f"ValueError in Kruskal-Wallis test: {e}")
        return pd.Series({'h_statistic': np.nan, 'p_value': np.nan, 'effect_size': np.nan})

    return pd.Series({'h_statistic': h_statistic, 'p_value': p_value, 'effect_size': effect_size})

def compute_likelihood_ratio_test_poisson(row, groups):
    """
    Compute Poisson Likelihood Ratio Test for a given row of count data.

    Args:
        row (pandas.Series): A row of count data from a DataFrame.
        groups (dict): A dictionary mapping group names to the corresponding column indices.

    Returns:
        pandas.Series: A series containing the LRT statistic, p-value and effect size.
                       Returns NaN for all values if the test cannot be performed.
    """
    group_data = {group: row[samples].dropna().values for group, samples in groups.items() if all(sample in row.index for sample in samples)}
  
    if len(group_data) < 2 or not all(len(data) > 0 for data in group_data.values()):
        logging.info("Insufficient data for Poisson LRT.")
        return pd.Series({'lr_statistic': np.nan, 'p_value': np.nan, 'effect_size': np.nan})

    all_data = np.concatenate(list(group_data.values()))
    
    # Calculate pooled counts
    K1, K2 = [np.sum(data) for data in group_data.values()]
    N1, N2 = [len(data) for data in group_data.values()]

    # Maximum likelihood estimates with small epsilon to avoid divide by zero
    epsilon = 1e-10
    theta_hat_1 = max(K1 / N1, epsilon)
    theta_hat_2 = max(K2 / N2, epsilon)
    theta_hat = max((K1 + K2) / (N1 + N2), epsilon)

    # Log-likelihoods
    ll_alt = (K1 * np.log(theta_hat_1) - theta_hat_1 * N1 +
              K2 * np.log(theta_hat_2) - theta_hat_2 * N2)
    ll_null = ((K1 + K2) * np.log(theta_hat) - theta_hat * (N1 + N2))

    # Likelihood ratio statistic
    lr_statistic = 2 * (ll_alt - ll_null)

    # Check if lr_statistic is valid
    if np.isnan(lr_statistic) or lr_statistic < 0:
        if lr_statistic < -1e-10:  # Allow for small negative values due to floating point errors
            logging.info(f"Invalid LR statistic: {lr_statistic}")
            return pd.Series({'lr_statistic': np.nan, 'p_value': np.nan, 'effect_size': np.nan})
        else:
            lr_statistic = 0  # Set very small negative values to zero

    # P-value (using chi-square distribution with 1 degree of freedom)
    p_value = 1 - chi2.cdf(lr_statistic, df=1)
    
    # Effect size (log ratio of rates)
    if theta_hat_1 == 0 or theta_hat_2 == 0:
        effect_size = np.nan
    else:
        effect_size = np.log(theta_hat_1 / theta_hat_2)

    return pd.Series({'lr_statistic': lr_statistic, 'p_value': p_value, 'effect_size': effect_size})

def perform_analysis(df, morph_groups, comparison, test_type, binarization_threshold=0):
    """
    Perform statistical analysis on the given dataframe using the specified test type.

    Args:
    df (dask.dataframe.DataFrame): The input dataframe containing the unitig data.
    morph_groups (dict): A dictionary mapping morph types to lists of sample names.
    comparison (str): The type of comparison to perform, either 'Group1vsGroup2' or 'overall'.
    test_type (str): The type of statistical test to perform ('poisson', 'bernoulli', 'kw', or 'chi2').
    binarization_threshold (float, optional): The threshold for binarizing the data. Used only for
                                              'bernoulli', 'kw', and 'chi2' tests. Default is 0.

    Returns:
    pandas.DataFrame or None: A dataframe containing the results of the statistical tests, including test
                              statistics, p-values, effect sizes, and adjusted p-values. Returns None if an
                              error occurs during the analysis.
    """
    if comparison != 'overall':
        desired_groups = [group.strip() for group in comparison.split('vs')]
        groups = {k: morph_groups[k] for k in morph_groups if k in desired_groups}
        if not groups:
            logging.error(f"Error: Specified groups in comparison ({desired_groups}) not found.")
            return None
    else:
        groups = morph_groups

    relevant_columns = [column for group in groups.values() for column in group]
    if 'Unitigs_id' in relevant_columns: relevant_columns.remove('Unitigs_id')

    df_slice = df[relevant_columns]
    
    try:
        if test_type == 'poisson':
            meta = pd.Series({
                'lr_statistic': np.float64,
                'p_value': np.float64,
                'effect_size': np.float64
            }, index=['lr_statistic', 'p_value', 'effect_size'])
            results = df_slice.apply(lambda row: compute_likelihood_ratio_test_poisson(row, groups), axis=1, meta=meta)
        else:
            # Binarize data for non-Poisson tests
            df_binarized = (df_slice > binarization_threshold).astype(int)
            logging.info("Data was binarized for non-Poisson test.")
            logging.info(f"Binarized data sample:\n{df_binarized.head()}")
            
            if test_type == 'kw':
                meta = pd.Series({
                    'h_statistic': np.float64,
                    'p_value': np.float64,
                    'effect_size': np.float64
                }, index=['h_statistic', 'p_value', 'effect_size'])
                results = df_binarized.apply(lambda row: compute_kruskal_wallis_test(row, groups), axis=1, meta=meta)
            elif test_type == 'chi2':
                meta = pd.Series({
                    'chi2_statistic': np.float64,
                    'p_value': np.float64,
                    'effect_size': np.float64
                }, index=['chi2_statistic', 'p_value', 'effect_size'])
                results = df_binarized.apply(lambda row: compute_chi_square_test(row, groups), axis=1, meta=meta)
            elif test_type == 'bernoulli':
                meta = pd.Series({
                    'lr_statistic': np.float64,
                    'p_value': np.float64,
                    'effect_size': np.float64
                }, index=['lr_statistic', 'p_value', 'effect_size'])
                results = df_binarized.apply(lambda row: compute_likelihood_ratio_test_bernoulli(row, groups), axis=1, meta=meta)
            else:
                raise ValueError("Invalid test type specified.")
    except Exception as e:
        logging.error(f"Error during analysis: {e}")
        return None

    results_df = results.compute()
    valid_rows = results_df.notna().all(axis=1)
    
    if valid_rows.sum() > 0:
        valid_p_values = results_df.loc[valid_rows, 'p_value']
        # Bonferroni correction
        adjusted_p_values = np.minimum(valid_p_values * valid_rows.sum(), 1.0)
        results_df['adjusted_p_value'] = np.nan
        results_df.loc[valid_rows, 'adjusted_p_value'] = adjusted_p_values
    else:
        results_df['adjusted_p_value'] = np.nan

    total_tests = len(results_df)
    valid_tests = valid_rows.sum()
    invalid_tests = total_tests - valid_tests
    
    logging.info(f"Total tests attempted: {total_tests}")
    logging.info(f"Valid tests (used in BH correction): {valid_tests} ({valid_tests/total_tests*100:.2f}%)")
    logging.info(f"Tests where the statistical test couldn't be performed: {invalid_tests} ({invalid_tests/total_tests*100:.2f}%)")
    
    if valid_tests > 0:
        significant_tests = (results_df.loc[valid_rows, 'adjusted_p_value'] < 0.05).sum()
        logging.info(f"Tests with adjusted p-value < 0.05: {significant_tests} ({significant_tests/valid_tests*100:.2f}% of valid tests)")
    else:
        logging.info("No valid tests to assess significance.")

    return results_df

def check_columns(file_path):
    """
    Check the number of columns in a CSV file.

    Args:
    file_path (str): Path to the CSV file.

    Returns:
    tuple: Number of columns and actual column names.
    """
    small_df = dd.read_csv(file_path, sample=1000000)  # only sample 1 MB of data
    actual_columns = small_df.columns
    num_columns = len(actual_columns)
    logging.info(f"Sample data check: The CSV file appears to have {num_columns} columns.")
    return num_columns, actual_columns

def check_file_format(file_path):
    """
    Check the format of the input file.

    Args:
    file_path (str): Path to the input file.

    Returns:
    tuple: (is_binary, num_columns, has_header)
    """
    with open(file_path, 'r') as f:
        first_line = f.readline().strip()
        second_line = f.readline().strip()

    # Check if the first line contains column names
    has_header = not first_line.split()[0].replace('.', '').isdigit()

    # Check if the data is binary (0s and 1s only)
    is_binary = all(val in ['0', '1'] for val in second_line.split()[1:])

    num_columns = len(first_line.split())

    logging.info(f"File format: Binary: {is_binary}, Columns: {num_columns}, Has Header: {has_header}")
    return is_binary, num_columns, has_header

def main(args):
    """
    Main function to set up the environment, process input, and save outputs.

    Args:
    args (Namespace): Argument namespace containing command line arguments.
    """
    logging.info("Starting main function...")

    # Check if output file already exists
    formatted_comparison = args.comparison.strip().replace(' ', '_').replace('/', '').replace('\\', '')
    filename = f"{formatted_comparison}_{args.test_type}_bin_thresh_{args.bin_thresh:.2f}.csv"
    output_file_path = os.path.join(args.output_dir, filename)

    if os.path.exists(output_file_path):
        logging.error(f"Error: Output file {output_file_path} already exists. Exiting to prevent overwriting.")
        sys.exit(1)  # Exit with error code 1

    # Proceed with the rest of the function if the file doesn't exist
    client = setup_dask_client(args.threads, args.memory)
    metadata = load_metadata(METADATA_FILE)
    morph_groups = map_runs_to_morph(metadata)

    # Load column names from the file
    with open(COLUMN_NAMES_FILE, 'r') as f:
        column_names = [line.strip() for line in f]
    column_names = [name.replace("data/unitigs/", "").replace(".unitigs.fa.gz", "") for name in column_names]

    # Ensure 'Unitigs_id' is the first column, but only if it's not already there
    if 'Unitigs_id' not in column_names:
        column_names.insert(0, 'Unitigs_id')
    elif column_names[0] != 'Unitigs_id':
        column_names.remove('Unitigs_id')
        column_names.insert(0, 'Unitigs_id')

    logging.info(f"Loaded {len(column_names)} column names from {COLUMN_NAMES_FILE}")

    # Check the format of the input file
    is_binary, num_columns, has_header = check_file_format(args.unitig_matrix)

    if is_binary:
        # For binary data
        dtype_dict = {col: float for col in column_names}  # Use float to handle values like 0.94
        dtype_dict['Unitigs_id'] = str

        df = dd.read_csv(args.unitig_matrix, 
                        dtype=dtype_dict,
                        sep=',')  # CSV format
        
        logging.info("Binary data detected. Read as CSV.")
    else:
        # For non-binary data
        dtype_dict = {col: float for col in column_names if col != "Unitigs_id"}
        dtype_dict['Unitigs_id'] = str

        df = dd.read_csv(args.unitig_matrix, 
                        names=column_names,
                        dtype=dtype_dict,
                        sep='\s+',  # Whitespace-separated
                        header=None)
        
        logging.info("Non-binary data detected. Read as whitespace-separated.")

    logging.info("Unitig data loaded into Dask DataFrame.")
    logging.info(f"DataFrame columns: {df.columns}")
    logging.info(f"DataFrame sample:\n{df.head()}")

    # Perform the analysis
    results = perform_analysis(df, morph_groups, args.comparison, args.test_type, args.bin_thresh)

    if results is not None:
        os.makedirs(args.output_dir, exist_ok=True)
        logging.info(f"Attempting to save results to: {output_file_path}")

        try:
            results.to_csv(output_file_path, index=False)
            logging.info(f"Results have been saved to {output_file_path}")
        except Exception as e:
            logging.error(f"Failed to save results: {e}")
    else:
        logging.warning("No results were generated.")

    client.close()
    logging.info("Dask client closed. Script completed.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform statistical tests on unitig data.")
    parser.add_argument("-t", "--threads", type=int, required=True, help="Number of threads to use")
    parser.add_argument("-m", "--memory", required=True, help="Amount of memory to use (e.g., '32GB')")
    parser.add_argument("-u", "--unitig_matrix", required=True, help="Path to the unitig matrix file")
    parser.add_argument("-c", "--comparison", required=True,
                    help="Type of comparison to perform, formatted as 'Group1vsGroup2' (e.g., 'AvsB') or 'overall'")
    parser.add_argument("-o", "--output_dir", required=True, help="Output directory to save the results")
    parser.add_argument("--test_type", choices=['poisson', 'kw', 'chi2', 'bernoulli'], default='poisson',
                        help="Type of statistical test to perform. 'poisson' for Poisson LRT, 'kw' for Kruskal-Wallis test, 'chi2' for Chi-square test, 'bernoulli' for Bernoulli LRT")
    parser.add_argument("--bin_thresh", type=float, default=0,
                        help="Threshold for binarization (between 0 and 1). Used for KW, Chi-squared, and Bernoulli tests. Values greater than this will be set to 1, else 0. Default is 0.")
    args = parser.parse_args()
    
    # Validate binarization threshold
    if not 0 <= args.bin_thresh <= 1:
        parser.error("Binarization threshold must be between 0 and 1.")

    # Ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)

    # Setup logging
    log_path = setup_logging(args.output_dir)

    # Log all parameters
    logging.info("Script started with the following parameters:")
    for arg, value in vars(args).items():
        logging.info(f"{arg}: {value}")

    # Log the paths of important files
    logging.info(f"Metadata file: {METADATA_FILE}")
    logging.info(f"Column names file: {COLUMN_NAMES_FILE}")
    logging.info(f"Log file created at: {log_path}")

    # Run the main function
    main(args)