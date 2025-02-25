import pandas as pd
from scipy.stats import pearsonr

def compute_correlation(cloneRich, groupby_cols, corr1, corr2, save=False, path=None, compute_corr_matrix=False):
    """
    Compute Pearson correlation between two variables, grouped by specific columns.
    
    Parameters
    ----------
    cloneRich : pd.DataFrame
        The input DataFrame containing the data.
    groupby_cols : list of str
        List of columns to group by (e.g., ['Cregion_simple', 'tissue']).
    corr1 : str
        The name of the first column to compute correlation for (e.g., 'BaggArea').
    corr2 : str
        The name of the second column to compute correlation for (e.g., 'gini_index').
    save : bool, optional
        Whether to save the output to CSV. Default is False.
    path : str, optional
        The file path to save the output CSV. Default is None.
    compute_corr_matrix : bool, optional
        Whether to compute and return the correlation matrix and p-value matrix. Default is False.
        
    Returns
    -------
    pd.DataFrame
        A DataFrame containing the correlation and p-value for each group.
    pd.DataFrame, optional
        A correlation matrix (if compute_corr_matrix is True).
    pd.DataFrame, optional
        A p-value matrix (if compute_corr_matrix is True).
    """
    
    # Initialize an empty dictionary to store the correlation results for each group
    grouped_correlations = {}
    
    # Loop through each group defined by 'groupby_cols' and compute Pearson correlation for 'corr1' and 'corr2'
    for group, group_data in cloneRich.groupby(groupby_cols):
        correlation, p_value = pearsonr(group_data[corr1], group_data[corr2])  # Pearson correlation
        grouped_correlations[group] = (correlation, p_value)  # Store the results for each group
    
    # Convert the dictionary of results into a DataFrame
    corrDf = pd.DataFrame(grouped_correlations).T  # Transpose to get groups as rows
    corrDf.columns = ['corr', 'Pvalue']  # Name the columns
    

    # If 'save' is True, write the DataFrame to a CSV file at the specified 'path'
    if save and path:
        corrDf.to_csv(path)
    
    # If compute_corr_matrix is True, compute the correlation matrix and p-value matrix
    if compute_corr_matrix:
        # Reset the index and set the new column names for the group keys
        corrDf = corrDf.reset_index(names=['chain', 'region'])  # Set new column names for the groupby columns
    
        corrmat = corrDf[['chain', 'region', 'corr']].copy().pivot_table(index='chain', columns='region')
        corrmat.columns = corrmat.columns.get_level_values(1)  # Flatten the multi-index of columns
        
        Pmat = corrDf[['chain', 'region', 'Pvalue']].copy().pivot_table(index='chain', columns='region')
        Pmat.columns = Pmat.columns.get_level_values(1)  # Flatten the multi-index of columns
        
        # Save the matrices if 'save' is True
        if save and path:
            corrmat.to_csv(f'{path}_correlation_matrix.csv')
            Pmat.to_csv(f'{path}_pvalue_matrix.csv')

        # Return the matrices
        return corrDf, corrmat, Pmat
    
    # Return the correlation DataFrame if no matrix is needed
    return corrDf