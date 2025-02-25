from .spatial_bcr_desc import compute_index

def compute_grouped_index(_Index_compute, index, groups, column_name, exclude_values=None, exclude_name=None, check_column=None):
    """
    Compute the specified index for each group and handle missing values.
    
    Parameters
    ----------
    _Index_compute : pd.DataFrame
        The input DataFrame containing the data.
    index : str
        The name of the index to compute (e.g., 'gini_index').
    groups : list of str
        List of columns to group by (e.g., ['sample', 'Cregion_simple']).
    column_name : str
        The column to apply the index function to (e.g., 'count').
    exclude_values : list, optional
        List of values to exclude in the result. Default is None.
    exclude_name : str, optional
        The name of the column to check for exclusion values. Default is None.
    check_column : str, optional
        The name of the column to check for NaN values. If not provided, defaults to the name of the index column.
        
    Returns
    -------
    pd.DataFrame
        A DataFrame with the computed index for each group.
    """
    
    
    if check_column is None:
        check_column = index
    # Drop rows with NaN values in the specified column (defaults to `index`)
    # Compute the index for each group
    tmp_df = _Index_compute.groupby(groups)[column_name].apply(lambda x: compute_index(index, x)).reset_index(name=check_column)
    tmp_df = tmp_df.dropna(subset=[check_column])
    
    # Exclude values if specified
    if exclude_values is not None:
        tmp_df = tmp_df[~tmp_df[exclude_name].isin(exclude_values)]
    
    return tmp_df

