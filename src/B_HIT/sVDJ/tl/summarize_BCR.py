import pandas as pd

def compute_clone_counts(rep_loc, groupby_cols, count_col, extra_cols, count_name='count', freq_name='freq', if_count=True, if_freq=True):
    """
    Compute clone counts and frequencies for specified groups.

    Parameters
    ----------
    rep_loc : pd.DataFrame
        Input DataFrame containing clone information.
    groupby_cols : list
        Columns to group by (e.g., sample, region, family ID).
    count_col : str
        Column to count occurrences of (e.g., 'family_id').
    extra_cols : list
        Additional columns to include before counting.
    count_name : str, optional
        Name for the count result column (default is 'count').
    freq_name : str, optional
        Name for the frequency result column (default is 'freq').
    if_count : bool, optional
        Whether to compute and return counts (default is True).
    if_freq : bool, optional
        Whether to compute and return frequencies (default is True).

    Returns
    -------
    pd.DataFrame
        A DataFrame containing counts, frequencies, or both.
    """
    
    # Step 1: Calculate the value counts (UmiCount)
    rep_agg = rep_loc[groupby_cols + extra_cols + [count_col]].value_counts().reset_index(name='UmiCount')
    
    # Step 2: Compute count per group
    if if_count:
        count_df = rep_agg.groupby(groupby_cols)[count_col].value_counts().reset_index(name=count_name)
        if not if_freq:  # Corrected this check
            return count_df
    
    # Step 3: Compute frequency per group
    if if_freq:
        freq_df = rep_agg.groupby(groupby_cols)[count_col].value_counts(normalize=True).reset_index(name=freq_name)
        if not if_count:  # Corrected this check
            return freq_df
    
    # Step 4: Merge count and frequency data
    merged_df = pd.merge(freq_df, count_df, on=groupby_cols + [count_col])
    
    return merged_df


def compute_richness(df, extra_cols, groupby_cols, richness_name, default_value=0, return_df=False):
    """
    Compute the richness of groups based on unique combinations of key columns.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame containing the data to compute richness.
    extra_cols : list
        List of columns which is not included in groupby_cols but needs to be selected (e.g. ['sample', 'family_id']).
    groupby_cols : list
        Columns that define groups.
    richness_name : str
        Name of the richness column.
    default_value : int, optional (default=0)
        Default value to return if a group is missing from the lookup.
    return_df : bool, optional (default=False)
        If True, return the richness DataFrame instead of mapped values.

    Returns
    -------
    pd.DataFrame or pd.Series
        If return_df=True, returns the richness DataFrame.
        If return_df=False, returns a Series with richness values mapped to df.
    """

    # Step 1: Remove duplicates to get unique combinations
    richness_df = df[extra_cols + groupby_cols].drop_duplicates()

    # Step 2: Compute richness (size of each group)
    richness_df = richness_df.groupby(groupby_cols).size().reset_index(name=richness_name)

    # Step 3: Remove groups with zero richness
    richness_df = richness_df[richness_df[richness_name] != 0].copy()

    # Step 4: Create a dictionary for fast lookup
    richness_dict = richness_df.set_index(groupby_cols)[richness_name].to_dict()

    if return_df:
        return richness_df  # Return full DataFrame
    else:
        return richness_dict


# cloneRich = rep_agg[['sample', 'Cregion_simple', 'family_id', 'Bcell_aggregate_label','Bagg_Anno_res', 'BaggArea' ]].drop_duplicates() 
# cloneRich = cloneRich.groupby(['Bcell_aggregate_label', 'Bagg_Anno_res', 'BaggArea', 'Cregion_simple']).size().reset_index(name='cloneFamilyRichness')
# cloneRich = cloneRich[cloneRich['cloneFamilyRichness']!=0].copy()

# def compute_clone_location_class(rep_loc, family_col, sample_col, region_col, x_col, y_col, agg_label_col):
#     """
#     Compute the classification of clone locations based on B cell aggregation.

#     Parameters:
#     rep_loc : pd.DataFrame
#         Input DataFrame containing clone information.
#     family_col : str
#         Column name for the clone family ID.
#     sample_col : str
#         Column name for the sample identifier.
#     region_col : str
#         Column name for the region classification.
#     x_col : str
#         Column name for the x-coordinate.
#     y_col : str
#         Column name for the y-coordinate.
#     agg_label_col : str
#         Column name for the B cell aggregation label.

#     Returns:
#     pd.DataFrame
#         Processed DataFrame with clone location classifications.
#     """

#     import numpy as np
#     import pandas as pd

#     # Step 1: Create a copy of the input DataFrame
#     rep_agg_all = rep_loc.copy()

#     # Step 2: Create a new column based on B cell aggregation status
#     rep_agg_all['Bcell_aggregate_label_2'] = np.array(rep_agg_all[agg_label_col])
#     rep_agg_all.loc[~rep_agg_all['Bcell_aggregate_label_2'].isna(), 'Bcell_aggregate_label_2'] = 'Aggregates'
#     rep_agg_all.loc[rep_agg_all['Bcell_aggregate_label_2'].isna(), 'Bcell_aggregate_label_2'] = 'Scattered'

#     # Step 3: Compute UMI counts
#     rep_agg_all_1 = rep_agg_all[[family_col, sample_col, region_col, x_col, y_col]].value_counts().reset_index(name='UmiCount')

#     # Step 4: Compute family_id counts per sample and region
#     _Index_compute_count = compute_clone_counts(rep_agg_all_1, [sample_col, region_col], family_col, [], count_name='bin20Count', freq_name=None, if_count=True, if_freq=False)
#     family_id_dict = _Index_compute_count.set_index([sample_col, region_col, family_col]).to_dict()['bin20Count']

#     # Step 5: Remove duplicates and map family_id_bin20Count
#     rep_agg_all_1 = rep_agg_all[[sample_col, region_col, family_col, 'Bcell_aggregate_label_2']].drop_duplicates()
#     rep_agg_all_1['family_id_bin20Count'] = rep_agg_all_1.apply(lambda row: family_id_dict.get((row[sample_col], row[region_col], row[family_col])), axis=1)

#     # Step 6: Pivot table to classify clones based on aggregation
#     df = rep_agg_all_1.pivot_table(index=[sample_col, region_col, family_col], columns=['Bcell_aggregate_label_2'])
#     df.columns = df.columns.get_level_values(1)
#     df['class'] = 'Not_InAgg'
#     df.loc[df['Scattered'].isna(), 'class'] = 'InAgg'
#     df.loc[(~df['Scattered'].isna()) & (~df['Aggregates'].isna()), 'class'] = 'Shared'
#     df = df.reset_index()

#     # Step 7: Map classification back to the original DataFrame
#     cloneLocClassDict = df[[sample_col, region_col, family_col, 'class']].set_index([sample_col, region_col, family_col]).to_dict()['class']
#     rep_agg = rep_loc[[family_col, sample_col, region_col, x_col, y_col]].value_counts().reset_index(name='UmiCount')
#     rep_agg['familyLocClass'] = rep_agg.apply(lambda row: cloneLocClassDict.get((row[sample_col], row[region_col], row[family_col])), axis=1)

#     # Step 8: Compute clone counts and frequencies by location class
#     _Index_compute = compute_clone_counts(rep_agg, [sample_col, region_col, 'familyLocClass'], family_col, [], count_name='count', freq_name='freq', if_count=True, if_freq=True)

#     return _Index_compute
