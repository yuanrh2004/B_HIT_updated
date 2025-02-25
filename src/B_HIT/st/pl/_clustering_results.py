import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import anndata as ad

def clustering_results(adata: ad.AnnData, cluster_column: str, color_dict: dict = None, spot_size: int = 50):
    """
    Plots the clustering results, displaying the cluster labels for each data point using the specified colors.
    
    Parameters
    ----------
    adata: AnnData object
        containing spatial information and clustering labels.
    - cluster_column: str
        the name of the column containing the cluster labels.(eg. the name of the column in use_rep in cluster_Auto_k)
    - color_dict: dict
        a dictionary specifying the color for each cluster label (optional).
    - spot_size: int
        the size of the data points (default is 50).
    """
    color_lst = [ 
    '#C9E9D9', # 0
    '#ffc0cb', # 1
    '#FEFDD3', # 2 #F8E2E1  ok 
    '#C1B1D2', # 3 ok 
    '#EEE5F8', # 4 ok
    '#4074AA', # 5 ok 
    '#EF9194', # 6 ok 
    '#D0353A', # 7 ok
    '#78A5C9', # 8 ok
    '#798af2', # 9 ok
    '#E4F3FC', # 10  cceeff ok
    '#8BD0D8', # 11 ok
    '#F8E886', # 12 f5cfdc
    '#A4C4D9'
]

    color_dict = {}
    cluster_lst = np.sort(list(set(adata.obs[cluster_column])))
    for i in range(len(cluster_lst)): 
        color_dict[cluster_lst[i]] = color_lst[i]
    # Create the plot
    fig, ax = plt.subplots(figsize=(8, 8))

    # Use Scanpy's spatial plotting function
    sc.pl.spatial(adata, color=[cluster_column], ax=ax, show=False, spot_size=spot_size, palette=color_lst)
    
    # Display the plot
    plt.title(f"Clustering results for {cluster_column}")  #这里的title可以改一下
    plt.show()
