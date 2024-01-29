from .utility import retrieve_extended_neighborhoods
import numpy as np


def egophily(adata, cluster_key, shortest_path_distances, radius, copy=False):
    """

    Parameters
    ----------
    adata :
    cluster_key :
    radius :

    Returns
    -------

    """
    extended_neighborhoods = retrieve_extended_neighborhoods(shortest_path_distances, radius)
    egophilies = np.zeros(len(extended_neighborhoods))
    for cell, extended_neighborhood in extended_neighborhoods.items():
        type_of_cell = adata[cluster_key].iloc[cell]
        cell_type_map = adata[cluster_key].iloc[extended_neighborhood]
        num_cells = len(cell_type_map)
        num_cells_of_same_type = cell_type_map.value_counts().loc[type_of_cell]
        egophilies[cell] = num_cells_of_same_type / num_cells
    if copy:
        return egophilies
    adata.obsm['egophily'] = egophilies
