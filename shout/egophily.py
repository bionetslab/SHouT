from .utility import *
import numpy as np
import time


def egophily(adata, cluster_key, radius, coord_type='generic', copy=False, shortest_path_distances=None,
             extended_neighborhoods=None):
    """

    Parameters
    ----------
    adata :
    cluster_key :
    radius :
    coord_type :
    copy :
    shortest_path_distances :
    extended_neighborhoods :

    Returns
    -------

    """
    if shortest_path_distances is None:
        _, shortest_path_distances = get_spatial_graph(adata, coord_type)
    if extended_neighborhoods is None:
        extended_neighborhoods = get_extended_neighborhoods(shortest_path_distances, radius)
    egophilies = np.zeros(len(extended_neighborhoods))
    cell_type_map = adata.obs[cluster_key]
    for cell, extended_neighborhood in extended_neighborhoods.items():
        type_of_cell = cell_type_map.iloc[cell]
        local_cell_type_map = cell_type_map.iloc[extended_neighborhood]
        num_cells = len(local_cell_type_map)
        num_cells_of_same_type = local_cell_type_map.value_counts().loc[type_of_cell]
        egophilies[cell] = num_cells_of_same_type / num_cells
    if copy:
        return egophilies
    adata.obs[f'egophily_{radius}'] = egophilies
