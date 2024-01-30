from .utility import *
from scipy.stats import entropy
import numpy as np


def global_entropy(adata, cluster_key, normalize=True, copy=False):
    """

    Parameters
    ----------
    adata :
    cluster_key :
    normalize :
    copy :

    Returns
    -------

    """
    cell_type_map = adata.obs[cluster_key]
    normalization_constant = np.log2(cell_type_map.nunique()) if normalize else 1
    num_cells = len(cell_type_map)
    cell_type_frequencies = (cell_type_map.value_counts() / num_cells).to_numpy()
    if copy:
        return entropy(cell_type_frequencies, base=2) / normalization_constant
    adata.uns['global_entropy'] = entropy(cell_type_frequencies, base=2) / normalization_constant


def local_entropy(adata, cluster_key, radius, normalize=True, coord_type='generic', copy=False,
                  shortest_path_distances=None, extended_neighborhoods=None):
    """

    Parameters
    ----------
    adata :
    cluster_key :
    radius :
    normalize :
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
    local_entropies = np.zeros(len(extended_neighborhoods))
    cell_type_map = adata.obs[cluster_key]
    normalization_constant = np.log2(cell_type_map.nunique()) if normalize else 1
    for cell, extended_neighborhood in extended_neighborhoods.items():
        local_cell_type_map = cell_type_map.iloc[extended_neighborhood]
        num_cells = len(local_cell_type_map)
        cell_type_frequencies = (local_cell_type_map.value_counts() / num_cells).to_numpy()
        local_entropies[cell] = entropy(cell_type_frequencies, base=2) / normalization_constant
    if copy:
        return local_entropies
    adata.obs[f'local_entropy_{radius}'] = local_entropies

