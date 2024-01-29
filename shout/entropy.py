from .utility import retrieve_extended_neighborhoods
from scipy.stats import entropy
import numpy as np


def local_entropy(adata, cluster_key, shortest_path_distances, radius, copy=False):
    """Computes local entropy scores.

    Parameters
    ----------
    adata : Annotated data object with spatial graph already computed.
    cluster_key : Key of column specifying cell types, states, or clusters.
    radius : Radius to be considered.

    Returns
    -------

    """
    extended_neighborhoods = retrieve_extended_neighborhoods(shortest_path_distances, radius)
    local_entropies = np.zeros(len(extended_neighborhoods))
    for cell, extended_neighborhood in extended_neighborhoods.items():
        cell_type_map = adata[cluster_key].iloc[extended_neighborhood]
        num_cells = len(cell_type_map)
        cell_type_frequencies = (cell_type_map.value_counts() / num_cells).to_numpy()
        local_entropies[cell] = entropy(cell_type_frequencies, base=2)
    if copy:
        return local_entropies
    adata.obsm['local_entropy'] = local_entropies


def global_entropy(adata, cluster_key, copy=False):
    """

    Parameters
    ----------
    adata :
    cluster_key :

    Returns
    -------

    """
    cell_type_map = adata[cluster_key]
    num_cells = len(cell_type_map)
    cell_type_frequencies = (cell_type_map.value_counts() / num_cells).to_numpy()
    if copy:
        return entropy(cell_type_frequencies, base=2)
    adata.uns['global_entropy'] = entropy(cell_type_frequencies, base=2)
