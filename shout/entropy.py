from .utility import *
from scipy.stats import entropy
import numpy as np
from scipy.sparse.csgraph import shortest_path


def global_entropy(adata, cluster_key, normalize=True, num_cell_types=None, copy=False):
    """Global entropy returns Shannon's entropy value of the network, normalized by total number of cell types present in input data.
        Shannon's entropy is a measure of heterogeneity. Value of 1 implies a randomly generated network while value of 0 implies a network with only one celltype/ cluster.

    Parameters
    ----------
    adata : anndata.AnnData (Mandatory parameter)
        Annotated data object containing spatial omics data with spatial coordinates stored in `adata`.obsm['spatial'].
        For more info, go to https://anndata.readthedocs.io/en/latest/.
    
    cluster_key : str (Mandatory parameter)
        adata.obs[cluster_key] contains key where clustering/ cell type annotations are stored.
    
    normalize : bool (Optional parameter | default True)
        If $normalize = True$, normalize by number of cell types in adata.obs[cluster_key] when calculating entropy scores.
    
    num_cell_types : int (Optional parameter | default None)
        If $num_cell_types=None$, num_cell_types is the number of cell types present in adata.obs[cluster_key].
    
    copy : bool (Optional parameter | default False)
        $copy = True$ returns all scores as a dict, $copy = False$ saves all scores as part of the input anndata object "adata".

    Returns
    -------

    """
    cell_type_map = adata.obs[cluster_key]
    if num_cell_types is None:
        num_cell_types = cell_type_map.nunique()
    normalization_constant = np.log2(num_cell_types) if normalize else 1
    num_cells = len(cell_type_map)
    cell_type_frequencies = (cell_type_map.value_counts() / num_cells).to_numpy()
    if copy:
        return entropy(cell_type_frequencies, base=2) / normalization_constant
    adata.uns['global_entropy'] = entropy(cell_type_frequencies, base=2) / normalization_constant


def local_entropy(adata, cluster_key, radius, normalize=True, num_cell_types=None, copy=False,
                  shortest_path_distances=None, extended_neighborhoods=None):
    """Local entropy returns Shannon's entropy value per node within the vicinity of input parameter "radius", normalized by total number of cell types present in input data.
        Shannon's entropy is a measure of heterogeneity. Value of 1 implies random distribution while value of 0 implies a homogeneous subnetwork, around the node.

    Parameters
    ----------
    adata : anndata.AnnData (Mandatory parameter)
        Annotated data object containing spatial omics data with spatial coordinates stored in `adata`.obsm['spatial'].
        For more info, go to https://anndata.readthedocs.io/en/latest/.
    
    cluster_key : str (Mandatory parameter)
        adata.obs[cluster_key] contains key where clustering/ cell type annotations are stored.
    
    radius : int (Mandatory parameter)
        n-hop neighbor over which local homophily scores are to be calculated.
    
    normalize : bool (Optional parameter | default True)
        If $normalize = True$, normalize by number of cell types in adata.obs[cluster_key] when calculating entropy scores.
    
    copy : bool (Optional parameter | default False)
        $copy = True$ returns all scores as a dict, $copy = False$ saves all scores as part of the input anndata object "adata".
    
    shortest_path_distances : scipy.sparse.csr_matrix (Optional parameter | condition <symmetrical matrix with zero diagonal> | default None)
        Shortest path distances between all pairs of cells as computed by ``.utility.get_spatial_graph()``.
        If $shortest_path_distances = None$, shortest_path_distances is obtained upon generation of spatial neighbor graph with Delaunay triangulation using ``squidpy.gr.spatial_neighbors()``.
    
    extended_neighborhoods : dict[int, list[int]] (Optional parameter | default None)
        A dictionary with a key for each cell in ``shortest_path_distances`` (index of corresponding row/column) and lists of indices of all cells whose shortest path distance from the key cell does not extend ``radius`` as values.
        If $extended_neighborhoods = None$, extended_neighborhoods is obtained upon generation of spatial neighbors graph with Delaunay triangulation using ``squidpy.gr.spatial_neighbors()``.

    Returns
    -------

    """
    if shortest_path_distances is None:
        adj_matrix = get_spatial_graph(adata)
        shortest_path_distances = shortest_path(adj_matrix)
    if extended_neighborhoods is None:
        extended_neighborhoods = get_extended_neighborhoods(shortest_path_distances, radius)
    local_entropies = np.zeros(len(extended_neighborhoods))
    cell_type_map = adata.obs[cluster_key]
    if num_cell_types is None:
        num_cell_types = cell_type_map.nunique()
    normalization_constant = np.log2(num_cell_types) if normalize else 1
    for cell, extended_neighborhood in extended_neighborhoods.items():
        local_cell_type_map = cell_type_map.iloc[extended_neighborhood]
        num_cells = len(local_cell_type_map)
        cell_type_frequencies = (local_cell_type_map.value_counts() / num_cells).to_numpy()
        local_entropies[cell] = entropy(cell_type_frequencies, base=2) / normalization_constant
    if copy:
        return local_entropies
    adata.obs[f'local_entropy_{radius}'] = local_entropies

