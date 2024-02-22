from .utility import *
import numpy as np
from scipy.sparse.csgraph import shortest_path


def egophily(adata, cluster_key, radius, copy=False, shortest_path_distances=None,
             extended_neighborhoods=None):
    """

    Parameters
    ----------
    adata : anndata.AnnData (Mandatory parameter)
        Annotated data object containing spatial omics data with spatial coordinates stored in `adata`.obsm['spatial'].
        For more info, go to https://anndata.readthedocs.io/en/latest/.
    
    cluster_key : str (Mandatory parameter)
        adata.obs[cluster_key] contains key where clustering/ cell type annotations are stored.
    
    radius : int (Mandatory parameter)
        n-hop neighbor over which local homophily scores are to be calculated.
    
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
