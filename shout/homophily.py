from .utility import *
import numpy as np
from scipy.sparse.csgraph import shortest_path


def global_homophily(adata, cluster_key, copy=False, adj_matrix=None, adj_matrix_homophilic=None):
    """ Global homophily is a score that measures the heterogeneity of a network by computing the ratio of the number of edges between nodes of the same type, and the total number of edges in the network.

    Parameters
    ----------
    adata : anndata.AnnData (Mandatory parameter)
        Annotated data object containing spatial omics data with spatial coordinates stored in `adata`.obsm['spatial'].
        For more info, go to https://anndata.readthedocs.io/en/latest/.
    
    cluster_key : str (Mandatory parameter)
        adata.obs[cluster_key] contains key where clustering/ cell type annotations are stored.
    
    copy : bool (Optional parameter | default False)
        $copy = True$ returns all scores as a dict, $copy = False$ saves all scores as part of the input anndata object "adata".
    
    adj_matrix : scipy.sparse.csr_matrix (Optional parameter | condition <symmetrical matrix of 0s and 1s, with zero diagonal> | default None)
        Binary adjacency matrix of spatial graph.
        symmetrical matrix where 1 represents an edge between cells, and 0 represents no edge between points.
        If $adj_matrix = None$, adjacency matrix is obtained upon generation of spatial neighbor graph with Delaunay triangulation using ``squidpy.gr.spatial_neighbors()``.
    
    adj_matrix_homophilic : scipy.sparse.csr_matrix (Optional parameter | condition <symmetrical matrix of 0s and 1s, with zero diagonal> | default None)
        Adjacency matrix retaining only homophilic edges, that is edges where both nodes are of the same cell type (or cluster).
        If $adj_matrix_homophilic = None$, adj_matrix_homophilic is obtained upon generation of spatial neighbor graph with Delaunay triangulation using ``squidpy.gr.spatial_neighbors()``.

    Returns
    -------
    global_homophily : list[float]  
        Returns value if $copy = True$, saves to adata.obs[f'local_homophily_{radius}'] if $copy = False$.

    """
    if adj_matrix is None:
        adj_matrix, _ = get_spatial_graph(adata, compute_shortest_path_distances=False)
    if adj_matrix_homophilic is None:
        adj_matrix_homophilic = get_homophilic_edges(adata, cluster_key, adj_matrix)
    if copy:
        return adj_matrix_homophilic.sum() / adj_matrix.sum()
    adata.uns['global_homophily'] = adj_matrix_homophilic.sum() / adj_matrix.sum()


def local_homophily(adata, cluster_key, radius, copy=False, adj_matrix=None,
                    adj_matrix_homophilic=None, shortest_path_distances=None, extended_neighborhoods=None):
    """ Local homophily is the ratio of the number of nodes of the same class/ cell type forming an edge, and the number nodes of different classes/ cell types forming an edge, within the vicinity of input parameter "radius".
        Local homophily quantifies heterogeneity by rewarding homogeneous edges and penalizing heterogeneous edges.
        In other words, when dissimilar nodes cluster together, we get a higher value of homophily.

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
    
    adj_matrix : scipy.sparse.csr_matrix (Optional parameter | condition <symmetrical matrix of 0s and 1s, with zero diagonal> | default None)
        Binary adjacency matrix of spatial graph.
        symmetrical matrix where 1 represents an edge between cells, and 0 represents no edge between points.
        If $adj_matrix = None$, adjacency matrix is obtained upon generation of spatial neighbor graph with Delaunay triangulation using ``squidpy.gr.spatial_neighbors()``.
    
    adj_matrix_homophilic : scipy.sparse.csr_matrix (Optional parameter | condition <symmetrical matrix of 0s and 1s, with zero diagonal> | default None)
        Adjacency matrix retaining only homophilic edges, that is edges where both nodes are of the same cell type (or cluster).
        If $adj_matrix_homophilic = None$, adj_matrix_homophilic is obtained upon generation of spatial neighbor graph with Delaunay triangulation using ``squidpy.gr.spatial_neighbors()``.
    
    shortest_path_distances : scipy.sparse.csr_matrix (Optional parameter | condition <symmetrical matrix with zero diagonal> | default None)
        Shortest path distances between all pairs of cells as computed by ``.utility.get_spatial_graph()``.
        If $shortest_path_distances = None$, shortest_path_distances is obtained upon generation of spatial neighbor graph with Delaunay triangulation using ``squidpy.gr.spatial_neighbors()``.
    
    extended_neighborhoods : dict[int, list[int]] (Optional parameter | default None)
        A dictionary with a key for each cell in ``shortest_path_distances`` (index of corresponding row/column) and lists of indices of all cells whose shortest path distance from the key cell does not extend ``radius`` as values.
        If $extended_neighborhoods = None$, extended_neighborhoods is obtained upon generation of spatial neighbors graph with Delaunay triangulation using ``squidpy.gr.spatial_neighbors()``.
    
    Returns
    -------
    f'local_homophily_{radius}' : list[float]  
        Returns value if $copy = True$, saves to adata.obs[f'local_homophily_{radius}'] if $copy = False$.

    """
    if adj_matrix is None:
        adj_matrix = get_spatial_graph(adata)
    if shortest_path_distances is None:
        shortest_path_distances = shortest_path(adj_matrix)
    if adj_matrix_homophilic is None:
        adj_matrix_homophilic = get_homophilic_edges(adata, cluster_key, adj_matrix)
    if extended_neighborhoods is None:
        extended_neighborhoods = get_extended_neighborhoods(shortest_path_distances, radius)
    local_homophilies = np.zeros(len(extended_neighborhoods))
    for cell, extended_neighborhood in extended_neighborhoods.items():
        sub_adj_matrix = adj_matrix[np.ix_(extended_neighborhood, extended_neighborhood)]
        sub_adj_matrix_homophilic = adj_matrix_homophilic[np.ix_(extended_neighborhood, extended_neighborhood)]
        sum_of_degrees = sub_adj_matrix.sum()
        if sum_of_degrees == 0:
            local_homophilies[cell] = 0
        else:
            local_homophilies[cell] = sub_adj_matrix_homophilic.sum() / sum_of_degrees
    if copy:
        return local_homophilies
    adata.obs[f'local_homophily_{radius}'] = local_homophilies


def get_homophilic_edges(adata, cluster_key, adj_matrix):
    """ This function returns a homophilic adjacency matrix, i.e., an adjacency matrix that retains edges where both of the nodes belong to the same set.

    Parameters
    ----------
    adata : anndata.AnnData (Mandatory parameter)
        Annotated data object containing spatial omics data with spatial coordinates stored in `adata`.obsm['spatial'].
        For more info, go to https://anndata.readthedocs.io/en/latest/.
    
    cluster_key : str (Mandatory parameter)
        adata.obs[cluster_key] contains key where clustering/ cell type annotations are stored.
    
    adj_matrix : scipy.sparse.csr_matrix (Optional parameter | condition <symmetrical matrix of 0s and 1s, with zero diagonal> | default None)
        Binary adjacency matrix of spatial graph.
        symmetrical matrix where 1 represents an edge between cells, and 0 represents no edge between points.
        If $adj_matrix = None$, adjacency matrix is obtained upon generation of spatial neighbor graph with Delaunay triangulation using ``squidpy.gr.spatial_neighbors()``.
    
    Returns
    -------
    adj_matrix_homophilic : scipy.sparse.csr_matrix (Optional parameter | condition <symmetrical matrix of 0s and 1s, with zero diagonal> | default None)
        Adjacency matrix retaining only homophilic edges, that is edges where both nodes are of the same cell type (or cluster).
        If $adj_matrix_homophilic = None$, adj_matrix_homophilic is obtained upon generation of spatial neighbor graph with Delaunay triangulation using ``squidpy.gr.spatial_neighbors()``.
    
    """
    adj_matrix_homophilic = adj_matrix.copy()
    support_adj_matrix = adj_matrix.nonzero()
    for edge in zip(support_adj_matrix[0], support_adj_matrix[1]):
        if adata.obs[cluster_key].iloc[edge[0]] != adata.obs[cluster_key].iloc[edge[1]]:
            adj_matrix_homophilic[edge] = 0
    return adj_matrix_homophilic
