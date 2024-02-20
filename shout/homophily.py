from .utility import *
import numpy as np
from scipy.sparse.csgraph import shortest_path


def global_homophily(adata, cluster_key, copy=False, adj_matrix=None, adj_matrix_homophilic=None):
    """

    Parameters
    ----------
    adata :
    cluster_key :
    copy :
    adj_matrix :
    adj_matrix_homophilic :

    Returns
    -------

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
    """

    Parameters
    ----------
    adata (Mandatory parameter | type <AnnData>) - Annotated data object containing spatial information. For more info, go to https://anndata.readthedocs.io/en/latest/.
    cluster_key (Mandatory parameter | type <str>) - adata.obs[cluster_key] contains key where clustering/ cell type annotations are stored.
    radius (Mandatory parameter | type <int> or <float>) - n-hop neighbor over which local homophily scores are to be calculated.
    copy (Optional parameter | type <bool> | default False) - $copy = True$ returns all scores as a dict, $copy = False$ saves all scores as part of the input anndata object "adata".
    adj_matrix (Optional parameter | type <symmetrical matrix of 0s and 1s, with zero diagonal> | default None) - symmetrical matrix where 1 represents an edge between cells, and 0 represents no edge between points. If $adj_matrix = None$, adjacency matrix is obtained upon generation of spatial neighbor graph with Delaunay triangulation.
    adj_matrix_homophilic (Optional parameter | type <symmetrical matrix of 0s and 1s, with zero diagonal> | default None) - Adjacency matrix retaining only homophilic edges, that is edges where both nodes are of the same cell type (or cluster). If None, adj_matrix_homophilic is obtained upon generation of spatial neighbor graph with Delaunay triangulation.
    shortest_path_distances (Optional parameter | type <symmetrical matrix with zero diagonal> | default None) - If None, shortest_path_distances is obtained upon generation of spatial neighbor graph with Delaunay triangulation.
    extended_neighborhoods (Optional parameter | type <dict[int, list[int]]> | default None) - A dictionary with a key for each cell, and value being the list of neighbors within less than "radius" distance of the cell. If None, extended_neighborhoods is obtained upon generation of spatial neighbors graph with Delaunay triangulation.

    Returns
    -------

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
    """

    Parameters
    ----------
    adata :
    cluster_key :
    adj_matrix :

    Returns
    -------

    """
    adj_matrix_homophilic = adj_matrix.copy()
    support_adj_matrix = adj_matrix.nonzero()
    for edge in zip(support_adj_matrix[0], support_adj_matrix[1]):
        if adata.obs[cluster_key].iloc[edge[0]] != adata.obs[cluster_key].iloc[edge[1]]:
            adj_matrix_homophilic[edge] = 0
    return adj_matrix_homophilic
