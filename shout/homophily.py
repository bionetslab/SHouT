from .utility import *
import numpy as np
import time


def global_homophily(adata, cluster_key, coord_type='generic', copy=False, adj_matrix=None, adj_matrix_homophilic=None):
    """

    Parameters
    ----------
    adata :
    cluster_key :
    coord_type :
    copy :
    adj_matrix :
    adj_matrix_homophilic :

    Returns
    -------

    """
    start=time.time()
    if adj_matrix is None:
        adj_matrix, _ = get_spatial_graph(adata, coord_type, compute_shortest_path_distances=False)
    if adj_matrix_homophilic is None:
        adj_matrix_homophilic = get_homophilic_edges(adata, cluster_key, adj_matrix)
    end=time.time()
    time_elapsed=end-start
    if copy:
        return adj_matrix_homophilic.sum() / adj_matrix.sum(), time_elapsed
    adata.uns['global_homophily'] = adj_matrix_homophilic.sum() / adj_matrix.sum()
    adata.uns['global_homophily_TIME'] = time_elapsed


def local_homophily(adata, cluster_key, radius, coord_type='generic', copy=False, adj_matrix=None,
                    adj_matrix_homophilic=None, shortest_path_distances=None, extended_neighborhoods=None):
    """

    Parameters
    ----------
    adata :
    cluster_key :
    radius :
    coord_type :
    copy :
    adj_matrix :
    adj_matrix_homophilic :
    shortest_path_distances :
    extended_neighborhoods :

    Returns
    -------

    """
    start=time.time()
    if shortest_path_distances is None or adj_matrix is None:
        shortest_path_distances, adj_matrix = get_spatial_graph(adata, coord_type)
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
    end=time.time()
    time_elapsed=end-start
    if copy:
        return local_homophilies, time_elapsed
    adata.obs[f'local_homophily_{radius}'] = local_homophilies
    adata.uns[f'local_homophily_{radius}_TIME'] = time_elapsed


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
