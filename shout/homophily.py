from .utility import retrieve_extended_neighborhoods
import numpy as np


def global_homophily(adata, adj_matrix, adj_matrix_homophilic, copy=False):
    """

    Parameters
    ----------
    adata :
    adj_matrix :
    adj_matrix_homophilic :

    Returns
    -------

    """
    if copy:
        return adj_matrix_homophilic.sum() / adj_matrix.sum()
    adata.uns['global_homophily'] = adj_matrix_homophilic.sum() / adj_matrix.sum()


def local_homophily(adata, shortest_path_distances, radius, adj_matrix, adj_matrix_homophilic, copy=False):
    """

    Parameters
    ----------
    adata :
    shortest_path_distances :
    radius :
    adj_matrix :
    adj_matrix_homophilic :

    Returns
    -------

    """
    extended_neighborhoods = retrieve_extended_neighborhoods(shortest_path_distances, radius)
    local_homophilies = np.zeros(len(extended_neighborhoods))
    for cell, extended_neighborhood in extended_neighborhoods.items():
        sub_adj_matrix = adj_matrix[np.ix_([extended_neighborhood], [extended_neighborhood])]
        sub_adj_matrix_homophilic = adj_matrix_homophilic[np.ix_([extended_neighborhood], [extended_neighborhood])]
        local_homophilies[cell] = sub_adj_matrix_homophilic.sum() / sub_adj_matrix.sum()
    if copy:
        return local_homophilies
    adata.obsm['local_homophily'] = local_homophilies
