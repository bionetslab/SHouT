from .entropy import *
from .homophily import *
from .egophily import *


def all_scores(adata, cluster_key, radius, normalize=True, coord_type='generic', copy=False):
    """

    Parameters
    ----------
    adata :
    cluster_key :
    radius :
    normalize :
    coord_type :
    copy :

    Returns
    -------

    """
    adj_matrix, shortest_path_distances = get_spatial_graph(adata, coord_type)
    extended_neighborhoods = get_extended_neighborhoods(shortest_path_distances, radius)
    adj_matrix_homophilic = get_homophilic_edges(adata, cluster_key, adj_matrix)
    scores = dict()
    scores['global_entropy'] = global_entropy(adata, cluster_key, normalize=normalize, copy=copy)
    scores[f'local_entropy_radius'] = local_entropy(adata, cluster_key, radius, normalize=normalize, copy=copy,
                                                    shortest_path_distances=shortest_path_distances,
                                                    extended_neighborhoods=extended_neighborhoods)
    scores['global_homophily'] = global_homophily(adata, cluster_key, coord_type=coord_type, copy=copy,
                                                  adj_matrix=adj_matrix, adj_matrix_homophilic=adj_matrix_homophilic)
    scores[f'local_homophily_{radius}'] = local_homophily(adata, cluster_key, radius, coord_type=coord_type, copy=copy,
                                                          adj_matrix=adj_matrix,
                                                          adj_matrix_homophilic=adj_matrix_homophilic,
                                                          extended_neighborhoods=extended_neighborhoods)
    scores[f'egophily_{radius}'] = egophily(adata, cluster_key, radius, coord_type=coord_type, copy=copy,
                                            shortest_path_distances=shortest_path_distances,
                                            extended_neighborhoods=extended_neighborhoods)
    if copy:
        return scores
