from .entropy import *
from .homophily import *
from .egophily import *
from .utility import *


def all_scores(adata, cluster_key, radii, normalize=True, num_cell_types=None, coord_type='generic', copy=False):
    """

    Parameters
    ----------
    adata :
    cluster_key :
    radii :
    normalize :
    num_cell_types : 
    coord_type :
    copy :

    Returns
    -------

    """
    adj_matrix, shortest_path_distances = get_spatial_graph(adata, coord_type)
    adj_matrix_homophilic = get_homophilic_edges(adata, cluster_key, adj_matrix)
    scores = dict()
    scores['global_entropy'] = global_entropy(adata, cluster_key, normalize=normalize, num_cell_types=num_cell_types, copy=copy)
    scores['global_homophily'] = global_homophily(adata, cluster_key, coord_type=coord_type, copy=copy,
                                                  adj_matrix=adj_matrix, adj_matrix_homophilic=adj_matrix_homophilic)
    for radius in radii:
        extended_neighborhoods = get_extended_neighborhoods(shortest_path_distances, radius)
        scores[f'local_entropy_{radius}', f'local_entropy_{radius}_TIME'] = local_entropy(adata, cluster_key, radius, normalize=normalize, num_cell_types=num_cell_types, 
                                                          copy=copy, shortest_path_distances=shortest_path_distances,
                                                          extended_neighborhoods=extended_neighborhoods)
        scores[f'local_homophily_{radius}', f'local_homophily_{radius}_TIME'] = local_homophily(adata, cluster_key, radius, coord_type=coord_type, copy=copy,
                                                              adj_matrix=adj_matrix,
                                                              adj_matrix_homophilic=adj_matrix_homophilic,
                                                              extended_neighborhoods=extended_neighborhoods)
        scores[f'egophily_{radius}', f'egophily_{radius}_TIME'] = egophily(adata, cluster_key, radius, coord_type=coord_type, copy=copy,
                                                shortest_path_distances=shortest_path_distances,
                                                extended_neighborhoods=extended_neighborhoods)
    if copy:
        return scores
