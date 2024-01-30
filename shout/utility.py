from squidpy.gr import spatial_neighbors
from scipy.sparse.csgraph import shortest_path
from scipy.sparse import csr_matrix


def compute_spatial_graph(adata, cluster_key, coord_type='generic'):
    """Computes spatial graph using Delaunay triangulation and annotates edges as homophilic or heterophilic.

    Parameters
    ----------
    adata :
    cluster_key :
    coord_type :

    Returns
    -------

    """
    adj_matrix, _ = spatial_neighbors(adata, delaunay=True, coord_type=coord_type, copy = True)
    shortest_path_distances = shortest_path(adj_matrix)
    # adj_matrix_homophilic = csr_matrix(adj_matrix)
    adj_matrix_homophilic=adj_matrix.copy()
    support_adj_matrix = adj_matrix.nonzero()
    for edge in zip(support_adj_matrix[0], support_adj_matrix[1]):
        if adata.obs[cluster_key].iloc[edge[0]] != adata.obs[cluster_key].iloc[edge[1]]:
            adj_matrix_homophilic[edge] = 0
    return adj_matrix, adj_matrix_homophilic, shortest_path_distances


def retrieve_extended_neighborhoods(shortest_path_distances, radius):
    """

    Parameters
    ----------
    adata :
    cell :
    radius :

    Returns
    -------

    """
    all_close_cell_pairs = (shortest_path_distances <= radius).nonzero()
    all_close_cell_pairs = list(zip(all_close_cell_pairs[0], all_close_cell_pairs[1]))
    extended_neighborhoods = {cell: [] for cell in range(shortest_path_distances.shape[0])}
    for close_cell_pair in all_close_cell_pairs:
        extended_neighborhoods[close_cell_pair[0]].append(close_cell_pair[1])
    return extended_neighborhoods
