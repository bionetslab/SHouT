from squidpy.gr import spatial_neighbors
from scipy.sparse.csgraph import shortest_path


def get_spatial_graph(adata, coord_type='generic', compute_shortest_path_distances=True):
    """Computes spatial graph using Delaunay triangulation and shortest path distances within this graph.

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data object containing spatial omics data with spatial coordinates stored in `adata`.obsm['spatial'].
    coord_type : str
        Type of coordinate system. Valid options are 'generic' and 'generic'.
    compute_shortest_path_distances : bool
        If `True`, shortest path distances are computed in addition to the spatial adjacency matrix.

    Returns
    -------
    adj_matrix : scipy.sparse.csr_matrix
        Binary adjacency matrix of spatial graph.
    shortest_path_distances : scipy.sparse.csr_matrix or None
        Shortest path distances between all pairs of cells or `None` if ``compute_shortest_path_distances`` is `False`.
    """
    adj_matrix, _ = spatial_neighbors(adata, delaunay=True, coord_type=coord_type, copy = True)
    if compute_shortest_path_distances:
        shortest_path_distances = shortest_path(adj_matrix)
        return adj_matrix, shortest_path_distances
    return adj_matrix, None


def get_extended_neighborhoods(shortest_path_distances, radius):
    """For each cell, computes the sets of cells whose shortest path distances do not extend a user-specified radius.

    Parameters
    ----------
    shortest_path_distances : scipy.sparse.csr_matrix
        Shortest path distances between all pairs of cells as computed by ``get_spatial_graph()``.
    radius : int
        User-specified radius.

    Returns
    -------
    extended_neighborhoods : dict[int, list[int]]
        A dictionary with a key for each cell in ``shortest_path_distances`` (index of corresponding row/column) and
        lists of indices of all cells whose shortest path distance from the key cell does not extend ``radius`` as
        values.
    """
    all_close_cell_pairs = (shortest_path_distances <= radius).nonzero()
    all_close_cell_pairs = list(zip(all_close_cell_pairs[0], all_close_cell_pairs[1]))
    extended_neighborhoods = {cell: [] for cell in range(shortest_path_distances.shape[0])}
    for close_cell_pair in all_close_cell_pairs:
        extended_neighborhoods[close_cell_pair[0]].append(close_cell_pair[1])
    return extended_neighborhoods
