from squidpy.gr import spatial_neighbors
import networkx as nx
from scipy import sparse
import pandas as pd
from sklearn.metrics.pairwise import pairwise_distances
import numpy as np


def get_spatial_graph(adata, graph_type):
    """Computes spatial graph using Delaunay triangulation and shortest path distances within this graph.

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data object containing spatial omics data with spatial coordinates stored in `adata`.obsm['spatial'].

    Returns
    -------
    adj_matrix : scipy.sparse.csr_matrix
        Binary adjacency matrix of spatial graph.
    """
    
    
    
    if graph_type=='delaunay_graph':
        spatial_neighbors(adata, coord_type='generic', delaunay=True)
    elif graph_type=='knn_graph':
        samples = adata.obsm['spatial']
        from sklearn.neighbors import NearestNeighbors
        neigh = NearestNeighbors(n_neighbors=2)
        neigh.fit(samples)
        nearest_neighbors=[]
        nearest_neighbor_distances=[]
        for i in samples:
            nearest_neighbor_distances.append(neigh.kneighbors([i])[0][0][1])
            nearest_neighbors.append(neigh.kneighbors([i])[1][0][1])
        G_nearest_neighbors=nx.Graph()
        G_nearest_neighbor_distances=nx.Graph()
        edges_nearest_neighbors=[(list(pd.DataFrame(nearest_neighbors).index)[i], list(pd.DataFrame(nearest_neighbors)[0])[i]) for i in range(0, len(list(pd.DataFrame(nearest_neighbors).index)))]
        edges_nearest_neighbor_distances=[(list(pd.DataFrame(nearest_neighbor_distances).index)[i], list(pd.DataFrame(nearest_neighbor_distances)[0])[i]) for i in range(0, len(list(pd.DataFrame(nearest_neighbor_distances).index)))]
        G_nearest_neighbors.add_edges_from(edges_nearest_neighbors)
        weighted_edges_knn_graph=[(list(pd.DataFrame(nearest_neighbors).index)[i], list(pd.DataFrame(nearest_neighbors)[0])[i], list(pd.DataFrame(nearest_neighbor_distances)[0])[i]) for i in range(0, len(list(pd.DataFrame(nearest_neighbor_distances).index)))]
        G_nearest_neighbor_distances.add_weighted_edges_from(weighted_edges_knn_graph)
        A_nearest_neighbors = nx.to_numpy_array(G_nearest_neighbors)
        A_nearest_neighbors=sparse.csr_matrix(A_nearest_neighbors)
        A_nearest_neighbor_distances = nx.to_numpy_array(G_nearest_neighbor_distances, weight='weight')
        A_nearest_neighbor_distances=sparse.csr_matrix(A_nearest_neighbor_distances)
        adata.obsp['spatial_connectivities']=A_nearest_neighbors
        adata.obsp['spatial_distances']=A_nearest_neighbor_distances
    
    elif graph_type=='dist_thresh_graph':
        X = adata.obsm['spatial']
        Y = adata.obsm['spatial']
        pairwise_distances_matrix=pairwise_distances(X, Y, metric='sqeuclidean')
        # # Pairwise distance matrix filtering for the computation of adjacency matrix:
        # # i. Filtering by absolute value:
        # # pairwise_distances_matrix_filtered = np.where(pairwise_distances_matrix > $a$, 0, pairwise_distances_matrix), where $a$: a real number
        # # ii. Filtering by measures of central tendency (mean, median):
        # # pairwise_distances_matrix_filtered = np.where(pairwise_distances_matrix > $a$, 0, pairwise_distances_matrix), where $a$: {np.mean(pairwise_distances_matrix), np.median(pairwise_distances_matrix)}
        # # iii. Filtering by lowest $a$ percentile of distance values [i.e., each node pair ${n1, n2}$ in the adjacency matrix pairwise_distances_matrix with Euclidean distance $d(n1,n2)$ betwen them smaller than the lowest $a$ percentile of all pairwise distances, shall be retained, all others dropped]:
        # # pairwise_distances_matrix_filtered = np.where(pairwise_distances_matrix > np.percentile(pairwise_distances_matrix, $a$), 0, pairwise_distances_matrix), where $a$: $R \in \[0,100\]$
        A_nearest_neighbor_distance_matrix = np.where(pairwise_distances_matrix > np.percentile(pairwise_distances_matrix, 5), 0, pairwise_distances_matrix)
        A_nearest_neighbor_distances = sparse.csr_matrix(A_nearest_neighbor_distance_matrix)
        A_nearest_neighbor_matrix = np.where(A_nearest_neighbor_distance_matrix > 0, 1, 0)
        A_nearest_neighbors = sparse.csr_matrix(A_nearest_neighbor_matrix)
        adata.obsp['spatial_connectivities']=A_nearest_neighbors
        adata.obsp['spatial_distances']=A_nearest_neighbor_distances
    
    
    adj_matrix, _ = spatial_neighbors(adata, delaunay=True, coord_type='generic', copy=True)
    return adj_matrix


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
