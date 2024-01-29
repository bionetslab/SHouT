

def local_entropy(adata, cluster_key, radius):
    """Computes local entropy scores.

    Parameters
    ----------
    adata : Annotated data object with spatial graph already computed.
    cluster_key : Key of column specifying cell types, states, or clusters.
    radius : Radius to be considered.

    Returns
    -------

    """
    pass


def global_entropy(adata, cluster_key):
    """

    Parameters
    ----------
    adata :
    cluster_key :

    Returns
    -------

    """
    pass


def _compute_local_cell_type_frequencies(adata, cluster_key, radius):
    """

    Parameters
    ----------
    adata :
    cluster_key :
    radius :

    Returns
    -------
    Array with local cell type frequencies for all cells (cells are rows, cell types are columns).

    """


def _compute_global_cell_type_frequencies(adata, cluster_key):
    """

    Parameters
    ----------
    adata :
    cluster_key :

    Returns
    -------
    Dictionary with cell types as key and frequencies as values.

    """