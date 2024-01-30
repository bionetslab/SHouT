#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 20:14:14 2024

@author: surya
"""

import squidpy as sq
import numpy as np
from scipy.sparse.csgraph import shortest_path
import pickle
# ---
from shout import utility, homophily
# ---

adata_pickle_path='/Volumes/time-mach-14pro/adata_cell_celltypesAnnotated_TCL_vs_PSO.pkl'
# ---
with open(adata_pickle_path, 'rb') as f:
    pickle_=pickle.load(f)
# ---
adata=pickle_[81]
# ---
adj_matrix, adj_matrix_homophilic, shortest_path_distances=utility.get_spatial_graph(adata, 'celltype', 'generic')
radius=2
extended_neighborhoods=utility.get_extended_neighborhoods(shortest_path_distances, radius)
# ---

# ===== Global homophily =====
homophily.global_homophily(adata, adj_matrix, adj_matrix_homophilic)

# ===== Local homophily =====
homophily.local_homophily(adata, shortest_path_distances, radius, adj_matrix, adj_matrix_homophilic)

# =============================================


