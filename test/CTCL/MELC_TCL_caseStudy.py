#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 13:00:45 2024

@author: surya
"""

import numpy as np
import pandas as pd
import pickle
import squidpy as sq
import seaborn as sns
import SHouT.shout
from IPython.display import set_matplotlib_formats
import random
import scanpy as sc
import time
import pickle
import matplotlib.pyplot as plt
import hdf5plugin
import os

adata_path='/Users/surya/Documents/GITHUB-REPOSITORIES/SHouT/test/data/'
cluster_key='celltype'
write_adata_path='/Users/surya/Documents/GITHUB-REPOSITORIES/SHouT/test/results/'

Pkl={}
for file in os.listdir(adata_path):
    filename = os.fsdecode(file)
    if filename.endswith(".h5ad"): 
        Pkl[filename.split('.')[0]]=sc.read_h5ad(adata_path+filename) # or, sc.read()
        continue
    else:
        continue

patients_good=[]
patients_bad=[]
for i in Pkl:
    adata=Pkl[i]
    try:
        start = time.time()
        shout.all_scores(adata, cluster_key=cluster_key, radii=[1,2,3,4,5], copy=False)
        end = time.time()
        adata.uns['SHouT_execution_time']=end-start
        # ---
        # ---
        _str_=write_adata_path+str(i)+'.h5ad'
        adata.write_h5ad(
        _str_,
        compression=hdf5plugin.FILTERS["zstd"]
        )
        # ---
        print('Patient#: '+str(i))
        patients_good.append(i)
    except:
        print('Patient# '+str(i)+': Not enough cells in image to generate spatial neighbors graph')
        patients_bad.append(i)
        



















