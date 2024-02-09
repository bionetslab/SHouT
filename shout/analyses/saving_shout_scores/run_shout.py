import sys
sys.path.insert(0, '../../../')
from shout import all_scores
# ---
import os
import scanpy as sc
import time
# ---

def read_h5ad_to_pkl(adata_path):
    Pkl={} # Make empty dict.
    for file in os.listdir(adata_path): # Go through all files in adata_path
        filename = os.fsdecode(file)
        if filename.endswith(".h5ad"): # Check to make sure we are only reading .h5ad files.
            Pkl[filename.split('.')[0]]=sc.read_h5ad(adata_path+filename) # or, sc.read()
            continue
        else:
            continue
    return Pkl

def compute_heterogeneity_scores(Pkl, cluster_key, radii, copy):
    patients_good=[]
    patients_bad=[]
    Het_scores={}
    Times={}
    for i in Pkl:
        adata=Pkl[i]
        start = time.time()
        # ---
        try:
            if copy==False:
                all_scores(adata, cluster_key=cluster_key, radii=radii, copy=copy)
                end = time.time()
                adata.uns['SHouT_execution_time']=end-start
            else:
                het_scores=all_scores(adata, cluster_key=cluster_key, radii=radii, copy=copy)
                Het_scores[i]=het_scores
                end = time.time()
                time_elapsed=end-start
                Times[i]=time_elapsed
            print(f'Patient# {i}: Computing local heterogeneity scores...')
            patients_good.append(i)
        except:
            print(f'Patient# {i}: Not enough cells in image to generate spatial neighbors graph.')
            patients_bad.append(i)
            continue
        # ---
        
        
        # ---
        _str_=str(i)+'.h5ad'
        adata.write_h5ad(
        _str_,
        # compression=hdf5plugin.FILTERS["zstd"]
        compression="gzip"
        )
        # ---
    if copy==False:
        return patients_good, patients_bad, None, None
    else:
        return patients_good, patients_bad, Het_scores, Times

def run(adata_path, cluster_key, radii, normalize, num_cell_types, coord_type, copy):
    Pkl=read_h5ad_to_pkl(adata_path) # Read .h5ad files, and save as dict.
    patients_good, patients_bad, het_scores, time_elapsed=compute_heterogeneity_scores(Pkl, cluster_key, radii, copy)
    return patients_good, patients_bad, het_scores, time_elapsed  
    