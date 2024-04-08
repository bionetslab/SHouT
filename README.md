![scAnalyzer](/SHouT.jpeg)

# SHouT
The source code for SHouT (**s**patial **h**eter**o**geneity q**u**antification **t**ool), which aids the analysis of spatial cell graphs with respect to spatial heterogeneity.


SHouT provides two sample-level scores, namely (1) global homophily (edge-based), and (2) global entropy (node-based).


Additionally, SHouT provides three node- (here, cell-) level scores, namely (1) local entropy (node-based), (2) local homophily (edge-based), and (3) local egophily (node-based).






For tutorial, go to .

For cutaneous T-cell lymphoma (CTCL) case study, go to .

## Heterogeneity scores
### I. Global heterogeneity scores

#### a. Global entropy

Global entropy--which is a measure of randomness of the network $G$={$N$, $E$}, where $N$ is the set of all nodes, and $E$ is the set of all edges--is formally defined as:

$$ E_{shannon} = \dfrac{\sum\limits_{k=1}^{K}(n_{c_k})}{N} log_{e} \dfrac{\sum\limits_{k=1}^{K}(n_{c_k})}{N} $$

where, $k=[1, K]$ are the cluster numbers corresponding to clusters $[c_1, c_K]$ present in all nodes within shortest distance $\leq$ radius $r$ from the current node

$n_{c_k}$ is the number of nodes present in cluster $c_k$

$N$: total no. of nodes in neighborhood (all nodes within shortest distance of $\leq$ $r$ from the current node)

**_NOTE:_**  We get one score in the range of [0, 1] per node in the network.


#### b. Global entropy score per network

In case of local entropy, SHout accepts as input a radius (or a set of radii), such that:

$$ E_{shannon} = \dfrac{\sum\limits_{k=1}^{K}(n_{c_k})}{N} log_{e} \dfrac{\sum\limits_{k=1}^{K}(n_{c_k})}{N} $$

where, $k=[1, K]$ are the cluster numbers corresponding to clusters $[c_1, c_K]$ present in all nodes in the network

$n_{c_k}$ is the number of nodes present in cluster $c_k$

$N$: total no. of nodes in the network

**_NOTE:_**  We get one score in the range of [0, 1] for the entire network.








- Agglomerative clustering performed using the [schist nested model](https://schist.readthedocs.io/en/latest/clustering_pbmc.html#clustering-pbmc).
- The schist nested model achieves a more fine-grained clustering as compared to Leiden or schist planted model, achieved over several clustering levels--thereby giving the user higher control over cluster analyses, and how many clusters to stop the clustering algorithm at.
###### Running the clustering algorithm:
Open command prompt/ terminal, then run:
```bash
python3 1_compute_list_of_essential_proteins_for_clustering.py <adata_pickle_path> --user_selected_cluster_level <user_selected_cluster_level>
```

The positional arguments are:
```
[1] adata_pickle_path                   Description: Specify path to anndata pickle data; type=str
```

The optional arguments are:
```
[1] --user_selected_cluster_level       Description: Specify cluster level to get essential proteins list from; type=str [read more below]

{any_integer_value, "complete"} accepted, default is 'complete':
    - If value is integer positive: That represents the cluster level
    - If value is integer negative: That represents the cluster number in reverse order (for example: -1 represents (n-1)-th cluster)
    - If value is 'complete': That represents the last (n-th) cluster. Please note that if your input is "complete", just type it in normally without inverted commas/ quotation marks -- the program is capable of reading strings directly without any quotation symbols.
    - For all other values (i. negative integers; ii. fractional values; iii. any strings other than "complete"; iv. integer value is out of range of [1, n] where n is the total number of clusters): the default is n.
```

The outputs are:
- List of essential proteins per cluster per patient
- Dotplot of potential marker genes vs. cell type clusters (example below for **patient id: 1** of the [TNBC MIBI dataset](https://www.science.org/doi/full/10.1126/sciadv.aax5851))
<br/><br/><br/><br/>
![dotplot_____PatientID-1_____ClusteringLevel-2](readme-images/dotplot_____PatientID-1_____ClusteringLevel-2.png)


#### b. Global entropy:
- In this section, we compute and visualize importances of essential proteins that drive clustering results.

###### Running the protein importance computation:
Open command prompt/ terminal, then run:
```bash
python3 2_visualize_importances_of_essential_proteins_for_clustering.py <adata_pickle_path> <essential_proteins_per_patient_pickle_path>
```

The positional arguments are:
```
[1] adata_pickle_path                                       Description: Specify path to anndata pickle data; type=str
[2] essential_proteins_per_patient_pickle_path              Description: Specify path to essential proteins dictionary generated from the clustering (previous step); type=str
```


The outputs are:
- Feature importance scores calculated using **permutation importance** on **one-vs-all classifier** for **patient id: 1** in the [TNBC MIBI dataset](https://www.science.org/doi/full/10.1126/sciadv.aax5851):
<br/><br/><br/><br/>
![feature_imp_scores_patient1_PermImportance_OneVsAllClassifier](readme-images/feature_imp_scores_patient1_PermImportance_OneVsAllClassifier.png)

- Feature importance scores calculated using **Gini index** on **random forest (RF) classifier** for **patient id: 1** in the [TNBC MIBI dataset](https://www.science.org/doi/full/10.1126/sciadv.aax5851):
<br/><br/><br/><br/>
![feature_imp_scores_patient1_Gini_RFClassifier](readme-images/feature_imp_scores_patient1_Gini_RFClassifier.png)

- Feature importance scores calculated using **permutation importance** on **random forest (RF) classifier** for **patient id: 1** in the [TNBC MIBI dataset](https://www.science.org/doi/full/10.1126/sciadv.aax5851):
<br/><br/><br/><br/>
![feature_imp_scores_patient1_PermImportance_RFClassifier](readme-images/feature_imp_scores_patient1_PermImportance_RFClassifier.png)

### II. Cluster analyses (deciding on the optimal cluster level)
#### i. Maximum bipartite matching between clusterings:
- The similarity scores between clusterings have been obtained for all patients, and their distributions have been plotted against cluster level in the schist nested model.

###### Running the protein importance computation:
Open command prompt/ terminal, then run:
```bash
python3 3_plot_of_robustness_vs_intercluster_similarity.py <adata_pickle_path> --method <method>
```

The positional arguments are:
```
[1] adata_pickle_path                   Description: Specify path to anndata pickle data; type=str
```

The optional arguments are:
```
[1] --method                            Description: Specify cluster agreement method between clusterings; type=str; options={"top-down", "bottom-up"}; default='top-down' [read more below]

Note: The schist hierarchical clustering method is agglomerative, meaning bottom-up.
```

The outputs are:
- Plot of robustness vs. inter-cluster similarity (example below for **patient id: 1** of the [TNBC MIBI dataset](https://www.science.org/doi/full/10.1126/sciadv.aax5851)):
<br/><br/><br/><br/>
(A) Top-down cluster mapping:
<br/><br/><br/><br/>
![Plot_of_robustness_vs_inter_cluster_similarity_method_topDown](readme-images/Plot_of_robustness_vs_inter_cluster_similarity_method_topDown.png)

<br/><br/><br/><br/>
(B) Bottom-up cluster mapping:
<br/><br/><br/><br/>
![Plot_of_robustness_vs_inter_cluster_similarity_method_bottomUp](readme-images/Plot_of_robustness_vs_inter_cluster_similarity_method_bottomUp.png)

<br/><br/><br/><br/>
<br/><br/><br/><br/>

- You can decide on the optimal cluster level based on these plots.
- For our case study on the [TNBC MIBI dataset](https://www.science.org/doi/full/10.1126/sciadv.aax5851) dataset, we will just set optimal cluster level $n = N-1$, where $N$ is the highest cluster level in the schist agglomerative model, where there is just one cluster.





#### ii. Measure no. of clusters per cluster level:

- This may further help some users decide how many cluster levels to set for the schist agglomerative clustering model.

###### Running the computation for no. of clusters per cluster level:
Open command prompt/ terminal, then run:
```bash
python3 4_plot_no_of_cluster_levels_vs_no_of_clusters_per_cluster_level.py <adata_pickle_path> <clusterings_patientLevel_dict_path>
```

The positional arguments are:
```
[1] adata_pickle_path                                                     Description: Specify path to anndata pickle data; type=str
[2] clusterings_patientLevel_dict_path                                    Specify path to clustering combinations calculated in the previous step (saved as 'clusterings_patientLevel_dict.pkl' in step **II.i. Maximum bipartite matching between clusterings**)
```

The outputs are:
- Plot of no. of clusters per cluster level across all samples in the [TNBC MIBI dataset](https://www.science.org/doi/full/10.1126/sciadv.aax5851):
<br/><br/><br/><br/>
![No_of_cluster_levels_vs_No_of_clusters_per_cluster_level](readme-images/No_of_cluster_levels_vs_No_of_clusters_per_cluster_level.png)






### III. Differential analyses (before cell type annotations)
#### i. Protein correlations
- Generate protein correlation matrix.
- Subsequently, perform MWU-test on protein correlation values between conditions.
- Retain correlations with p-values<0.05 as important protein-protein correlations.

###### Obtaining the essential protein correlations:
Open command prompt/ terminal, then run:
```bash
python3 5_DA_protein_correlations.py <adata_pickle_path> <dependent_variable_name>
```

The positional arguments are:
```
[1] adata_pickle_path                   Description: Specify path to anndata pickle data; type=str
[2] dependent_variable_name             Description: Specify the name of the dependent variable; type=str
    
Note: The aforementioned variable <dependent_variable_name> should be present under observation metadata (obsm) of the anndata onject containing gene/ protein expression.
```

The outputs are:
- Bar plot of significant protein correlation p-values
<br/><br/><br/><br/>
![protein_correlation_pvalues_impProteins_allPatients.png](readme-images/protein_correlation_pvalues_impProteins_allPatients.png)



#### ii. Check multiple protein profiles
- Count all cells containing all N-combinaions of proteins present above the threshold.
- Subsequently, perform MWU-test on no. of cells protein profiles expressed in, across conditions.
- Retain protein profiles with p-values<0.05 as significant protein profiles.

###### Obtaining the essential protein correlations:
Open command prompt/ terminal, then run:
```bash
python3 6_DA_multiple_protein_profiles.py <adata_pickle_path> <dependent_variable_name> --N <N> --threshold <threshold>
``` 

The positional arguments are:
```
[1] adata_pickle_path                   Description: Specify path to anndata pickle data; type=str
[2] dependent_variable_name             Description: Specify the name of the dependent variable; type=str
    
Note: The aforementioned variable <dependent_variable_name> should be present under observation metadata (obsm) of the anndata onject containing gene/ protein expression.
```

The optional arguments are:
```
[1] --N                         Description: Specify number of proteins to be analyzed within protein profile; type=int, default=2
[2] --threshold                 Description: Specify minimum value over which protein expressions are considered for further analyses; type=float, default=0.5
```

The outputs are:
- Bar plot of significant, multiple protein coexpression p-values
<br/><br/><br/><br/>
![multi_protein_coexpression_pvalues_impProteins_allPatients.png](readme-images/multi_protein_coexpression_pvalues_impProteins_allPatients.png)
<br/><br/><br/><br/>
Nothing found for N=2, threshold=0.5, across all patients in the [TNBC MIBI dataset](https://www.science.org/doi/full/10.1126/sciadv.aax5851).

##
##
##

# Under development:

### IV. Differential analyses (after cell type annotations)
#### i. Neighborhood enrichment, MWU-test
#### ii. Co-occurence, MWU-test
#### iii. Centrality scores
#### iv. Ripley's statistics
#### v. Spatial autocorrelation
