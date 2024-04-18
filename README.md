![scAnalyzer](/SHouT.jpeg)

# SHouT
The source code for SHouT (**s**patial **h**eter**o**geneity q**u**antification **t**ool), which aids the analysis of spatial cell graphs with respect to spatial heterogeneity.


SHouT provides two sample-level scores, namely (1) global homophily (edge-based), and (2) global entropy (node-based).


Additionally, SHouT provides three node- (here, cell-) level scores, namely (1) local entropy (node-based), (2) local homophily (edge-based), and (3) egophily (node-based).


For tutorial, go to https://github.com/bionetslab/SHouT/blob/main/tutorial.ipynb.

For cutaneous T-cell lymphoma (CTCL) case study, go to .

## Heterogeneity scores

- SHouT starts by computing sample-specific spatial neighborhood graphs $G=(V,E,\lambda_V)$ from the pre-processed imaging data, where $V\subseteq\mathcal{C}$ is the set of cells for the sample under consideration, the set $E$ contains an edge $cc^\prime$ for two cells $c,c^\prime\in V$ if $c$ and $c^\prime$ are spatially adjacent (computed with Squidpy's \texttt{spatial\_neighbors} function with the parameter \texttt{delaunay} set to \texttt{True}), and $\lambda_V$ denotes the restriction of the cell type label function $\lambda$ to $V$.

- Subsequently, the two global scores, and three node-specific scores are calculated.

### I. Global heterogeneity scores

#### a. Global entropy

Global entropy is formally defined as:

$$ H(G)=-\log(|T|)^{-1}\cdot\sum_{t\in T}p_G(t)\cdot\log(p_G(t))\in[0,1]\text{,} $$

where $p_G(t)=|\{c\in V\mid \lambda_V(c)=t\}|/|V|$ is the fraction of cells in $V$ that are of type $t$. Large values of $H(G)$ indicate that cell type heterogeneity is high for the sample represented by $G$. The second global score\,---\,global homophily\,---\,is defined as the fraction of edges.


#### b. Global homophily

Global homophily is formally defined as:

$$ h(G)=|E|^{-1}\cdot\sum_{cc^\prime\in E}[\lambda_V(c)=\lambda_V(c^\prime)]\in[0,1] $$

where, $h(G)$ is the fraction of edges in the spatial graph $G$ that connects cells of the same type $([.]: {True, False} \rightarrow {0,1}$ is the Iverson bracket, i.e., $[True]=1$ and $[False]=0)$. Large values of $h(G)$ indicate that cells tend to be adjacent to cells of the same type in the sample represented by $G$.

**_NOTE:_**  We get one global score in the range of [0, 1] per radius for the entire network.



### I. Local heterogeneity scores

- In case of local heterogeneity scores, SHout accepts as input a radius $r$ (or a set of radii).

- It then calculates local scores within the $r$-hop neighborhood $N_r(c)=\{c^\prime\in V\mid d_G(c,c^\prime\leq r\}$ of an individual cell $c\in V$. Here, $r$ is a hyper-parameter and $d_G:V\times V\to\mathbb{N}$ is the shortest path distance.

#### a. Local entropy

Local entropy is formally defined as:

$$ H_r(c)=-\log(|T|)^{-1}\cdot\sum_{t\in T}p_{N_r(c)}(t)\cdot\log(p_{N_r(c)}(t))\in[0,1] $$

Here, as opposed to global entropy, cell type fractions $p_{N_r(c)}(t)=|\{c\in N_r(c)\mid \lambda_V(c)=t\}|/|N_r(c)|$ are computed only with respect to the cells contained in the $r$-hop neighborhood of $c$.


#### b. Local homophily

Local homophily is formally defined as:

$$ h_r(c)=|E_{N_r(c)}|^{-1}\cdot\sum_{c^\prime c^{\prime\prime}\in E_{N_r(c)}}[\lambda_V(c^\prime)=\lambda_V(c^{\prime\prime})]\in[0,1] $$

Here, as opposed to global homophily, we only consider the subset of edges $E_{N_r(c)}=\{c^\prime c^{\prime\prime}\in E\mid c^\prime,c^{\prime\prime}\in N_r(c)\}$ that connect two cells contained in the $r$-hop neighborhood of $c$.


#### c. Egophily

Egophily is defined as:

$$e_r(c)=p_{N_r(c)}(\lambda_V(c))\in[0,1]$$

where, egophily $e_r(c)$ is the fraction of cells within the $r$-hop neighborhood of $c$ that have the same cell type as $c$.



## Installation

Install conda environment as follows (there also exists a requirements.txt)
```bash
conda create --name shout
conda activate shout
pip install scipy==1.10.1 numpy==1.23.5 squidpy==1.3.0
```


## Cutaneous T cell Lymphoma (CTCL) case study
The CTCL case study can be accessed at: [https://github.com/bionetslab/ctcl_case_study](https://github.com/bionetslab/ctcl_case_study)


## Citing the work

Please cite the paper as follows:
- Will be updated once available.


















