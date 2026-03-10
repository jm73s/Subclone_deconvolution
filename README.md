# Subclone Deconvolution

Input: Variant data from various samples of a bacterial infection that descends from a common ancestor (tsv with variant data and columns: position	depth	ref	alt	DP4	allele_freq	PV4	sample_name)

Output: Identifies subclone populations and their phylogenetic structure

## Executables

'''run_pyclone_iterative.py ''': clusters into subclones
'''run_pairtree.py''': infers phylogeny

## Clustering into subclones

Use Dirichlet process mixture models to infer clonal population structure, using pyclone-vi ()

run_pyclone_iterative.py 

### Visualise clustering

pyclone_viz.ipynb

## Infer phylogenetic structure 

Use pairtree to infer evolutionary relationship between subclones and build a tree, the algorithm conisists of 2 parts

1) Compute pairwise relation tensor over all subclones. This provides the probability over each of four possible evolutionary relationships between each pair: different branches, child, parent, cocluster.

2) Use the pairwise relation tensor to sample possible phylogenies, assigning a likelihood to each. As each phylogeny is sampled, Pairtree computes the subclone frequency of each mutation (or cluster of mutations) within the tree, balancing the need to fit the observed mutation data while still obeying tree constraints.

run_pairtree.sh

### Visualise tree corrected subclone frequencies

pairtree_visualise.ipynb

