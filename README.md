# Subclone Deconvolution

Input: Short read illumina data from various samples of a bacterial infection that descends from a common ancestor

Output: Identifies subclone populations and their phylogenetic structure

## Variant calling:

Allign to reference and call variants

## Clustering into subclones

Use Dirichlet process mixture models to infer clonal population structure, using pyclone-vi ()

run_pyclone_iterative.py 

## Infer phylogenetic structure 

Use pairtree to infer evolutionary relationship between subclones and build a tree, the algorithm conisists of 2 parts

1) Compute pairwise relation tensor over all subclones. This provides the probability over each of four possible evolutionary relationships between each pair: different branches, child, parent, cocluster.

2) Use the pairwise relation tensor to sample possible phylogenies, assigning a likelihood to each. As each phylogeny is sampled, Pairtree computes the subclone frequency of each mutation (or cluster of mutations) within the tree, balancing the need to fit the observed mutation data while still obeying tree constraints.

run_pairtree.sh


