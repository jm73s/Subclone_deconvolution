# Subclone Deconvolution

Input: Short read illumina data from various samples of a bacterial infection that descends from a common ancestor

Output: Identifies subclone populations and their phylogenetic structure

## Variant calling:

Allign to reference and call variants

## Clustering into subclones

Use Dirichlet process mixture models to infer clonal population structure, using pyclone-vi ()

run_pyclone_iterative.py 

## Infer phylogenetic structure 

Use pairtree to infer evolutionary relationship between subclones and build a tree

run_pairtree.sh


