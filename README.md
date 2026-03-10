# Subclone Deconvolution

Input: Variant data from various samples of a bacterial infection that descends from a common ancestor (tsv with variant data and columns: position	depth	ref	alt	DP4	allele_freq	PV4	sample_name)

Output: Identifies subclone populations and their phylogenetic structure


## Clustering into subclones

Use Dirichlet process mixture models to infer clonal population structure, using [pyclone-vi](https://github.com/Roth-Lab/pyclone-vi)

`run_pyclone_iterative.py` 

### Visualise clustering

`pyclone_viz.ipynb`

## Infer phylogenetic structure 

Use [pairtree](https://github.com/morrislab/pairtree) to infer evolutionary relationship between subclones and build a tree, the algorithm conisists of 2 parts

1) Compute pairwise relation tensor over all subclones. This provides the probability over each of four possible evolutionary relationships between each pair: different branches, child, parent, cocluster.

2) Use the pairwise relation tensor to sample possible phylogenies, assigning a likelihood to each. As each phylogeny is sampled, Pairtree computes the subclone frequency of each mutation (or cluster of mutations) within the tree, balancing the need to fit the observed mutation data while still obeying tree constraints.

`run_pairtree.sh`

### Visualise tree corrected subclone frequencies

`pairtree_visualise.ipynb`

## Executables


### `run_pyclone_iterative.py`

Clusters mutations into subclones.

```bash
usage: run_pyclone_iterative.py [-h] [--max_cluster MAX_CLUSTER]
                                [--restarts RESTARTS]
                                [--num_threads NUM_THREADS]
                                [--std_thresh STD_THRESH]
                                [--fixed_thresh FIXED_THRESH]
                                [--neglegible_thresh NEGLEGIBLE_THRESH]
                                [--window_size WINDOW_SIZE]
                                [--recomb_threshold RECOMB_THRESHOLD]
                                pat
```

**Positional arguments**
```
pat
    Patient name
```

**Options**
```
-h, --help
    Show this help message and exit

--max_cluster MAX_CLUSTER
    Maximum number of clusters

--restarts RESTARTS
    Number of restarts for PyClone

--num_threads NUM_THREADS
    Number of threads for PyClone fit

--std_thresh STD_THRESH
    Standard deviation threshold for variability

--fixed_thresh FIXED_THRESH
    Threshold above which a cluster is considered fixed

--neglegible_thresh NEGLEGIBLE_THRESH
    Threshold below which a cluster is considered negligible

--window_size WINDOW_SIZE
    Sliding window size for recombination filtering

--recomb_threshold RECOMB_THRESHOLD
    Threshold for mutations in a window to consider recombination
```
2) `run_pairtree.py`: infers phylogeny

