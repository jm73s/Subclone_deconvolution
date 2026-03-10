#!/usr/bin/env python3

import pandas as pd
import subprocess
import numpy as np
import argparse
import time

# run pyclone iteratively for a patient until there are no more mutations or clusters to remove

def run_pyclone(pat,max_cluster, restarts, iteration, num_threads=20, seed=15):
    """
    Run pyclone for a patient with the given parameters.
    input: pat: patient name, max_cluster: maximum number of clusters, restarts: number of restarts for pyclone,
           num_threads: number of threads for pyclone, seed: random seed for pyclone
    output: pyclone output file in h5 and tsv format in pyclone pat folder
    iteration: current iteration number to use in file naming
    """
    subprocess.run(
        [
            "pyclone-vi", "fit",
            "--seed", str(seed),
            "-i", f"{pat}/mut_dyn/pyclone2/pyclone_in_{pat}_iter{iteration}.tsv",
            "-o", f"{pat}/mut_dyn/pyclone2pyclone_out_{pat}_iter{iteration}.h5",
            "-c", str(max_cluster),
            "-d", "beta-binomial",
            "-r", str(restarts),
            "-t", str(num_threads)
        ],
        check=True
    )
    subprocess.run(
        [
            "pyclone-vi", "write-results-file",
            "-i", f"{pat}/mut_dyn/pyclone2/pyclone_out_{pat}_iter{iteration}.h5",
            "-o", f"{pat}/mut_dyn/pyclone2/pyclone_out_{pat}_iter{iteration}.tsv"
        ],
        check=True
    )


# remove fixed and too low clusters
def cluster_filter(pat, std_thresh, fixed_thresh, neglegible_thresh, var_thresh,iteration, true_mut_df):
    '''
    Clean clusters by removing those that are fixed (probably differencence with pao1) and those that are mostly negligible (noise or contamination)
    after ignoring timepoints with high variability.
    pat: patient to consider
    std: threshold for variability to ignore timepoints
    fixed_thresh: threshold above which a cluster is considered fixed
    neglegible_thresh: threshold below which a cluster is considered negligible if only goes above once or less
    iteration: current iteration number 
    true_mut_df: dataframe with true mutation dynamics for the patient
    '''


    #load pyclone df
    mut_df=pd.read_csv(f'{pat}/mut_dyn/pyclone2/pyclone_out_{pat}_iter{iteration}.tsv',sep='\t')
    mut_df_out=mut_df.copy() # for output later
    # add sample dates
    mut_df['date'] = pd.to_datetime(mut_df['sample_id'].str.split('_').str[0], format='%d-%m-%Y', errors='coerce')
    mut_df = mut_df.dropna(subset=['date'])
    mut_df.sort_values(by='date',inplace=True)

    #assign clusters to true mut
    true_mut_df = true_mut_df.merge(
    mut_df[['mutation_id','cluster_id']].drop_duplicates(),
    left_on='position',
    right_on='mutation_id',
    how='left'
    ).drop('position', axis=1)
    true_mut_df['cluster_id']=true_mut_df['cluster_id'].astype('Int64') # allow NaN
    # remove non clustered mutations
    true_mut_df=true_mut_df[true_mut_df['cluster_id'].notnull()]
    # add predicted pyclone cellular prevalence
    true_mut_df=pd.merge(true_mut_df,mut_df[['mutation_id','date','cellular_prevalence']],on=['mutation_id','date'],how='inner')


    # calc std for each cluster at each timepoint from true mut df
    std = true_mut_df.groupby(['cluster_id','date'])['allele_freq'].std().reset_index()
    std.rename(columns={'allele_freq':'std'},inplace=True)
    # add to true mut df
    true_mut_df_clean=pd.merge(true_mut_df,std,on=['cluster_id','date'],how='inner')
    # remove entries with std bigger than threshold
    true_mut_df_clean=true_mut_df_clean[true_mut_df_clean['std']<std_thresh]
    # find clusters where predicted cellular prevalence never bellow fixed_thresh after removing high std entries
    #fixed_clusters=true_mut_df_clean.groupby('cluster_id').apply(lambda x: (x['cellular_prevalence'] < fixed_thresh).any() == False)
    # find clusters where predicted cellular prevalence falls at most once bellow fixed_thresh after removing high std entries
    fixed_clusters=true_mut_df_clean.groupby('cluster_id').apply(lambda x: (x['cellular_prevalence'] < fixed_thresh).sum() < 2)
    # find clusters where predicted cellular prevalence only above neglegible_thresh once after removing high std entries
    neglegible_clusters=true_mut_df_clean.groupby('cluster_id').apply(lambda x: (x['cellular_prevalence'] > neglegible_thresh).sum() < 2)
    # find clusters where predicted cellular prevalence never further away than var_thresh from mean more than once, after removing high std entries
    mean=true_mut_df_clean.groupby('cluster_id')['cellular_prevalence'].mean()
    constant_clusters = true_mut_df_clean.groupby('cluster_id').apply(lambda x: (
    ((x['cellular_prevalence'] > mean.loc[x.name] + var_thresh) |
     (x['cellular_prevalence'] < mean.loc[x.name] - var_thresh))
    .sum() < 2))
    
    remove_clusters=fixed_clusters[fixed_clusters].index.tolist()+ neglegible_clusters[neglegible_clusters].index.tolist()+ constant_clusters[constant_clusters].index.tolist()

    # output cleaned pyclone data
    # remove clusters
    mut_df_out=mut_df_out[~mut_df_out['cluster_id'].isin(remove_clusters)]
    # relable clusters to account for removed ones and align with pairtree numbering
    # get sorted unique cluster labels
    clusters_sorted = sorted(mut_df_out['cluster_id'].unique())
    # build mapping: old_label -> new_label (1..12)
    mapping = {old: new for new, old in enumerate(clusters_sorted, start=1)}
    # relabel
    mut_df_out['cluster_id'] = mut_df_out['cluster_id'].map(mapping)
    # write to file
    with open(f'{pat}/mut_dyn/pyclone2/pyclone_out_{pat}_iter{iteration}_cluster.tsv','w') as f:
            mut_df_out.to_csv(f,sep='\t',index=False)

    return sorted(remove_clusters)

def sliding_recomb_filter(pat, window_size, threshold,iteration, skip=1):

    # load initial pyclone df
    mut_df=pd.read_csv(f'{pat}/mut_dyn/pyclone2/pyclone_out_{pat}_iter{iteration}_cluster.tsv',sep='\t')
    mut_df['date']=mut_df['sample_id'].str.split('_').str[0]
    mut_df['date']=pd.to_datetime(mut_df['date'],format='%d-%m-%Y', errors='coerce')
    mut_df = mut_df.dropna(subset=['date'])
    
    recomb_pos=set()
    for cluster_id in mut_df['cluster_id'].unique():

        cluster_df=mut_df[mut_df['cluster_id']==cluster_id]
        cluster_df['position']=cluster_df['mutation_id'].str.split('_',expand=True)[1].astype(int)
        # get pos for cluster
        c_pos=cluster_df['position'].unique() 

        # array with length up to biggest pos, with true where mut and false where not
        numline=np.array((c_pos.max()+1)*[False])
        numline[c_pos]=True


        # sliding window across numline and count number of mutations in each window, if above threshold, add to recomb_pos
        for i in range(1, c_pos.max()-window_size+1, skip):
            count=numline[i:i+window_size].sum()
            if count > threshold:
                recomb_pos.update(np.where(numline[i:i+window_size])[0]+i)
        
    # output cleaned pyclone data
    # remove muts in recombination regions
    remove_muts=[f'pos_{pos}' for pos in list(recomb_pos)]
    mut_df=mut_df[~mut_df['mutation_id'].isin(remove_muts)]
    # write to file
    with open(f'{pat}/mut_dyn/pyclone2/pyclone_out_{pat}_iter{iteration}_cluster_recomb.tsv','w') as f:
            mut_df.to_csv(f,sep='\t',index=False)

    return recomb_pos

def gen_pyclone_df(pat, output_name, iteration):

    # load mutation dynamics df
    pat_dyn = pd.read_csv(f'{pat}/mut_dyn/mut_evolution_dates.txt',sep='\t')
    pat_dyn['date_lane']=pat_dyn['date'].str.cat(pat_dyn['lane'], sep='_')
    # remove positions where allele_freq is always over 0.9 
    #always_high = pat_dyn.groupby('position')['allele_freq'].transform(lambda x: (x > 0.9).all())
    # remove positions where allele_freq is always over 0.9 (excpet maybe at 1 pos)
    always_high = pat_dyn.groupby('position')['allele_freq'].transform(lambda x: (x < 0.9).sum() < 2)
    # remove positions where allele_freq is above 5% 1 time or less
    only_once_over_5 = pat_dyn.groupby('position')['allele_freq'].transform(lambda x: (x > 0.05).sum() < 2)
    pat_dyn = pat_dyn.loc[~always_high & ~only_once_over_5].copy()
    df_pyclone=pat_dyn[['position','date_lane','DP4']].copy()
    rename_dict={'position':'mutation_id','date_lane':'sample_id'}
    df_pyclone.rename(columns=rename_dict,inplace=True)
    df_pyclone['mutation_id']='pos_'+df_pyclone['mutation_id'].astype(str)
    df_pyclone['ref_counts']=df_pyclone['DP4'].str.split(',', expand=True)[[0,1]].astype(int).sum(axis=1)
    df_pyclone['alt_counts']=df_pyclone['DP4'].str.split(',', expand=True)[[2,3]].astype(int).sum(axis=1)
    df_pyclone['major_cn']=1
    df_pyclone['minor_cn']=0
    df_pyclone['normal_cn']=1
    df_pyclone['tumour_content']=1
    df_pyclone['error_rate']=0.0003 # phred 35
    df_pyclone.drop(columns='DP4',inplace=True)

    # if using existing pyclone output, extract mutation_id column to filter df_pyclone to only those mutation_ids and apply recomb filter
    if iteration>0:
        usecols=['mutation_id']
        df_pyclone_old = pd.read_csv(f'{pat}/mut_dyn/pyclone2/pyclone_out_{pat}_iter{iteration-1}_cluster_recomb.tsv', sep='\t',usecols=usecols) # skip first row with cluster info if needed
        mutation_ids=df_pyclone_old['mutation_id'].unique()
        df_pyclone = df_pyclone[df_pyclone.mutation_id.isin(mutation_ids)].copy() # filter df_pyclone to only mutation_ids in existing pyclone output
    
    df_pyclone.to_csv(f'{pat}/mut_dyn/pyclone2/{output_name}',sep='\t',index=False)


def iterative_pyclone(pat, max_cluster, restarts, num_threads, std_thresh, fixed_thresh, neglegible_thresh, var_thresh,window_size, recomb_threshold):
    """
    Run pyclone iteratively for a patient until there are no more mutations or clusters to remove.
    After each run, filter clusters and optionally remove those with high recombination.
    pat: patient name
    max_cluster: maximum number of clusters for pyclone
    restarts: number of restarts for pyclone
    num_threads: number of threads for pyclone fit
    std_thresh: threshold for variability to ignore timepoints in cluster filtering
    fixed_thresh: threshold above which a cluster is considered fixed in cluster filtering
    neglegible_thresh: threshold below which a cluster is considered negligible if only goes above once or less in cluster filtering
    var_thresh: Threshold below which a cluster is considered constant
    window_size: size of sliding window for recombination filtering
    recomb_threshold: number of mutations in window above which it is considered a recombination region
    """
    iteration=0


    # load called snp df for the patient to use in cluster filtering
    true_mut_df=pd.read_csv(f'{pat}/mut_dyn/mut_evolution_dates.txt',sep='\t')
    true_mut_df['date']=pd.to_datetime(true_mut_df['date'],format='%d-%m-%Y')
    true_mut_df['position']='pos_'+true_mut_df['position'].astype(str)


    while True:
        print(f'Iteration {iteration}')

        # generate pyclone input file from mutation dynamics and existing pyclone output, then run pyclone
        print('Running PyClone...')
        gen_pyclone_df(pat, f'pyclone_in_{pat}_iter{iteration}.tsv', iteration)
        start_time = time.time()
        run_pyclone(pat,max_cluster,restarts,iteration,num_threads=num_threads)
        end_time = time.time()
        print(f"Execution time for PyClone: {end_time - start_time:.2f} seconds")

        print('Filtering clusters and recombination regions...')
        start_time = time.time()
        # remove fixed and negligible clusters, then remove clusters with high recombination
        removed_clusters=cluster_filter(pat,std_thresh,fixed_thresh,neglegible_thresh, var_thresh,iteration,true_mut_df)
        print(f'Removed clusters: {removed_clusters}')
        end_time = time.time()
        print(f"Execution time for cluster_filter: {end_time - start_time:.2f} seconds")
        
        start_time = time.time()
        recomb_pos=sliding_recomb_filter(pat,window_size,recomb_threshold,iteration)
        print(f'Positions in recombination regions: {len(recomb_pos)}')
        end_time = time.time()
        print(f"Execution time for recomb_filter: {end_time - start_time:.2f} seconds")
        
        #  if no clusters or positions removed, stop iterating
        if len(removed_clusters)==0 and len(recomb_pos)==0:
            break

        iteration+=1

def main():
    parser = argparse.ArgumentParser(description="Run PyClone iteratively for a patient.")
    parser.add_argument("pat", type=str, help="Patient name")
    parser.add_argument("--max_cluster", type=int, default=40, help="Maximum number of clusters")
    parser.add_argument("--restarts", type=int, default=50, help="Number of restarts for PyClone")
    parser.add_argument("--num_threads", type=int, default=20, help="Number of threads for PyClone fit")
    parser.add_argument("--std_thresh", type=float, default=0.05, help="Standard deviation threshold for variability")
    parser.add_argument("--fixed_thresh", type=float, default=0.9, help="Threshold above which a cluster is considered fixed")
    parser.add_argument("--neglegible_thresh", type=float, default=0.05, help="Threshold below which a cluster is considered negligible")
    parser.add_argument("--var_thresh", type=float, default=0.025, help="Threshold below which a cluster is considered constant")
    parser.add_argument("--window_size", type=int, default=500, help="Sliding window size for recombination filtering")
    parser.add_argument("--recomb_threshold", type=int, default=5, help="Threshold for mutations in a window to consider recombination")
    
    args = parser.parse_args()
    
    iterative_pyclone(
        pat=args.pat,
        max_cluster=args.max_cluster,
        restarts=args.restarts,
        num_threads=args.num_threads,
        std_thresh=args.std_thresh,
        fixed_thresh=args.fixed_thresh,
        neglegible_thresh=args.neglegible_thresh,
        var_thresh=args.var_thresh,
        window_size=args.window_size,
        recomb_threshold=args.recomb_threshold
    )
    
if __name__ == "__main__":
    main()
