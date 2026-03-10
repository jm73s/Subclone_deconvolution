import subprocess


def build_pairtree_inputs(pat,iter):
    def compress_position(cols,df): # compress all cols for each position into a single string
        grouped = df.groupby('position')[cols].apply(lambda x: pd.Series({col: ','.join(x[col].astype(str)) for col in cols}), include_groups=False)
        return grouped.reset_index()



    # load snp df
    df_snp=pd.read_csv(f'{pat}/mut_dyn/mut_evolution_dates.txt',sep='\t',usecols=['position','depth','DP4','lane','date','allele_freq'])
    df_snp['lane']=df_snp['lane']+'_'+df_snp['date']
    # filter snps
    # only keep snps in pyclone output
    pyclone_df=pd.read_csv(f'{pat}/mut_dyn/pyclone/pyclone_out_{pat}_iter{iter}.tsv',sep='\t',usecols=['mutation_id','cluster_id'],skiprows=1)
    pyclone_pos=pyclone_df['mutation_id'].map(lambda x: int(x.split('_')[1])).unique()
    df_snp=df_snp[df_snp['position'].isin(pyclone_pos)].copy()
    ## remove positions where allele_freq is always over 0.9
    #always_high = df_snp.groupby('position')['allele_freq'].transform(lambda x: (x > 0.9).all())
    ## remove positions where allele_freq is above 5% 1 time or less
    #only_once_over_5 = df_snp.groupby('position')['allele_freq'].transform(lambda x: (x > 0.05).sum() < 2)
    #df_snp = df_snp.loc[~always_high & ~only_once_over_5].copy()
    ## generate new info

    #df_snp['var_reads']=df_snp['DP4'].map(lambda x: int(x.split(',')[2])+int(x.split(',')[3]))
    #df_snp['var_read_prob']=1
    #df_snp.rename(columns={'depth':'total_reads'},inplace=True)
    # vectorized DP4 parsing -> var_reads
    dp4 = df_snp['DP4'].str.split(',', expand=True)
    dp4 = dp4.astype(int)
    df_snp['var_reads'] = dp4.iloc[:, 2] + dp4.iloc[:, 3]
    df_snp['var_read_prob'] = 1
    df_snp.rename(columns={'depth': 'total_reads'}, inplace=True)

    # create ssm
    # sort so that we have the same lane order for each position
    df_snp.sort_values(['position','lane'],inplace=True)
    
    # create new df
    # compress columns for each pos
    df_ssm=compress_position(['var_reads','var_read_prob','total_reads'], df_snp)
    df_ssm.rename(columns={'position':'name'},inplace=True)

    # create name col and edit id col
    df_ssm['name']=df_ssm['name'].map(lambda x: f'pos_{x}')
    df_ssm['id']=df_ssm.index.map(lambda x: f's{x}') # need unique ids for pairtree s0,s1,...
    df_ssm.to_csv(f'{pat}/mut_dyn/pairtree/ssm_{pat}.tsv',sep='\t',index=False)

    # create json
    json_dict={}
    # add samples to json
    json_dict['samples']=df_snp['lane'].unique().tolist()
    # add clusters to json
    pyclone_df.drop_duplicates(subset='mutation_id', inplace=True)
    name_to_id = dict(zip(df_ssm['name'], df_ssm['id']))
    pyclone_df['pairtree_id'] = pyclone_df['mutation_id'].map(name_to_id)
#    pyclone_df['pairtree_id']=pyclone_df['mutation_id'].map(lambda x: df_ssm[df_ssm['name']==x]['id'].values[0])
    clusters=[]
    for cluster_id in sorted(pyclone_df['cluster_id'].unique()):
        clusters.append(pyclone_df[pyclone_df.cluster_id==cluster_id]['pairtree_id'].tolist())
    json_dict['clusters']=clusters
    json_dict['garbage']= []
    with open(f'{pat}/mut_dyn/pairtree/json_{pat}.json','w') as f:
        json.dump(json_dict,f,indent=2)

def run_pairtree(pat):
        subprocess.run(["bash", "pairtree.sh", pat],check=True)
    