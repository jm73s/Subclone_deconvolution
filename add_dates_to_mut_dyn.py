import argparse
import pandas as pd
import numpy as np

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Add dates to mutation dynamics tsv file")
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='tsv file with mutation dynamics')
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='output name for file')
    
    args = parser.parse_args()
    input_file=args.input
    output_file=args.output

    # read in pat_sample.txt
    pat_sample=pd.read_csv("/lustre/scratch127/pam/teams/team348/jm73/pseudo/sequencing_progress/pat_sample.txt",sep='\t')
    # read in mut dyn file
    mut_dyn=pd.read_csv(input_file,sep='\t')
    mut_dyn['date']=mut_dyn['lane'].map(lambda x: pat_sample[pat_sample['ID']==x].Date.values[0] if x in pat_sample.ID.values else np.nan)
    mut_dyn.to_csv(output_file,sep='\t',index=False)
