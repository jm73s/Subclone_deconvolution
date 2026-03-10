import argparse
import pandas as pd
import numpy as np

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Find mutations that are fixed in a patient")
#    parser.add_argument('-v', '--vcf-in', type=str, nargs='+', required=True,
#                        help='1 or more VCF/BCF input files to search')
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='tsv file for a pat containing allele freq for each pos and sample')
    # parser.add_argument('-f', '--file_ext', type=str, required=False, default='csv', choices=['csv', 'tab'],
    #                     help='File extension (default=csv)')
    parser.add_argument('-o', '--output', type=str, required=False, default="fixed_muts.out",
                        help='A single output name for file (default=Mutation_search_VCF.out)')

    args = parser.parse_args()
    input_file=args.input
    output_file=args.output

    pos_fixed=[]
    mut_evo=pd.read_csv(input_file,sep='\t') 
    # for each pos get np array of evolution of allele freq through samples and find those alwas > 0.90
    for pos in mut_evo.position.unique():
        freqs=np.array(mut_evo[mut_evo.position==pos].allele_freq)
        if freqs.min() > 0.90:
            pos_fixed.append(int(pos))

    with open(output_file, 'w') as f:
        for item in pos_fixed:
            f.write(f"AE004091\t{item}\n")
        