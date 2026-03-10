#conda activate pyclone-vi

pat=$1


rm ${pat}/mut_dyn/pyclone/pyclone_out_${pat}_filtered.*


pyclone-vi fit --seed 15 -i ${pat}/mut_dyn/pyclone/pyclone_in_${pat}_filtered.tsv -o ${pat}/mut_dyn/pyclone/pyclone_out_${pat}_filtered.h5 -c 40 -d beta-binomial -r 50
pyclone-vi write-results-file -i ${pat}/mut_dyn/pyclone/pyclone_out_${pat}_filtered.h5 -o ${pat}/mut_dyn/pyclone/pyclone_out_${pat}_filtered.tsv