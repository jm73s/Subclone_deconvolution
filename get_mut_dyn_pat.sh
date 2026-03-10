pat=$1
clean=$2

module load badger/bcftools/1.22
conda activate minority_variant_caller

# get validated samples
mapfile -t val_samples < ${pat}/validated_samples2.txt

mkdir ${pat}/mut_dyn

# get mvars and cvars locations
for i in "${!val_samples[@]}"; do
    if [ "$clean" == "clean" ]; then    
        pat_vars[$((2*i))]="/lustre/scratch127/pam/teams/team348/jm73/pseudo/analysis/pat_by_pat/${pat}/clean_cvars/${val_samples[i]}_cvars_clean.vcf" 
    else
        pat_vars[$((2*i))]="/lustre/scratch127/pam/teams/team348/jm73/pseudo/mapping/new_maps_raw/filtered_variants/cvars/${val_samples[i]}.cvars.vcf"
    fi
    pat_vars[$((2*i + 1))]="/lustre/scratch127/pam/teams/team348/jm73/pseudo/mapping/new_maps_raw/filtered_variants/mvars/${val_samples[i]}.mvars.vcf"
done

## get common positions in filtered vcf
rm ${pat}/mut_dyn/all_pos.txt
get_all_mutation_positions.py -v ${pat_vars[@]} -o ${pat}/mut_dyn/all_pos.txt

## get raw vcfs
for i in "${!val_samples[@]}"; do
                pat_vars_raw[i]="/lustre/scratch127/pam/teams/team348/jm73/pseudo/mapping/new_maps_raw/${val_samples[i]}/vcf/${val_samples[i]}.vcf.gz"
done
## get evolution of vars through sample
if [ "$clean" == "clean" ]; then
    rm ${pat}/mut_dyn/mut_evolution_clean.txt
    Mutation_search_VCF.py -v ${pat_vars_raw[@]} -i ${pat}/mut_dyn/all_pos.txt -o ${pat}/mut_dyn/mut_evolution_clean.txt
else
    rm ${pat}/mut_dyn/mut_evolution.txt
    Mutation_search_VCF.py -v ${pat_vars_raw[@]} -i ${pat}/mut_dyn/all_pos.txt -o ${pat}/mut_dyn/mut_evolution.txt
fi

# add dates 
if [ "$clean" == "clean" ]; then
    rm ${pat}/mut_dyn/mut_evolution_clean_dates.txt
    python3 add_dates_to_mut_dyn.py -i ${pat}/mut_dyn/mut_evolution_clean.txt -o ${pat}/mut_dyn/mut_evolution_clean_dates.txt
else
    rm ${pat}/mut_dyn/mut_evolution_dates.txt
    python3 add_dates_to_mut_dyn.py -i ${pat}/mut_dyn/mut_evolution.txt -o ${pat}/mut_dyn/mut_evolution_dates.txt
fi