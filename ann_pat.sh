pat=$1

mapfile -t val_samples < ${pat}/validated_samples.txt


#mvars
for sam in ${val_samples[@]}
do
	vcf="/lustre/scratch127/pam/users/jm73/pseudo/mapping/new_maps_raw/filtered_variants/mvars/${sam}.mvars.vcf"
	awk -F'\t' 'BEGIN{OFS="\t"}{ if($1 == "AE004091") $1 = "Chromosome"; print }' $vcf |
	snpEff ann Pseudomonas_aeruginosa_pao1 -no-intergenic -no-downstream -no-upstream -no-utr -no-intron -s ${pat}/ann_mvars/stats/${sam}_mvars_stats.html |
        cat >> ${pat}/ann_mvars/ann_$(basename $vcf)
done

#cvars
for sam in ${val_samples[@]}
do
        vcf="/lustre/scratch127/pam/teams/team348/jm73/pseudo/analysis/pat_by_pat/${pat}/clean_cvars/${sam}_cvars_clean.vcf"
        awk -F'\t' 'BEGIN{OFS="\t"}{ if($1 == "AE004091") $1 = "Chromosome"; print }' $vcf |
        snpEff ann Pseudomonas_aeruginosa_pao1 -no-intergenic -no-downstream -no-upstream -no-utr -no-intron -s ${pat}/ann_cvars/stats/${sam}_cvars_stats.html |
        cat >> ${pat}/ann_cvars/ann_$(basename $vcf)
done



#for v in ${pat}/results/vcf_filter_out/*.vcf
#do
#	awk -F'\t' 'BEGIN{OFS="\t"}{ if($1 == "AE004091") $1 = "Chromosome"; print }' $v |
#	snpEff ann Pseudomonas_aeruginosa_pao1 -no-intergenic -no-downstream -no-upstream -no-utr -no-intron -s ${pat}/results/vcf_filter_out/stats/$(basename $v)_stats.html |
#	cat >> ${pat}/results/vcf_filter_out/ann_vcf/ann_$(basename $v)
#done
