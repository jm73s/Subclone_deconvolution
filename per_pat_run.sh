#! /bin/bash
#BSUB -o run_all_pats.o
#BSUB -e run_all_pats.e
#BSUB -q long
#BSUB -n 10
#BSUB -M60000
#BSUB -R "select[mem>60000] rusage[mem=60000] span[hosts=1]"

#conda activate minority_variant_caller
#module load badger/bcftools/1.21 
conda activate pyclone-vi

script=run_pyclone.sh

## get list of pats with validated samples
#rm all_pats.txt
#for pat in *_*
#do
#	if test -f ${pat}/validated_samples2.txt; then
#		echo ${pat} >> all_pats.txt
#	fi
#done

parallel -j 4 "bash ${script} {} " :::: all_pats_left.txt



