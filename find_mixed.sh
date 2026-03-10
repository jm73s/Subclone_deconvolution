rm mixed.tsv
echo -e "pat\tmixed_timepoints\tmixed" >> mixed.tsv
mapfile -t pats < all_pats.txt
parallel python3 check_mixed.py -p {} -k 300 -w 0.05 -s 5 ::: "${pats[@]}" >> mixed.tsv