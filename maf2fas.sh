mkdir -p temp_${SLURM_ARRAY_TASK_ID}
#ls maf/*Sct0y12_*.maf > all_maf_files.txt
while IFS= read -r line; do
        scaffold=$(echo $line | sed 's/maf\/hal25_bulbul_cnee_\(Sct0y12_[0-9]*\).*/\1/')
        echo $scaffold
        grep -n "^a" $line | cut -f1 -d: > temp_${SLURM_ARRAY_TASK_ID}/$scaffold-alns.txt
        echo "$(paste temp_${SLURM_ARRAY_TASK_ID}/$scaffold-alns.txt <(sed 1,1d temp_${SLURM_ARRAY_TASK_ID}/$scaffold-alns.txt))" > temp_${SLURM_ARRAY_TASK_ID}/$scaffold-alns.txt
        echo "$(awk '{print ($1+1 "     " $2-2)}' temp_${SLURM_ARRAY_TASK_ID}/$scaffold-alns.txt)" > temp_${SLURM_ARRAY_TASK_ID}/$scaffold-alns.txt
        fin_val=$(wc -l < $line)
        fin_val=$((fin_val-1))
        echo "$(sed "s/\-[2]/$fin_val/" temp_${SLURM_ARRAY_TASK_ID}/$scaffold-alns.txt)" > temp_${SLURM_ARRAY_TASK_ID}/$scaffold-alns.txt
        while IFS="     " read -r start stop; do
                stop2=$(($stop+1))
                sed -n "$start,${stop}p;${stop2}q" $line > temp_${SLURM_ARRAY_TASK_ID}/$scaffold-$start-$stop.fa
                pos_start=$(awk 'NR==1{print $3}' temp_${SLURM_ARRAY_TASK_ID}/$scaffold-$start-$stop.fa)
                pos_intv=$(awk 'NR==1{print $4}' temp_${SLURM_ARRAY_TASK_ID}/$scaffold-$start-$stop.fa)
                pos_stop=$(($pos_start+$pos_intv-1))
                sed -i "s/^s\t\([A-za-z_]*\).[A-Za-z0-9_,.\/-]*\t[0-9]*\t[0-9]*\t[+-]\t[0-9]*\t\([ATCG]*\)/>\1\n\2/" temp_${SLURM_ARRAY_TASK_ID}/$scaffold-$start-$stop.fa
                mv temp_${SLURM_ARRAY_TASK_ID}/$scaffold-$start-$stop.fa temp_${SLURM_ARRAY_TASK_ID}/$scaffold-$pos_start-$pos_stop.fa
        done < "temp_${SLURM_ARRAY_TASK_ID}/$scaffold-alns.txt"
        grep -P "${scaffold}\t" ../cnee_bulbul_annot_sized.bed > temp_${SLURM_ARRAY_TASK_ID}/$scaffold-cnee.bed
        python ./organize_fasta.py --bed temp_${SLURM_ARRAY_TASK_ID}/$scaffold-cnee.bed --fasta_header $scaffold --out_folder fasta --temp_folder temp_${SLURM_ARRAY_TASK_ID}
        #rm temp_${SLURM_ARRAY_TASK_ID}/*.*
        find temp_${SLURM_ARRAY_TASK_ID}/ -type f -delete
done < "all_maf_split_${SLURM_ARRAY_TASK_ID}"        
