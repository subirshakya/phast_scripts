while getopts i:o:n:b: flag; do
	case "${flag}" in
		i) input=${OPTARG};;
		o) output=${OPTARG};;
		n) name=${OPTARG};;
		b) bed=${OPTARG};;
	esac
done

if [ ! -d $input/temp${SLURM_ARRAY_TASK_ID} ]; then
	mkdir -p $input/temp${SLURM_ARRAY_TASK_ID};
fi

if [ ! -d $output ]; then
	mkdir -p $output;
fi

if [ ! -f $input/all_maf_files.txt ]; then
	ls $input/*.maf > $input/all_maf_files.txt
fi

while IFS= read -r line; do
        scaffold=$(echo $line | sed "s/^.*${name}_\([A-Za-z0-9_\.]*\).*/\1/" | sed "s/\.maf//")
        echo $scaffold
        grep -n "^a" $line | cut -f1 -d: > $input/temp_${SLURM_ARRAY_TASK_ID}/$scaffold-alns.txt
        echo "$(paste $input/temp_${SLURM_ARRAY_TASK_ID}/$scaffold-alns.txt <(sed 1,1d $input/temp_${SLURM_ARRAY_TASK_ID}/$scaffold-alns.txt))" > $input/temp_${SLURM_ARRAY_TASK_ID}/$scaffold-alns.txt
        echo "$(awk '{print ($1+1 "     " $2-2)}' $input/temp_${SLURM_ARRAY_TASK_ID}/$scaffold-alns.txt)" > $input/temp_${SLURM_ARRAY_TASK_ID}/$scaffold-alns.txt
        fin_val=$(wc -l < $line)
        fin_val=$((fin_val-1))
        echo "$(sed "s/\-[2]/$fin_val/" $input/temp_${SLURM_ARRAY_TASK_ID}/$scaffold-alns.txt)" > $input/temp_${SLURM_ARRAY_TASK_ID}/$scaffold-alns.txt
        while IFS="     " read -r start stop; do
                stop2=$(($stop+1))
                sed -n "$start,${stop}p;${stop2}q" $line > $input/temp_${SLURM_ARRAY_TASK_ID}/$scaffold-$start-$stop.fa
                pos_start=$(awk 'NR==1{print $3}' $input/temp_${SLURM_ARRAY_TASK_ID}/$scaffold-$start-$stop.fa)
                pos_intv=$(awk 'NR==1{print $4}' $input/temp_${SLURM_ARRAY_TASK_ID}/$scaffold-$start-$stop.fa)
                pos_stop=$(($pos_start+$pos_intv-1))
                sed -i "s/^s\t\([A-za-z_]*\).[A-Za-z0-9_,.\/-]*\t[0-9]*\t[0-9]*\t[+-]\t[0-9]*\t\([ATCG]*\)/>\1\n\2/" temp_${SLURM_ARRAY_TASK_ID}/$scaffold-$start-$stop.fa
                mv $input/temp_${SLURM_ARRAY_TASK_ID}/$scaffold-$start-$stop.fa $input/temp_${SLURM_ARRAY_TASK_ID}/$scaffold-$pos_start-$pos_stop.fa
        done < "$input/temp_${SLURM_ARRAY_TASK_ID}/$scaffold-alns.txt"
        grep -P "${scaffold}\t" $bed > $input/temp_${SLURM_ARRAY_TASK_ID}/$scaffold-cnee.bed
        python ./organize_fasta.py --bed $input/temp_${SLURM_ARRAY_TASK_ID}/$scaffold-cnee.bed --fasta_header $scaffold --out_folder $output --temp_folder $input/temp_${SLURM_ARRAY_TASK_ID}
        #rm $input/temp_${SLURM_ARRAY_TASK_ID}/*.*
        find $input/temp_${SLURM_ARRAY_TASK_ID}/ -type f -delete
done < "$input/all_maf_files.txt"        
