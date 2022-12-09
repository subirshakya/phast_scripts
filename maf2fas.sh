usage()
{
    echo "usage: <command> options:<i|o|n|b|s|j>"
        echo
        echo "i: input folder"
        echo "o: output folder"
        echo "n: filename to append (default = file)"
        echo "b: bed file"
        echo "s: reference species"
        echo "j: number of parallel nodes (defaults = 1)"
}

faster () {
        input=$1
        line=$2
        scaffold=$3
        IFS="     " read -r start stop <<< "$4"
        stop2=$(($stop+1))
        sed -n "$start,${stop}p;${stop2}q" $line > $input/temp_$scaffold/$scaffold-$start-$stop.fa
        pos_start=$(awk 'NR==1{print $3}' $input/temp_$scaffold/$scaffold-$start-$stop.fa)
        pos_intv=$(awk 'NR==1{print $4}' $input/temp_$scaffold/$scaffold-$start-$stop.fa)
        pos_stop=$(($pos_start+$pos_intv-1))
        sed -i "s/^s\t\([A-za-z_0-9]*\).[A-Za-z0-9_,.\/:\-]*\t[0-9]*\t[0-9]*\t[+-]\t[0-9]*\t\([ATCG]*\)/>\1\n\2/" $input/temp_$scaffold/$scaffold-$start-$stop.fa
        mv $input/temp_$scaffold/$scaffold-$start-$stop.fa $input/temp_$scaffold/$scaffold-$pos_start-$pos_stop.fa
        }
export -f faster

name="file"
parallel=1
no_args="true"
while getopts hi:o:n:b:s:j: flag; do
        case "${flag}" in
                i) input=${OPTARG};;
                o) output=${OPTARG};;
                n) name=${OPTARG};;
                b) bed=${OPTARG};;
                s) refsp=${OPTARG};;
                j) parallel=${OPTARG};;
                h) usage
                exit 0;;
                (*)
                        usage
                        exit 0;;
        esac
        no_args="false"
done

[[ "$no_args" == "true" ]] && { usage; exit 1; }

if [ ! -d $output ]; then
        mkdir -p $output;
fi

if [ ! -f $input/all_maf_files.txt ]; then
        ls $input/*.maf > $input/all_maf_files.txt
fi

while IFS= read -r line; do
        scaffold=$(echo $line | sed "s/^.*${name}_\([A-Za-z0-9_\.]*\).*/\1/" | sed "s/\.maf//")
                echo $scaffold
                if [ ! -d $input/temp_$scaffold ]; then
                                mkdir -p $input/temp_$scaffold;
                fi
        grep -n "^a" $line | cut -f1 -d: > $input/temp_$scaffold/$scaffold-alns.txt
        echo "$(paste $input/temp_$scaffold/$scaffold-alns.txt <(sed 1,1d $input/temp_$scaffold/$scaffold-alns.txt))" > $input/temp_$scaffold/$scaffold-alns.txt
        echo "$(awk '{print ($1+1 "     " $2-2)}' $input/temp_$scaffold/$scaffold-alns.txt)" > $input/temp_$scaffold/$scaffold-alns.txt
        fin_val=$(wc -l < $line)
        fin_val=$((fin_val-1))
        echo "$(sed "s/\-[2]/$fin_val/" $input/temp_$scaffold/$scaffold-alns.txt)" > $input/temp_$scaffold/$scaffold-alns.txt
        parallel -j $parallel faster $input $line $scaffold :::: "$input/temp_$scaffold/${scaffold}-alns.txt"
        grep -P "${scaffold}\t" $bed > $input/temp_$scaffold/$scaffold-cnee.bed
        python ./organize_fasta.py --bed $input/temp_$scaffold/$scaffold-cnee.bed --fasta_header $scaffold --out_folder $output --temp_folder $input/temp_$scaffold/ --ref_species $refsp
        rm -r $input/temp_$scaffold
done < "$input/all_maf_files.txt"
