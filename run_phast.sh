wget https://cgl.gi.ucsc.edu/data/cactus/363-avian-2020.hal

#phastcons
#get cds
grep -P "\tCDS\t" P_atriceps_annot_sorted_fixed.gff > P_atriceps_cds.gff
gff3ToGenePred -useName P_atriceps_cds.gff P_atriceps_cds.gp
genePredToBed P_atriceps_cds.gp P_atriceps_cds.bed

#run as array of 5 files to get 4d sites from CDS
split -l 24000 --numeric-suffixes=1 ../../P_atriceps_cds.bed Patri_cds_
hal4dExtract --conserved 363-avian-2020.hal Brachypodius_atriceps ./4d/beds/Patri_cds_0${SLURM_ARRAY_TASK_ID} ./4d/4d_beds/Patri_4d_0${SLURM_ARRAY_TASK_ID}.bed

#join beds to single bed
cat *.bed > Patri_4d.bed

#get 4D sites in maf format for 25 genomes from hal file using bedfile
hal2mafMP.py  --numProc 6 --refGenome Brachypodius_atriceps --noDupes --noAncestors --refTargets Patri_4d.bed --targetGenomes Brachypodius_atriceps,Pycnonotus_jocosus,Sylvia_atricapilla,Zosterops_lateralis,Phylloscopus_trochilus,Hirundo_rustica,Sturnus_vulgaris,Taeniopygia_guttata,Geospiza_fortis,Corvus_brachyrhynchos,Manacus_manacus,Smithornis_capensis,Acanthisitta_chloris,Melopsittacus_undulatus,Halcyon_senegalensis,Aquila_chrysaetos,Phalacrocorax_harrisi,Nipponia_nippon,Aptenodytes_forsteri,Cuculus_canorus,Caloenas_nicobarica,Gallus_gallus,Eudromia_elegans,Dromaius_novaehollandiae,Cairina_moschata ../../363-avian-2020.hal /hal25_bulbul_4d.maf

#convert to ss
msa_view --in-format MAF --out-format SS --unordered-ss 4d/hal25_bulbul_4d.maf > 4d/hal25_bulbul_4d.ss

#phylo_models to get neutral rate
phyloFit --tree "(((Caloenas_nicobarica,(Cuculus_canorus,(((Halcyon_senegalensis,Aquila_chrysaetos),(Melopsittacus_undulatus,(Acanthisitta_chloris,((Corvus_brachyrhynchos,((Hirundo_rustica,(Phylloscopus_trochilus,((Sylvia_atricapilla,Zosterops_lateralis),(Pycnonotus_jocosus,Brachypodius_atriceps)))),(Sturnus_vulgaris,(Taeniopygia_guttata,Geospiza_fortis)))),(Manacus_manacus,Smithornis_capensis))))),((Phalacrocorax_harrisi,Nipponia_nippon),Aptenodytes_forsteri)))),(Cairina_moschata,Gallus_gallus)),(Dromaius_novaehollandiae,Eudromia_elegans))" --init-random --subst-mod SSREV --sym-freqs --log bulbul_neut.log --msa-format SS --out-root bulbul_nonconserved-4d  hal25_bulbul_4d.ss

#this will give a .mod file with rates for neutral processes

#get maf of 25 genomes from hal file for entire alignment
hal2mafMP.py  --numProc 6 --refGenome Brachypodius_atriceps --noDupes --noAncestors --targetGenomes Brachypodius_atriceps,Pycnonotus_jocosus,Sylvia_atricapilla,Zosterops_lateralis,Phylloscopus_trochilus,Hirundo_rustica,Sturnus_vulgaris,Taeniopygia_guttata,Geospiza_fortis,Corvus_brachyrhynchos,Manacus_manacus,Smithornis_capensis,Acanthisitta_chloris,Melopsittacus_undulatus,Halcyon_senegalensis,Aquila_chrysaetos,Phalacrocorax_harrisi,Nipponia_nippon,Aptenodytes_forsteri,Cuculus_canorus,Caloenas_nicobarica,Gallus_gallus,Eudromia_elegans,Dromaius_novaehollandiae,Cairina_moschata 363-avian-2020.hal hal25_bulbul_all.maf

#convert to ss
msa_view --in-format MAF --out-format SS --unordered-ss hal25_bulbul_all.maf > hal25_bulbul_all.ss

#get .mod file with rates for conserved sequences
phastCons --expected-length 45 --target-coverage 0.3 --rho 0.4 --estimate-rho bulbul_conserved --no-post-probs --msa-format SS ../hal25_bulbul_all.ss bulbul_nonconserved-4d.mod

#calculate phastcons scores for conserved element set
phastCons --expected-length 45 --target-coverage 0.3 --most-conserved ../cnee_bulbul_sized_maf-fixed.bed --msa-format MAF ../hal25_bulbul_all.maf bulbul_conserved.mod,bulbul_nonconserved-4d.mod > bulbul_scores.wig

#get alignments of cnee sequences in maf format for 25 genomes from hal file using bedfile target split by chromosome
hal2mafMP.py  --numProc 6 --refGenome Brachypodius_atriceps --splitBySequence --noDupes --noAncestors --refTargets ../../cnee_bulbul_sized_maf-fixed.bed --targetGenomes Brachypodius_atriceps,Pycnonotus_jocosus,Sylvia_atricapilla,Zosterops_lateralis,Phylloscopus_trochilus,Hirundo_rustica,Sturnus_vulgaris,Taeniopygia_guttata,Geospiza_fortis,Corvus_brachyrhynchos,Manacus_manacus,Smithornis_capensis,Acanthisitta_chloris,Melopsittacus_undulatus,Halcyon_senegalensis,Aquila_chrysaetos,Phalacrocorax_harrisi,Nipponia_nippon,Aptenodytes_forsteri,Cuculus_canorus,Caloenas_nicobarica,Gallus_gallus,Eudromia_elegans,Dromaius_novaehollandiae,Cairina_moschata ../../363-avian-2020.hal hal25_bulbul_cnee.maf

#convert maf to fasta using custom script
maf2fasta_array.sbatch

#concatenate fastas into one file and generate a bed file with position using custom script
python concatinaTOR.py --fasta_folder=fasta/ --write_interval=10000 --output_folder=concat_fasta/ --n_taxa=25 --per_taxa=0.75
#this program writes one file for each taxa that needs to be concatenated together to create a single fasta file
cat *.fa > hal25_bulbul_all.fa

#run phyloacc
phyloacc.py -a hal25_bulbul_all.fa -b charset_0.txt -i charset_cnee_loci.txt -m bulbul_nonconserved-4d_fixed.mod -o phyloacc -t "Brachypodius_atriceps;Pycnonotus_jocosus" -g "Sylvia_atricapilla;Zosterops_lateralis;Phylloscopus_trochilus;Hirundo_rustica;Sturnus_vulgaris;Taeniopygia_guttata;Geospiza_fortis;Corvus_brachyrhynchos;Manacus_manacus;Smithornis_capensis;Acanthisitta_chloris;Melopsittacus_undulatus;Halcyon_senegalensis;Aquila_chrysaetos;Phalacrocorax_harrisi;Nipponia_nippon;Aptenodytes_forsteri;Cuculus_canorus;Caloenas_nicobarica;Gallus_gallus;Eudromia_elegans;Dromaius_novaehollandiae;Cairina_moschata" -n 4 -batch 100 -j 2 -part "shared"

snakemake -p -s /n/holyscratch01/informatics/sshakya/HAL/4d/phyloacc/phyloacc-job-files/snakemake/run_phyloacc.smk --configfile /n/holyscratch01/informatics/sshakya/HAL/4d/phyloacc/phyloacc-job-files/snakemake/phyloacc-config.yaml --profile /n/holyscratch01/informatics/sshakya/HAL/4d/phyloacc/phyloacc-job-files/snakemake/profiles/slurm_profile --dryrun

#snakemake --delete-all-output -p -s /n/holyscratch01/informatics/sshakya/HAL/4d/phyloacc/phyloacc-job-files/snakemake/run_phyloacc.smk --configfile /n/holyscratch01/informatics/sshakya/HAL/4d/phyloacc/phyloacc-job-files/snakemake/phyloacc-config.yaml --profile /n/holyscratch01/informatics/sshakya/HAL/4d/phyloacc/phyloacc-job-files/snakemake/profiles/slurm_profile --dryrun

#cactus
singularity exec --cleanenv ~/cactus_v2.1.1-gpu.sif cactus-prepare /n/holyscratch01/informatics/sshakya/genomes/sylvia_cactus_seq.txt --outDir /n/holyscratch01/informatics/sshakya/genomes/cactus_snakemake/snakemake_out --jobStore /n/holyscratch01/informatics/sshakya/genomes/cactus_snakemake/snakemake_temp --gpu

working_dir: /n/holyscratch01/informatics/sshakya/genomes/cactus_snakemake
cactus_path: ~/cactus_v2.1.1-gpu.sif
input_file: /n/holyscratch01/informatics/sshakya/genomes/cactus_snakemake/snakemake_out/sylvia_cactus_seq.txt
output_dir: /n/holyscratch01/informatics/sshakya/genomes/cactus_snakemake/snakemake_out/
tmp_dir: /n/holyscratch01/informatics/sshakya/genomes/cactus_snakemake/snakemake_temp/

snakemake -p -s ~/GenomeAnnotation/genome_alignment/cactus-gpu-snakemake/cactus_gpu.smk --configfile /n/holyscratch01/informatics/sshakya/genomes/cactus_snakemake/sylvia_config.yaml --profile ~/GenomeAnnotation/genome_alignment/cactus-gpu-snakemake/profiles/slurm_profile/ --dryrun
