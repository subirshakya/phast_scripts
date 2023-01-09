Scripts to get data ready for phastcons and phyloacc

runphast.sh has steps to go from hal to doing a phyloacc analysis

Once you get a folder of maf alignments per scaffold, use maf2fas.py to convert from maf to fasta alignments (fasta alignments will be based on input bed file and will be named with fourth column of bed file). It takes five inputs:\
--maf > folder with maf files\
--bed > bed containing elements to extract and convert to fasta; fourth column is name of output fasta\
--out_folder > output folder name\
--ref_species > the name of the reference taxa used for the maf alignment\
--chromosome > the chromosome or scaffold name

The program maf2fas.py uses the AlignIO module of the Bio package. However, the MafIO module of AlignIO assigns species_name.scaffold_name as a separate entity if they come from two different scaffolds. The edited MafIO_edited.py code indexes the maf files based only on the species name and then produces only one sequence per taxa (concatenating anything that might be spread across two or more scaffolds). You can replace the original MafIO.py code in your anaconda > libs > Bio > AlignIO folder to make this work (I would rename the old MafIO.py to something else to recover that file if needed)

concatinaTOR.py concatinates the fasta alignments into one single alignment (replaces different sized and missing blocks with Ns). It takes the following inputs:\
--fasta_folder > folder with fastas (fastas can in folders within this folder)\
--bed_file > bed file with names of elements\
--write-interval > how frequently to write to fasta block (to save memory)\
--output_folder > name of output folder\
--n_taxa > number of total taxa to expect from your alignments\
--per_taxa > filter any alignments that contain less than this fraction of taxa (e.g. 0.1 would omit any alignmets that contain less than 0.1 * n_taxa taxa)\
--recursive > if you want to read fasta files within multiple folders within the fasta_folder

insight.R is used to get data into format for running INSIGHT (http://compgen.cshl.edu/INSIGHT/) to use for looking at selection
