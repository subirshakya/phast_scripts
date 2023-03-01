Scripts to get data ready for phastcons and phyloacc

runphast.sh has steps to go from hal to doing a phyloacc analysis

Once you get a folder of maf alignments per scaffold, use maf2fas.py to convert from maf to fasta alignments (fasta alignments will be based on input bed file and will be named with fourth column of bed file). It can take seven inputs:
| Options | Description |
| --------- | ----------- |
| --maf *filename(s)* | name of maf file/files; multiple files should be separated by commans |
| --bed *filename* | bed file containing elements to extract and convert to fasta; fourth column is name of output fasta; can parse out a 12 columned extended bed too |
| --out_folder *foldername* | output folder name; default is working directory |
| --ref_species *species* | the name of the reference taxa used for the maf alignment |
| --out_species *species.txt* | (optional) provide a text file with species to keep on each line |
| --outmaf | (optional) use flag if you want to output a combined maf file instead of a fasta file |
| --outmafname *name* | (optional) use to change name of output name of combined maf file; default is combined_maf |

The program maf2fas.py uses the AlignIO module of the Bio package. However, the MafIO module of AlignIO assigns species_name.scaffold_name as a separate entity if they come from two different scaffolds. The edited MafIO_edited.py code in the lib folder indexes the maf files based only on the species name and then produces only one sequence per taxa (concatenating anything that might be spread across two or more scaffolds). It also allows recovery of other metadata from a maf file.

concatinaTOR.py concatinates the fasta alignments into one single alignment (replaces different sized and missing blocks with Ns). It takes the following inputs:
| Options | Description |
| --------- | ----------- |
| --fasta_folder *foldername* | folder with fastas (fastas can in folders within this folder) |
| --bed_file *filename* | bed file with names of elements |
| --write-interval *number* | how frequently to write to fasta block (to save memory); default is 5 |
| --output_folder *foldername* | name of output folder |
| --n_taxa *number* | number of total taxa to expect from your alignments |
| --per_taxa *float* | filter any alignments that contain less than this fraction of taxa (e.g. 0.1 would omit any alignmets that contain less than 0.1 * n_taxa taxa) |
| --recursive | if you want to read fasta files within multiple folders within the fasta_folder |

insight.R is used to get data into format for running INSIGHT (http://compgen.cshl.edu/INSIGHT/) to use for looking at selection
