Scripts to get data ready for phastcons and phyloacc

runphast.sh has steps to go from hal to doing a phyloacc analysis

Once you get a folder of maf alignments per scaffold, use maf2fas.sh (which uses organize_fasta.py) to convert from maf to fasta alignments (fasta alignments will be based on input bed file and will be named with fourth column of bed file)

concatinaTOR.py concatinates the fasta alignments into one single alignment (replaces different sized and missing blocks with Ns)

insight.R is used to get data into format for running INSIGHT (http://compgen.cshl.edu/INSIGHT/) to use for looking at selection
