import argparse
import pandas as pd
from Bio import AlignIO
import os

parser = argparse.ArgumentParser()
parser.add_argument("--maf",
                    help="maf file")
parser.add_argument("--bed",
                    help="bed for given scaffold")
parser.add_argument("--out_folder",
                    help="output_folder")
parser.add_argument("--ref_species",
                    help="reference species")
parser.add_argument("--chromosome",
                    help="chromosome/scaffold name")
args = parser.parse_args()

idx = AlignIO.MafIO.MafIndex(args.maf + "index", args.maf, args.ref_species)

bed_df = pd.read_csv(args.bed, sep='\t', header=None)

if not os.path.exists(args.out_folder):
   os.makedirs(args.out_folder)

for index, row in bed_df.iterrows():
    start, stop, ce_name = row[1], row[2], row[3]
    #print (ce_name)
    multiple_alignment = idx.get_spliced([start], [stop], strand=1)
    for alignment in multiple_alignment:
        alignment.id = alignment.id.split(".")[0]
    AlignIO.write(multiple_alignment, os.path.join(args.out_folder, ce_name+".fa"), "fasta")
