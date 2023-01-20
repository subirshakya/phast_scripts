import argparse
import pandas as pd
from Bio import AlignIO
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), 'lib'))
import MafIO_edited

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

print("# generating maf index...")
idx = MafIO_edited.MafIndexEdit(args.maf + "index", args.maf, args.ref_species)

print("# reading bed file...")
bed_df = pd.read_csv(args.bed, sep='\t', header=None)

print("# creating output folder: " + args.out_folder)
if not os.path.exists(args.out_folder):
   os.makedirs(args.out_folder)

print("# beginning fasta conversion...")
for index, row in bed_df.iterrows():
    start, stop, ce_name = row[1], row[2], row[3]
    #print (ce_name)
    multiple_alignment = idx.get_spliced([start], [stop], strand=1)
    for alignment in multiple_alignment:
        alignment.id = alignment.id.split(".")[0]
    AlignIO.write(multiple_alignment, os.path.join(args.out_folder, ce_name+".fa"), "fasta")
