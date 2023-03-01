import argparse
import pandas as pd
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), 'lib'))
import MafIO_edited

class mafidx:
    def __init__(self, maf, ref_sp):
        self.index = MafIO_edited.MafIndexEdit(maf + "index", maf, ref_sp)
        
parser = argparse.ArgumentParser()
parser.add_argument("--maf",
                    help="maf file or a group of maf files separated by a comma")
parser.add_argument("--bed",
                    help="bed for given scaffold")
parser.add_argument("--out_folder",
                    help="output_folder",
                    default="./")
parser.add_argument("--ref_species",
                    help="reference species")
parser.add_argument("--out_species",
                    help="text file for species to keep, one per line",
                    default=[])                    
parser.add_argument("--outmaf", action='store_true', help="Output a combined maf file instead of fasta")
parser.add_argument("--outmafname",
                    help="output name for maf file",
                    default="combined_maf") 
args = parser.parse_args()

if len(args.out_species) > 0:
    with open(args.out_species, 'r') as species_file:
        all_species = [x.strip("\n") for x in species_file.readlines()]
    species_file.close()

print("# generating maf index...")
mafs = [x.strip() for x in args.maf.split(",")]
mafsindex = []
for maf_file in mafs:
    mafsindex.append(mafidx(maf_file, args.ref_species))

print("# reading bed file...")
bed_df = pd.read_csv(args.bed, sep='\t', header=None)
new_df = pd.DataFrame(columns=["Scaffold", "Start","End","Name"])
if len(bed_df.columns) >= 12:
    for index, row in bed_df.iterrows():
        blockslen = row[10].split(",")
        blocksstart = row[11].split(",")
        blocksstop = [(int(blockslen[x]) + int(blocksstart[x])) for x in range(0,len(blockslen))]
        data = {'Scaffold': row[0],
                'Start': [(int(row[1]) + int(blocksstart[x])) for x in range(0,len(blockslen)-1)],
                'End': [(int(row[1]) + int(blocksstop[x])) for x in range(0,len(blockslen)-1)],
                'Name': [(row[3] + blocksstart[x]) for x in range(0,len(blockslen)-1)]}
        new_df = pd.concat([new_df, pd.DataFrame(data)])
    bed_df = new_df

print("# creating output folder: " + args.out_folder)
if not os.path.exists(args.out_folder):
   os.makedirs(args.out_folder)

if args.outmaf:
    print('# combining maf files..')
    fin_alignment = []
    for index, row in bed_df.iterrows():
        start, stop, ce_name = row[1], row[2], row[3] 
        aln_dict = {}
        mult = MultipleSeqAlignment([])
        for maf_index_file in mafsindex:
            multiple_alignment = maf_index_file.index.get_spliced([start], [stop], args.ref_species, strand=1)
            for alignment in multiple_alignment:
                species = alignment.id.split(".")[0]
                if len(args.out_species) > 0:
                    if (species not in aln_dict.keys()) & (species in all_species):
                        aln_dict[species] = alignment
                else:
                    if species not in aln_dict.keys():
                        aln_dict[species] =  alignment
        for key in aln_dict.keys():
            if key == args.ref_species:
                aln_dict[key].id = aln_dict[key].name
                mult.append(aln_dict[key])
                
        for key in aln_dict.keys():
            if key != args.ref_species:   
                aln_dict[key].id = aln_dict[key].name
                mult.append(aln_dict[key])   
        fin_alignment.append(mult)
    AlignIO.write(fin_alignment, os.path.join(args.out_folder, args.outmafname+".maf"), "maf")
    
else:
    print("# beginning fasta conversion...")
    for index, row in bed_df.iterrows():
        start, stop, ce_name = row[1], row[2], row[3]
        aln_dict = {}
        mult = MultipleSeqAlignment([])
        for maf_index_file in mafsindex:
            multiple_alignment = maf_index_file.index.get_spliced([start], [stop], args.ref_species, strand=1)
            for alignment in multiple_alignment:
                species = alignment.id
                if len(args.out_species) > 0:
                    if (species not in aln_dict.keys()) & (species in all_species):
                        aln_dict[species] = alignment
                else:
                    if species not in aln_dict.keys():
                        aln_dict[species] =  alignment
        for key in aln_dict.keys():
            if key == args.ref_species:
                mult.append(aln_dict[key])               
        for key in aln_dict.keys():
            if key != args.ref_species:   
                mult.append(aln_dict[key])   
        AlignIO.write(mult, os.path.join(args.out_folder, ce_name+".fa"), "fasta")
    
