import argparse
import pandas as pd
import glob
from Bio import SeqIO
from Bio.Seq import Seq
import os

parser = argparse.ArgumentParser()
parser.add_argument("--bed", 
                    help="bed for given scaffold")
parser.add_argument("--fasta_header",
                    help="header for fasta")
parser.add_argument("--out_folder",
                    help="output_folder")
parser.add_argument("--temp_folder",
                    help="location of temp folder")
args = parser.parse_args()

bed_df = pd.read_csv(args.bed, sep='\t', header=None)
for index, row in bed_df.iterrows():
    start, stop, ce_name = row[1], row[2], row[3]
    curr_pos = start
    start_dict = []
    while curr_pos < stop:
        try:
            my_dict = SeqIO.to_dict(SeqIO.parse(glob.glob(args.temp_folder+"/"+args.fasta_header+"-"+str(curr_pos)+"-"+"*.fa")[0], "fasta"))
            if len(start_dict) > 0:
                set_my = set(my_dict.keys())
                set_start = set(start_dict.keys())
                for taxa in (set_my & set_start):
                    seq_size = len(my_dict[taxa].seq)
                    base_size = len(start_dict[taxa].seq)
                    start_dict[taxa].seq = Seq(str(start_dict[taxa].seq) + str(my_dict[taxa].seq))
                for taxa in (set_my - set_start):
                    start_dict[taxa] = my_dict[taxa]
                    start_dict[taxa].seq = Seq(("N"*base_size) + str(my_dict[taxa].seq))
                for taxa in (set_start - set_my):
                    start_dict[taxa].seq = Seq(str(start_dict[taxa].seq) + "N"*seq_size)
            else:
                start_dict = my_dict            
            curr_pos += len(my_dict["Brachypodius_atriceps"].seq)
        except:
            curr_pos += 1    
            for taxa in start_dict:
                start_dict[taxa].seq = Seq(str(start_dict[taxa].seq) + "N")
    try:            
        with open(os.path.join(args.out_folder,ce_name+"_"+args.fasta_header+"_"+str(start)+"-"+str(stop)+".fa"), 'w') as handle:
            SeqIO.write(start_dict.values(), handle, 'fasta')
    except:
        print("Failed:", ce_name+"_"+args.fasta_header+"_"+str(start)+"-"+str(stop))
