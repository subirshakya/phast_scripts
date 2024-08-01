import argparse
import glob
from Bio import SeqIO
from Bio.Seq import Seq
import os
import pandas as pd

def list2bed(input_list, output_folder):
    with open(os.path.join(args.output_folder,"charset.txt"), 'a') as handle:
        for each_bed in input_list:
            handle.writelines(str(each_bed[0])+"\t"+str(each_bed[1]-1)+"\t"+str(each_bed[2])+"\t"+str(each_bed[3])+"\n")
    handle.close()
    
def dict2fa(input_dict, output_folder, end=False):
    for taxa in input_dict:
        if not os.path.exists(os.path.join(args.output_folder,start_dict[taxa].name+".fa")):
            with open(os.path.join(args.output_folder,start_dict[taxa].name+".fa"), 'a') as handle:
                handle.writelines(">"+start_dict[taxa].name+"\n")
                handle.writelines(str(start_dict[taxa].seq))
                if end == True:
                    handle.writelines("\n")
            handle.close()    
        else:
            with open(os.path.join(args.output_folder,start_dict[taxa].name+".fa"), 'a') as handle:
                handle.writelines(str(start_dict[taxa].seq))
                if end == True:
                    handle.writelines("\n")
            handle.close() 
            

parser = argparse.ArgumentParser()
parser.add_argument("--fasta_folder", required=True,
                    help="header for fasta")
parser.add_argument("--bed_file", help="bed file with element details")
parser.add_argument("--write_interval", type=int, default=5,
                    help="how frequently (files) to write")
parser.add_argument("--output_folder", required=True,                   
                    help="output_folder")
parser.add_argument("--n_taxa", type=int, default=1)
parser.add_argument("--per_taxa", type=float, default=0.5)
parser.add_argument("--recursive", action='store_true', default=False)              
args = parser.parse_args()

start = 1
stop = 0
if args.recursive==True:
    all_files = glob.glob(args.fasta_folder+"/**/*.fa", recursive=True)
else:
    all_files = glob.glob(args.fasta_folder+"*.fa")
counter = 1

bed_file = pd.read_csv(args.bed_file, sep='\t', header=None)

#keep only files present in bed file
files_base = [os.path.basename(x).split(".")[0] for x in all_files]
all_files_base = set(bed_file[3])
keep_files = [all_files[i] for i,x in enumerate(files_base) if x in all_files_base]

#initialize
start_dict = {}
bed_index = []
 
if not os.path.exists(args.output_folder):
    os.makedirs(args.output_folder)
else:
    print("Folder exists, stuff will be appended if files exist")
    
for each_file in keep_files:
    if counter == args.write_interval:
        list2bed(bed_index, args.output_folder)
        dict2fa(start_dict, args.output_folder)
        bed_index = []
        for key in start_dict:
            start_dict[key].seq = Seq("") 
        counter = 1
    my_dict = SeqIO.to_dict(SeqIO.parse(each_file, "fasta"))
    id = os.path.splitext(os.path.basename(each_file))[0]
    if len(my_dict) >= (args.per_taxa*args.n_taxa):          
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
            pass
        stop += len(my_dict[list(my_dict.keys())[0]].seq)
        bed_row = bed_file[bed_file[3].str.contains(id)]
        bed_index.append((bed_row.iloc[0][0], start, stop, id))
        start = stop+1
        counter += 1
    else:
        print ("Omitting "+id+":too few taxa in alignment; "+str(len(my_dict))+" of "+str(args.n_taxa))  
            
list2bed(bed_index, args.output_folder)
dict2fa(start_dict, args.output_folder, end=True)
