import argparse
import glob
from Bio import SeqIO
from Bio.Seq import Seq
import os

def list2bed(input_list, output_folder):
    with open(os.path.join(args.output_folder,"charset.txt"), 'a') as handle:
        for each_bed in input_list:
            handle.writelines(str(each_bed[0])+"\t"+str(each_bed[1])+"\t"+str(each_bed[2])+"\t"+str(each_bed[3])+"\n")
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
parser.add_argument("--write_interval", type=int, default=5,
                    help="how frequently (files) to write")
parser.add_argument("--output_folder", required=True,                   
                    help="output_folder")
parser.add_argument("--omit_missing", default=True)
parser.add_argument("--n_taxa", type=int, default=1)
parser.add_argument("--per_taxa", type=float, default=0.5)                    
args = parser.parse_args()

start = 1
stop = 0
all_files = glob.glob(args.fasta_folder+"*")
counter = 1
start_dict = {}
bed_index = []
if not os.path.exists(args.output_folder):
    os.makedirs(args.output_folder)
for each_file in all_files: 
    if counter == args.write_interval:
        list2bed(bed_index, args.output_folder)
        dict2fa(start_dict, args.output_folder)
        bed_index = []
        for key in start_dict:
            start_dict[key].seq = Seq("") 
        counter = 1
    else:
        my_dict = SeqIO.to_dict(SeqIO.parse(each_file, "fasta"))
        filename = os.path.splitext(os.path.basename(each_file))[0].split("_")
        id = filename[0]
        scaffold = filename[1] + "_" + filename[2]
        if args.omit_missing == True:
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
                bed_index.append((scaffold, start, stop, id))
                start = stop+1
                counter += 1
            else:
                print ("Omitting "+scaffold+":too few taxa in alignment")   
list2bed(bed_index, args.output_folder)
dict2fa(start_dict, args.output_folder, end=True)
