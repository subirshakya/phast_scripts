library(PopGenome)
library(plyr)
library(rphast)
library(stringr)
library(ape)
library(geiger)

###INPUTS
#bed file should have four columns Scaffold Name, Start, Stop, and ID
#fasta files in folder should be for alignments of elements including flanking region
#fasta files should be named based on ID in bed file
#focal taxa should be in phylo_tree
bed_file <- "C:/Subir/Projects/Patriceps/CNEE/insight/bulbul_accl.bed" 
vcf_file <- "C:/Subir/Projects/Patriceps/CNEE/patri12_bi_nomissing.recode.vcf.gz"
fasta_folder <- "C:/Subir/Projects/Patriceps/CNEE/insight/fasta_flank/"
phylo_tree <- "C:/Subir/Projects/Patriceps/CNEE/insight/bulbul_4d.tree"
focal_taxa <- "Brachypodius_atriceps"
sister_taxa <- "Pycnonotus_jocosus"
n_indv <- 12
###

###OUTPUTS
#end path in /
output_folder <- "C:/Subir/Projects/Patriceps/CNEE/insight/"
output_name <- "flank_HH18"
###

###OPTIONS
flanking_extent <- 150 #how far away from element to look for neutral flanking sites
flanking_buffer <- 50 #do you want any gap between element and flanking region
diploid <- TRUE
###

get_pop_stats <- function(vcf_file, val1, val2, val3){
  #returns dataframe with bed index, ref and alt snps, snp position, and frequency
  nucls <- c("T", "C", "G", "A", "N", "-") #based on how popgenome codes
  flank_start <- val2-150
  if (flank_start <= 0){flank_start = 1}
  flank_end <- val3+150
  read_vcf <- readVCF(vcf_file, 1000, val1, flank_start, flank_end, approx = F, out=paste("temp", "_test", sep=""))
  vcf_stats <- neutrality.stats(read_vcf, FAST=F)
  vcf_stats <- detail.stats(vcf_stats, biallelic.structure = T)
  vcf_stats <- diversity.stats(vcf_stats, pi = T)
  if (vcf_stats@n.biallelic.sites == 0){
    site_file <- NULL
  }
  else{
    site_file <- data.frame(chrom=val1, start=val2, flank_start=flank_start,
                            end=val3, flank_end=flank_end, theta=(vcf_stats@theta_Watterson/vcf_stats@n.sites)[1],
                            snp_id=vcf_stats@region.data@biallelic.sites[[1]],
                            snp_major=nucls[vcf_stats@region.data@biallelic.substitutions[[1]][2,]],
                            snp_minor=nucls[vcf_stats@region.data@biallelic.substitutions[[1]][1,]],
                            snp_freq=as.vector(vcf_stats@region.stats@minor.allele.freqs[[1]])*length(vcf_stats@region.data@populations2[[1]][[1]]))
    site_file$snp_index <- site_file$snp_id - site_file$flank_start
  }
  return (site_file)
}

#read inputs
accl_bed <- read.table(bed_file, sep = "\t", col.names = c("Scaffold", "Start", "Stop", "ID"))
accl_bed$flank_left_start <- accl_bed$Start - flanking_extent
accl_bed$flank_left_stop <- accl_bed$Start - flanking_buffer
accl_bed$flank_right_start <- accl_bed$Stop + flanking_buffer
accl_bed$flank_right_stop <- accl_bed$Stop + flanking_extent

all_fastas <- list.files(fasta_folder, full.names = T)

tree <- read.newick.tree(phylo_tree)
tree_ape <- read.tree(phylo_tree)
child_mrca <- tips(tree_ape, getMRCA(tree_ape, c(focal_taxa, sister_taxa)))
drop_child <- setdiff(child_mrca, c(focal_taxa, sister_taxa))

#open output files
fileConn <- paste(output_folder,output_name,"_conserved.txt", sep = "")
fileFlank <- paste(output_folder,output_name,"_flank.txt", sep = "")
if (diploid){
  n_indv <- 2*n_indv
}
write(paste("samples", n_indv, sep = " "), fileConn)
write(paste("samples", n_indv, sep = " "), fileFlank)

alignment <- all_fastas[9]
for (alignment in all_fastas[1:10]){
  #get full alignment
  align <- read.msa(alignment, format = "FASTA")
  if (length(drop_child) > 0){
    align <- align[-(which(align$name %in% drop_child)),]
  }
  #get ce name
  ce_align <- gsub(fasta_folder, "", alignment, fixed = T)
  ce_align <- gsub(".fa[sta]{0,3}", "", ce_align)
  index_main_bed <- which(accl_bed$ID %in% ce_align)

  #subset flank alignment
  flanking_length <- flanking_extent - flanking_buffer
  align_flank <- align[,c(1:flanking_length,(nchar(align[1,]$seq)-flanking_length+1):(nchar(align[1,]$seq)))]
  
  #calculate lambda
  phylo_vals <- phyloFit(align_flank, tree)
  accl_bed[index_main_bed, "lambda"] <- gsub(paste(focal_taxa,  ":", sep=""), "",
                                                                                  str_match(phylo_vals$tree, paste(focal_taxa,  ":[0-9]*.[0-9]*", sep=""))[1],
                                                                                  fixed = T)
  
  #mask focal sequence and calculate probability of bases
  target_seq <- align[focal_taxa,]$seq
  align[focal_taxa,] <- "N"
  postprob <- postprob.msa(x = align, tm = phylo_vals, every.site = T)
  target_branch <- rownames(postprob)[which(rownames(postprob) %in% c(paste(focal_taxa,sister_taxa,sep="-"), 
                                                                      paste(sister_taxa,focal_taxa,sep="-")))]
  
  #read vcf file for region
  site_file <- get_pop_stats(vcf_file,
                             accl_bed[index_main_bed, "Scaffold"], 
                             accl_bed[index_main_bed, "Start"],
                             accl_bed[index_main_bed, "Stop"])
  if (is.null(site_file) == FALSE){
    write(paste("block", 
                paste(site_file$chrom[1], ":", site_file$flank_start[1], "-", site_file$flank_end[1], sep=""),
                "theta", site_file$theta[1],
                "lambda", accl_bed[index_main_bed, "lambda"], sep = "\t"), fileConn, append = T)
    write(paste("block", 
                paste(site_file$chrom[1], ":", site_file$flank_start[1], "-", site_file$flank_end[1], sep=""),
                "theta", site_file$theta[1],
                "lambda", accl_bed[index_main_bed, "lambda"], sep = "\t"), fileFlank, append = T)
    
    #prob file is 3D, first D is node name, second is base name (A,C,T,G), third is nucleotide position
    #get prob of sites
    for (site in seq(1, nrow(site_file))){
      if (site_file[site,"snp_id"] < site_file$start[1] | site_file[site,"snp_id"] > site_file$end[1]){
        write(paste("site",
                    paste(site_file[site,"chrom"], ":", site_file[site,"snp_id"], sep = ""),
                    "P",
                    postprob[target_branch,
                            site_file[site,"snp_major"], site_file[site,"snp_index"]],
                    postprob[target_branch,
                             site_file[site,"snp_minor"], site_file[site,"snp_index"]],
                    n_indv - site_file[site,"snp_freq"], site_file[site,"snp_freq"], sep="\t"), fileFlank, append = T)
        }
    }
    for (pos in seq(site_file$start[1], (site_file$end[1]))){
      if (pos %in% site_file$snp_id){
        pos_index <- site_file[which(pos %in% site_file$snp_id),]
        write(paste("site",
                    paste(pos_index["chrom"], ":", pos, sep = ""),
                    "P",
                    postprob[target_branch,
                             pos_index$snp_major, pos_index$snp_index],
                    postprob[target_branch,
                             pos_index$snp_minor, pos_index$snp_index],
                    n_indv - pos_index$snp_freq, pos_index$snp_freq, sep="\t"), fileConn, append = T)
      }
      else{
        if (substr(target_seq, pos - site_file$flank_start[1] - 1, pos-site_file$flank_start[1] - 1) %in% c("N","-")){
          write(paste("site",
                      paste(site_file$chrom[1], ":", pos, sep = ""),
                      "M",
                      1,
                      sep="\t"), fileConn, append = T)
        }
        else{
          write(paste("site",
                      paste(site_file$chrom[1], ":", pos, sep = ""),
                      "M",
                      postprob[target_branch,
                               substr(target_seq, pos - site_file$flank_start[1] - 1, pos-site_file$flank_start[1] - 1), 
                               pos - site_file$flank_start[1] - 1],
                      sep="\t"), fileConn, append = T)
        }
      }
    }
  }
}
