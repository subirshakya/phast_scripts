#phastcons
#get genepred bed
gff3ToGenePred -useName galgal7b.gff galgal7b_cds.gp
genePredToBed galgal7b_cds.gp galgal7b_cds.bed
split -l 10000 --numeric-suffixes=1 galgal7b_cds.bed ./4d/galgal7b_

#/n/holyscratch01/informatics/sshakya/genomes_phyloacc
~/hal/bin/hal4dExtract --conserved all_birds_new.hal Gallus_gallus7 ./4d/galgal7b_0${SLURM_ARRAY_TASK_ID} ./4d/galgal7b_4d_0${SLURM_ARRAY_TASK_ID}.bed

#merge all beds 
cat galgal7b_4d* > galgal7b_4d.bed #67664

#split across macro, micro, Z and W scaffolds
grep -f ../macro_chrom.txt galgal7b_4d.bed > galgal7b_4d_macro.bed #44426
grep -f ../micro_chrom.txt galgal7b_4d.bed > galgal7b_4d_micro.bed #19897
grep -f ../Z_chrom.txt galgal7b_4d.bed > galgal7b_4d_Z.bed #3104
grep -f ../W_chrom.txt galgal7b_4d.bed > galgal7b_4d_W.bed #237

#get 4D sites in maf format using bedfile
hal2mafMP.py  --numProc 6 --refGenome Gallus_gallus7 --noDupes --noAncestors --refTargets 4d/galgal7b_4d_macro.bed all_birds_new.hal 4d/all_birds_4d_macro.maf
hal2mafMP.py  --numProc 6 --refGenome Gallus_gallus7 --noDupes --noAncestors --refTargets 4d/galgal7b_4d_micro.bed all_birds_new.hal 4d/all_birds_4d_micro.maf
hal2mafMP.py  --numProc 6 --refGenome Gallus_gallus7 --noDupes --noAncestors --refTargets 4d/galgal7b_4d_Z.bed all_birds_new.hal 4d/all_birds_4d_Z.maf
hal2mafMP.py  --numProc 6 --refGenome Gallus_gallus7 --noDupes --noAncestors --refTargets 4d/galgal7b_4d_W.bed all_birds_new.hal 4d/all_birds_4d_W.maf

#convert to ss; ss did not work as expected so avoided
#msa_view --in-format MAF --out-format SS --unordered-ss 4d/all_birds_4d_macro.maf > 4d/all_birds_4d_macro.ss
#msa_view --in-format MAF --out-format SS --unordered-ss 4d/all_birds_4d_micro.maf > 4d/all_birds_4d_micro.ss
#msa_view --in-format MAF --out-format SS --unordered-ss 4d/all_birds_4d_Z.maf > 4d/all_birds_4d_Z.ss
#msa_view --in-format MAF --out-format SS --unordered-ss 4d/all_birds_4d_W.maf > 4d/all_birds_4d_W.ss

#run phylofit to get conserved rate
#tree:
#hal ((((((((((((((((((Certhia_americana,(Ficedula_albicollis,Sturnus_vulgaris)),(Taeniopygia_guttata,((((Diglossa_brunneiventris,(Emberiza_elegans,(Setophaga_coronata,Molothrus_ater))),Serinus_canaria),Motacilla_alba),Passer_domesticus))),((((((Sylvia_atricapilla,(Pycnonotus_jocosus,Brachypodius_melanocephalos)),Phylloscopus_trochilus),((Progne_subis,Riparia_riparia),Hirundo_rustica)),Acrocephalus_scirpaceus),Eremophila_alpestris),Parus_major)),(Corvus_cornix,Lanius_collurio)),(Lichenostomus_cassidix,Malurus_cyaneus)),(Myiozetetes_cayanensis,Chiroxiphia_lanceolata)),Acanthisitta_chloris),Melopsittacus_undulatus),Falco_peregrinus),(Athene_cunicularia,((Bucorvus_abyssinicus,((Colaptes_auratus,Pogoniulus_pusillus),((Chloroceryle_aenea,(Halcyon_senegalensis,Actenoides_hombroni)),Merops_nubicus))),Trogon_surrucura))),Aquila_chrysaetos),(((Sterna_hirundo,Pluvialis_apricaria),Phoenicopterus_ruber),(Theristicus_caerulescens,(Macronectes_giganteus,(((((Eudyptes_moseleyi,Eudyptes_sclateri),((Eudyptes_chrysolophus,Eudyptes_schlegeli),((Eudyptes_chrysocome,Eudyptes_filholi),Eudyptes_pachyrhynchus))),Megadyptes_antipodes),(((Spheniscus_demersus,Spheniscus_magellanicus),Spheniscus_humboldti),((Eudyptula_minor,Eudyptula_albosignata),Eudyptula_novaehollandiae))),(((Pygoscelis_antarcticus,Pygoscelis_papua),Pygoscelis_adeliae),(Aptenodytes_patagonicus,Aptenodytes_forsteri))))))),Porphyrio_hochstetteri),((Cuculus_canorus,Tauraco_erythrolophus),Columba_livia)),((Calypte_anna,Apus_apus),Caprimulgus_europaeus)),((Gallus_gallus4,Gallus_gallus7),Anas_platyrhynchos)),(((Dromaius_novaehollandiae,Eudromia_elegans),Rhea_americana),Struthio_camelus)),Gavialis_gangeticus);

phyloFit --tree "((((((((((((((((((Certhia_americana,(Ficedula_albicollis,Sturnus_vulgaris)),(Taeniopygia_guttata,((((Diglossa_brunneiventris,(Emberiza_elegans,(Setophaga_coronata,Molothrus_ater))),Serinus_canaria),Motacilla_alba),Passer_domesticus))),((((((Sylvia_atricapilla,(Pycnonotus_jocosus,Brachypodius_melanocephalos)),Phylloscopus_trochilus),((Progne_subis,Riparia_riparia),Hirundo_rustica)),Acrocephalus_scirpaceus),Eremophila_alpestris),Parus_major)),(Corvus_cornix,Lanius_collurio)),(Lichenostomus_cassidix,Malurus_cyaneus)),(Myiozetetes_cayanensis,Chiroxiphia_lanceolata)),Acanthisitta_chloris),Melopsittacus_undulatus),Falco_peregrinus),(Athene_cunicularia,((Bucorvus_abyssinicus,((Colaptes_auratus,Pogoniulus_pusillus),((Chloroceryle_aenea,(Halcyon_senegalensis,Actenoides_hombroni)),Merops_nubicus))),Trogon_surrucura))),Aquila_chrysaetos),(((Sterna_hirundo,Pluvialis_apricaria),Phoenicopterus_ruber),(Theristicus_caerulescens,(Macronectes_giganteus,(((((Eudyptes_moseleyi,Eudyptes_sclateri),((Eudyptes_chrysolophus,Eudyptes_schlegeli),((Eudyptes_chrysocome,Eudyptes_filholi),Eudyptes_pachyrhynchus))),Megadyptes_antipodes),(((Spheniscus_demersus,Spheniscus_magellanicus),Spheniscus_humboldti),((Eudyptula_minor,Eudyptula_albosignata),Eudyptula_novaehollandiae))),(((Pygoscelis_antarcticus,Pygoscelis_papua),Pygoscelis_adeliae),(Aptenodytes_patagonicus,Aptenodytes_forsteri))))))),Porphyrio_hochstetteri),((Cuculus_canorus,Tauraco_erythrolophus),Columba_livia)),((Calypte_anna,Apus_apus),Caprimulgus_europaeus)),((Gallus_gallus4,Gallus_gallus7),Anas_platyrhynchos)),(((Dromaius_novaehollandiae,Eudromia_elegans),Rhea_americana),Struthio_camelus)),Gavialis_gangeticus)" --init-random --subst-mod SSREV --sym-freqs --log 4d/all_birds_4d_macro.log --msa-format MAF --out-root 4d/all_birds_4d_macro_non-conserved 4d/all_birds_4d_macro.maf
phyloFit --tree "((((((((((((((((((Certhia_americana,(Ficedula_albicollis,Sturnus_vulgaris)),(Taeniopygia_guttata,((((Diglossa_brunneiventris,(Emberiza_elegans,(Setophaga_coronata,Molothrus_ater))),Serinus_canaria),Motacilla_alba),Passer_domesticus))),((((((Sylvia_atricapilla,(Pycnonotus_jocosus,Brachypodius_melanocephalos)),Phylloscopus_trochilus),((Progne_subis,Riparia_riparia),Hirundo_rustica)),Acrocephalus_scirpaceus),Eremophila_alpestris),Parus_major)),(Corvus_cornix,Lanius_collurio)),(Lichenostomus_cassidix,Malurus_cyaneus)),(Myiozetetes_cayanensis,Chiroxiphia_lanceolata)),Acanthisitta_chloris),Melopsittacus_undulatus),Falco_peregrinus),(Athene_cunicularia,((Bucorvus_abyssinicus,((Colaptes_auratus,Pogoniulus_pusillus),((Chloroceryle_aenea,(Halcyon_senegalensis,Actenoides_hombroni)),Merops_nubicus))),Trogon_surrucura))),Aquila_chrysaetos),(((Sterna_hirundo,Pluvialis_apricaria),Phoenicopterus_ruber),(Theristicus_caerulescens,(Macronectes_giganteus,(((((Eudyptes_moseleyi,Eudyptes_sclateri),((Eudyptes_chrysolophus,Eudyptes_schlegeli),((Eudyptes_chrysocome,Eudyptes_filholi),Eudyptes_pachyrhynchus))),Megadyptes_antipodes),(((Spheniscus_demersus,Spheniscus_magellanicus),Spheniscus_humboldti),((Eudyptula_minor,Eudyptula_albosignata),Eudyptula_novaehollandiae))),(((Pygoscelis_antarcticus,Pygoscelis_papua),Pygoscelis_adeliae),(Aptenodytes_patagonicus,Aptenodytes_forsteri))))))),Porphyrio_hochstetteri),((Cuculus_canorus,Tauraco_erythrolophus),Columba_livia)),((Calypte_anna,Apus_apus),Caprimulgus_europaeus)),((Gallus_gallus4,Gallus_gallus7),Anas_platyrhynchos)),(((Dromaius_novaehollandiae,Eudromia_elegans),Rhea_americana),Struthio_camelus)),Gavialis_gangeticus)" --init-random --subst-mod SSREV --sym-freqs --log 4d/all_birds_4d_micro.log --msa-format MAF --out-root 4d/all_birds_4d_micro_non-conserved 4d/all_birds_4d_micro.maf
phyloFit --tree "((((((((((((((((((Certhia_americana,(Ficedula_albicollis,Sturnus_vulgaris)),(Taeniopygia_guttata,((((Diglossa_brunneiventris,(Emberiza_elegans,(Setophaga_coronata,Molothrus_ater))),Serinus_canaria),Motacilla_alba),Passer_domesticus))),((((((Sylvia_atricapilla,(Pycnonotus_jocosus,Brachypodius_melanocephalos)),Phylloscopus_trochilus),((Progne_subis,Riparia_riparia),Hirundo_rustica)),Acrocephalus_scirpaceus),Eremophila_alpestris),Parus_major)),(Corvus_cornix,Lanius_collurio)),(Lichenostomus_cassidix,Malurus_cyaneus)),(Myiozetetes_cayanensis,Chiroxiphia_lanceolata)),Acanthisitta_chloris),Melopsittacus_undulatus),Falco_peregrinus),(Athene_cunicularia,((Bucorvus_abyssinicus,((Colaptes_auratus,Pogoniulus_pusillus),((Chloroceryle_aenea,(Halcyon_senegalensis,Actenoides_hombroni)),Merops_nubicus))),Trogon_surrucura))),Aquila_chrysaetos),(((Sterna_hirundo,Pluvialis_apricaria),Phoenicopterus_ruber),(Theristicus_caerulescens,(Macronectes_giganteus,(((((Eudyptes_moseleyi,Eudyptes_sclateri),((Eudyptes_chrysolophus,Eudyptes_schlegeli),((Eudyptes_chrysocome,Eudyptes_filholi),Eudyptes_pachyrhynchus))),Megadyptes_antipodes),(((Spheniscus_demersus,Spheniscus_magellanicus),Spheniscus_humboldti),((Eudyptula_minor,Eudyptula_albosignata),Eudyptula_novaehollandiae))),(((Pygoscelis_antarcticus,Pygoscelis_papua),Pygoscelis_adeliae),(Aptenodytes_patagonicus,Aptenodytes_forsteri))))))),Porphyrio_hochstetteri),((Cuculus_canorus,Tauraco_erythrolophus),Columba_livia)),((Calypte_anna,Apus_apus),Caprimulgus_europaeus)),((Gallus_gallus4,Gallus_gallus7),Anas_platyrhynchos)),(((Dromaius_novaehollandiae,Eudromia_elegans),Rhea_americana),Struthio_camelus)),Gavialis_gangeticus)" --init-random --subst-mod SSREV --sym-freqs --log 4d/all_birds_4d_Z.log --msa-format MAF --out-root 4d/all_birds_4d_Z_non-conserved 4d/all_birds_4d_Z.maf
phyloFit --tree "((((((((((((((((((Certhia_americana,(Ficedula_albicollis,Sturnus_vulgaris)),(Taeniopygia_guttata,((((Diglossa_brunneiventris,(Emberiza_elegans,(Setophaga_coronata,Molothrus_ater))),Serinus_canaria),Motacilla_alba),Passer_domesticus))),((((((Sylvia_atricapilla,(Pycnonotus_jocosus,Brachypodius_melanocephalos)),Phylloscopus_trochilus),((Progne_subis,Riparia_riparia),Hirundo_rustica)),Acrocephalus_scirpaceus),Eremophila_alpestris),Parus_major)),(Corvus_cornix,Lanius_collurio)),(Lichenostomus_cassidix,Malurus_cyaneus)),(Myiozetetes_cayanensis,Chiroxiphia_lanceolata)),Acanthisitta_chloris),Melopsittacus_undulatus),Falco_peregrinus),(Athene_cunicularia,((Bucorvus_abyssinicus,((Colaptes_auratus,Pogoniulus_pusillus),((Chloroceryle_aenea,(Halcyon_senegalensis,Actenoides_hombroni)),Merops_nubicus))),Trogon_surrucura))),Aquila_chrysaetos),(((Sterna_hirundo,Pluvialis_apricaria),Phoenicopterus_ruber),(Theristicus_caerulescens,(Macronectes_giganteus,(((((Eudyptes_moseleyi,Eudyptes_sclateri),((Eudyptes_chrysolophus,Eudyptes_schlegeli),((Eudyptes_chrysocome,Eudyptes_filholi),Eudyptes_pachyrhynchus))),Megadyptes_antipodes),(((Spheniscus_demersus,Spheniscus_magellanicus),Spheniscus_humboldti),((Eudyptula_minor,Eudyptula_albosignata),Eudyptula_novaehollandiae))),(((Pygoscelis_antarcticus,Pygoscelis_papua),Pygoscelis_adeliae),(Aptenodytes_patagonicus,Aptenodytes_forsteri))))))),Porphyrio_hochstetteri),((Cuculus_canorus,Tauraco_erythrolophus),Columba_livia)),((Calypte_anna,Apus_apus),Caprimulgus_europaeus)),((Gallus_gallus4,Gallus_gallus7),Anas_platyrhynchos)),(((Dromaius_novaehollandiae,Eudromia_elegans),Rhea_americana),Struthio_camelus)),Gavialis_gangeticus)" --init-random --subst-mod SSREV --sym-freqs --log 4d/all_birds_4d_W.log --msa-format MAF --out-root 4d/all_birds_4d_W_non-conserved 4d/all_birds_4d_W.maf

#get cnee
grep -P "\tgene\t" galgal7b.gff > galgal7b_cds.gff
bedtools intersect -a cnee_gal7.bed -b galgal7b_cds.gff -v > cnee_gal7_nonexonic.bed

#get maf for elements of each scaffold, $1 takes file with scaffold name and $2 takes scaffold type (Z, W, macro, micro)
while IFS= read -r line; do
  mkdir -p cnee/maf_$2
  grep -P $line cnee_gal7_nonexonic.bed > cnee/${line}.bed
  hal2mafMP.py --numProc 6 --refGenome Gallus_gallus7 --noDupes --noAncestors --refTargets cnee/${line}.bed all_birds_new.hal cnee/maf_$2/all_birds_${line}.maf
  rm cnee/${line}.bed
done < "$1"

#get fasta for each maf; note: i broke this up further to parallelize this
./maf2fasta.sh -i cnee/maf_W/ -o cnee/fasta_W/ -n all_birds -b cnee/cnee_gal7_nonexonic_annot_W.bed -s Gallus_gallus7
./maf2fasta.sh -i cnee/maf_Z/ -o cnee/fasta_Z/ -n all_birds -b cnee/cnee_gal7_nonexonic_annot_Z.bed -s Gallus_gallus7
./maf2fasta.sh -i cnee/maf_micro/ -o cnee/fasta_micro/ -n all_birds -b cnee/cnee_gal7_nonexonic_annot_micro.bed -s Gallus_gallus7
./maf2fasta.sh -i cnee/maf_macro/ -o cnee/fasta_macro/ -n all_birds -b cnee/cnee_gal7_nonexonic_annot_macro.bed -s Gallus_gallus7

#get mafs for all species except targets
hal2mafMP.py  --numProc 6 --refGenome Gallus_gallus7 --noDupes --noAncestors --splitBySequence --targetGenomes Gavialis_gangeticus,Struthio_camelus,Rhea_americana,Eudromia_elegans,Dromaius_novaehollandiae,Anas_platyrhynchos,Gallus_gallus7,Phoenicopterus_ruber,Columba_livia,Tauraco_erythrolophus,Cuculus_canorus,Caprimulgus_europaeus,Apus_apus,Calypte_anna,Porphyrio_hochstetteri,Pluvialis_apricaria,Sterna_hirundo,Macronectes_giganteus,Theristicus_caerulescens,Aquila_chrysaetos,Athene_cunicularia,Trogon_surrucura,Bucorvus_abyssinicus,Merops_nubicus,Pogoniulus_pusillus,Colaptes_auratus,Falco_peregrinus,Melopsittacus_undulatus,Acanthisitta_chloris,Myiozetetes_cayanensis,Chiroxiphia_lanceolata,Malurus_cyaneus,Lichenostomus_cassidix,Lanius_collurio,Corvus_cornix,Eremophila_alpestris,Parus_major,Acrocephalus_scirpaceus,Phylloscopus_trochilus,Sylvia_atricapilla,Certhia_americana,Sturnus_vulgaris,Ficedula_albicollis,Taeniopygia_guttata,Passer_domesticus,Motacilla_alba,Serinus_canaria,Emberiza_elegans,Molothrus_ater,Setophaga_coronata,Diglossa_brunneiventris all_birds_new.hal phastcons/non_targets_all_birds

#get mod files and remove targets from phylogenetic tree
#done using keep.tips in R::ape

#run phastcons first to estimate rho 
#For Z scaffold (GC=0.41)
phastCons --expected-length 45 --target-coverage 0.3 --rho 0.4 --estimate-rho non_targets_4d_NC_052572.1 --no-post-probs --msa-format MAF non_targets_all_birds_NC_052572.1.maf non_targets_4d_Z_non-conserved.mod
#For micro scaffolds (GC=0.52)
#For macro sacffolds (GC=0.42)
#micro run in serial, macro run in parallel

#use phastboots to combine runs for conserved rates for each group; usage: phyloBoot --read-mods '*macro_cons.txt' --output-average non_targets_4d_macro_conserved.mod

#caculate phastcons scores
phastCons --expected-length 45 --target-coverage 0.3 --most-conserved cons_beds/conserved_Z_NC_052572.1.bed --msa-format MAF non_targets_all_birds_NC_052572.1.maf non_targets_4d_Z_conserved.mod,non_targets_4d_Z_non-conserved.mod > cons_wig/conserved_Z_NC_052572.1.wig
#repeat for micro and macro scaffolds; each run in serial

#repeated phastcons with Tgut reference using same phastcons mod files as from previous

#take conserved_all_phastcon.bed (in galgal7) and conserved_tgut_galgal7_d5.bed (in tgut > halliftover to galgal7; merged with d 5); both filtered for 20 bp & 5000 bp
#merge these two phastcon datasets
cat conserved_all_phastcon_filt.bed conserved_tgut_galgal7_filt.bed | bedtools sort -i "stdin" | bedtools merge -d 5 -i "stdin" > conserved_phast2.bed

#get uniques in 3 datasets that is absent in conserved_phast2
bedtools subtract -A -a size_remapped_cnees.galGal4NCBI.bed -b conserved_phast2_filt.bed > ratite_unique.bed
bedtools subtract -A -a size_remapped_vocal.bed -b conserved_phast2_filt.bed > vocal_unique.bed
bedtools subtract -A -a size_remapped_Hum2chick_CNEE_NCBI.bed -b conserved_phast2_filt.bed > UCSC_unique.bed

#get conserved area under uniqueness for each of the 3 unique sets
bedtools intersect -wao -a ratite_unique.bed -b ancRep_separate_models_rev.bw.conserved_galgal7_merge.bed > ratite_unique_cov.bed
...

#use R to aggregate and get fractions
#combine conserved_phast2 with the other three unique datasets, merge (for any overlap among the uniques) and filter (+remove the W chromosome scaffolds)
cat conserved_phast2_filt.bed ratite_unique.bed vocal_unique.bed UCSC_unique.bed | bedtools sort -i "stdin" | bedtools merge -i "stdin" > final_working_conserved.bed

#run PhyloP to get conservation values of each element
#first need to get maf alignments of each scaffold with the respective elements
while IFS= read -r line; do
  mkdir -p final_dataset/maf_${2}_sort
  mkdir -p final_dataset/beds/
  grep -P $line final_dataset/final_working_conserved_filt.bed > final_dataset/beds/${line}.bed
  hal2mafMP.py --numProc 6 --refGenome Gallus_gallus7 --noDupes --noAncestors --refTargets final_dataset/beds/${line}.bed all_birds_new.hal final_dataset/maf_${2}_sort/all_birds_${line}.maf
done < "$1"

#use phylop to calculate phylop scores using respective conserved rates
phyloP --method LRT --mode CONACC -i MAF --features final_dataset/beds/NC_052572.1.bed 4d/all_birds_4d_Z_non-conserved.mod final_dataset/maf_Z_sort/all_birds_NC_052572.1.maf > final_dataset/all_birds_phylop_Z_NC_052572.1.out

#combine phylop values of each scaffold

#annotate file using bedtools and then R
bedtools annotate -counts -i final_working_conserved_filt.bed -files coding_exon_noutr.bed intron.bed conserved_all_phastcon_filt.bed conserved_tgut_galgal7_filt.bed ../size_remapped_Hum2chick_CNEE_NCBI_merged.bed ../size_remapped_cnees.galGal4NCBI.bed ../size_remapped_vocal.bed atac_all.bed > final_working_conserved_filt_annot.bed

#get fasta files for each element
python maf2fasta.py --maf maf_Z_sort/all_birds_NC_052572.1.maf --bed beds2/NC_052572.1.bed --out_folder fasta_Z/ --ref_species Gallus_gallus7 --chromosome NC_052572.1

#concatenate fasta files 
python ConcatinaTOR.py --fasta_folder fasta_Z --bed_file cnee_Z.txt --write_interval 1000 --output_folder fasta_concat/fasta_Z_concat/ --omit_missing --n_taxa 79 --per_taxa 0.10 --recursive

#run phyloacc
#update the mod file to include ancestral node names on tree, for some reason the default way in phyloacc does not work
#can use phyloacc.py -m all_birds_4d_Z_non-conserved.mod -n 1 --labeltree to get a nodelabelled tree but still need to update the labels a little bit
((((((((((((((((((Certhia_americana:0.0599534,(Ficedula_albicollis:0.0450892,Sturnus_vulgaris:0.0459971)Anc1:0.0149191)Anc22:0.00292677,(Taeniopygia_guttata:0.0489198,((((Diglossa_brunneiventris:0.0218787,(Emberiza_elegans:0.026383,(Setophaga_coronata:0.0214966,Molothrus_ater:0.0169871)Anc2:0.00645918)Anc23:0.00122959)Anc36:0.00955365,Serinus_canaria:0.0338919)Anc43:0.00682624,Motacilla_alba:0.0196812)Anc47:0.00485456,Passer_domesticus:0.04192)Anc51:0.00544704)Anc55:0.0112305)Anc59:0.028802,((((((Sylvia_atricapilla:0.0399216,(Pycnonotus_jocosus:0.0202878,Brachypodius_melanocephalos:0.0181029)Anc3:0.0310306)Anc24:0.0202979,Phylloscopus_trochilus:0.0324266)Anc37:0.0331933,((Progne_subis:0.00677215,Riparia_riparia:0.021334)Anc4:0.0227219,Hirundo_rustica:0.0182412)Anc25:0.0408118)Anc44:0.0117071,Acrocephalus_scirpaceus:0.0523735)Anc48:0.00971589,Eremophila_alpestris:0.0704565)Anc52:0.00435551,Parus_major:0.0364579)Anc56:0.00489059)Anc61:0.010584,(Corvus_cornix:0.023756,Lanius_collurio:0.0291614)Anc5:0.0165007)Anc63:0.00352419,(Lichenostomus_cassidix:0.0639511,Malurus_cyaneus:0.0721923)Anc6:0.021343)Anc65:0.0467941,(Myiozetetes_cayanensis:0.0368835,Chiroxiphia_lanceolata:0.0427223)Anc7:0.0236695)Anc66:0.0103418,Acanthisitta_chloris:0.098451)Anc67:0.0342691,Melopsittacus_undulatus:0.0763502)Anc68:0.00500232,Falco_peregrinus:0.0875945)Anc69:0.00394511,(Athene_cunicularia:0.0763883,((Bucorvus_abyssinicus:0.0968766,((Colaptes_auratus:0.0580542,Pogoniulus_pusillus:0.046023)Anc8:0.0784409,((Chloroceryle_aenea:0.0695274,(Halcyon_senegalensis:0.0259379,Actenoides_hombroni:0.027348)Anc9:0.0149656)Anc26:0.0377283,Merops_nubicus:0.121897)Anc38:0.00766596)Anc45:0.00457213)Anc49:0.00463949,Trogon_surrucura:0.125886)Anc53:0.000454279)Anc57:0.0158584)Anc70:0.00436846,Aquila_chrysaetos:0.0567426)Anc71:0.000635533,(((Sterna_hirundo:0.0490513,Pluvialis_apricaria:0.0362526)Anc10:0.00607084,Phoenicopterus_ruber:0.0337325)Anc27:0.0179159,(Theristicus_caerulescens:0.0334916,(Macronectes_giganteus:0.0300689,(((((Eudyptes_moseleyi:0.00603947,Eudyptes_sclateri:0.00595532)Anc11:0.0018661,((Eudyptes_chrysolophus:0.103578,Eudyptes_schlegeli:0.00355659)Anc12:0.0751101,((Eudyptes_chrysocome:0.00236426,Eudyptes_filholi:0.000490251)Anc13:0.00788461,Eudyptes_pachyrhynchus:0.0156512)Anc28:7.63159e-05)Anc39:0.00386822)Anc46:0.0027673,Megadyptes_antipodes:0.0122433)Anc50:0.0100572,(((Spheniscus_demersus:4.45702e-05,Spheniscus_magellanicus:0.0157264)Anc14:0.00114764,Spheniscus_humboldti:0.0037014)Anc29:0.00732771,((Eudyptula_minor:0.0148284,Eudyptula_albosignata:0.0133185)Anc15:0.00385652,Eudyptula_novaehollandiae:0.00225995)Anc30:0.00769011)Anc40:0.000297749)Anc54:0.00634312,(((Pygoscelis_antarcticus:0.0157893,Pygoscelis_papua:0.0178978)Anc16:0.00043673,Pygoscelis_adeliae:0.0079677)Anc31:0.0093239,(Aptenodytes_patagonicus:0.00795276,Aptenodytes_forsteri:0.00191802)Anc17:0.0115336)Anc41:0.0051001)Anc58:0.0232999)Anc60:0.00569994)Anc62:0.00690889)Anc64:0.0067635)Anc72:0.00795601,Porphyrio_hochstetteri:0.143894)Anc73:0.0015046,((Cuculus_canorus:0.138061,Tauraco_erythrolophus:0.0826187)Anc18:0.00422666,Columba_livia:0.115953)Anc32:0.0029288)Anc74:0.00532867,((Calypte_anna:0.129922,Apus_apus:0.0961625)Anc19:0.0136657,Caprimulgus_europaeus:0.0847003)Anc33:0.0039344)Anc75:0.0286412,((Gallus_gallus4:0.0144136,Gallus_gallus7:0.00747648)Anc20:0.15473,Anas_platyrhynchos:0.109813)Anc34:0.0170549)Anc76:0.0323037,(((Dromaius_novaehollandiae:0.0264611,Eudromia_elegans:0.217703)Anc21:0.00861333,Rhea_americana:0.0743868)Anc35:0.0142153,Struthio_camelus:0.0605634)Anc42:0.0269703)Anc77:0.109503,Gavialis_gangeticus:0.109503)root

phyloacc.py --overwrite -a cnee_Z.fa -b charset_fixed.txt -i charset_id.txt -m all_birds_4d_Z_non-conserved.mod -o output -t "Aptenodytes_patagonicus;Aptenodytes_forsteri;Pygoscelis_adeliae;Pygoscelis_papua;Pygoscelis_antarcticus;Megadyptes_antipodes;Eudyptula_novaehollandiae;Eudyptula_albosignata;Eudyptula_minor;Spheniscus_demersus;Spheniscus_humboldti;Spheniscus_magellanicus;Eudyptes_pachyrhynchus;Eudyptes_sclateri;Eudyptes_chrysolophus;Eudyptes_schlegeli;Eudyptes_chrysocome;Eudyptes_filholi;Eudyptes_moseleyi;Halcyon_senegalensis;Chloroceryle_aenea;Actenoides_hombroni;Progne_subis;Riparia_riparia;Hirundo_rustica;Brachypodius_melanocephalos;Pycnonotus_jocosus" -g "Gavialis_gangeticus" -n 4 -batch 100 -j 20 -part "shared"
