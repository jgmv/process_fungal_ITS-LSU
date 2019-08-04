#!/usr/bin/env bash


### prepare environment
# create folders for scripts and data
mkdir -p scripts
mkdir -p fasta_files
cp $1 fasta_files/its.fasta
cp $2 fasta_files/lsu.fasta


### download necessary scripts
# custom bash scripts
if test -f ./scripts/bash_seq_analysis/otuList.sh
then
  echo "./scripts/bash_seq_analysis found"
else
  git clone https://github.com/jgmv/bash_seq_analysis.git scripts/bash_seq_analysis

fi

for i in $(ls scripts/bash_seq_analysis/*.sh)
do
  . $i
done


# ITSx, if not locally installed
if command -v ITSx >/dev/null 2>&1
then
  echo "ITSx found"
elif command -v ./scripts/ITSx/ITSx >/dev/null 2>&1
then
  echo "./scripts/ITSx/ITSx found"
else
  git clone https://github.com/ncbi/ITSx.git scripts/ITSx
fi

# BeforePhylo.pl
if command -v BeforePhylo >/dev/null 2>&1
then
  echo "BeforePhylo found"
elif test -f ./scripts/BeforePhylo/BeforePhylo.pl
then
  echo "./scripts/BeforePhylo/BeforePhylo.pl found"
else
  git clone https://github.com/qiyunzhu/BeforePhylo.git scripts/BeforePhylo
fi


### ITS sequence identification
# download last version of UNITE ITS reference dataset, mothur release
# check https://unite.ut.ee/repository.php
if test -f ./data/UNITE_sh_dynamic.tax
then
  echo "UNITE database found"
else
  mkdir -p data
  wget -P data https://files.plutof.ut.ee/public/orig/56/25/5625BDC830DC246F5B8C7004220089E032CC33EEF515C76CD0D92F25BDFA9F78.zip
  unzip data/*.zip
  rm data/*.zip
  mv data/UNITE*_sh_dynamic.fasta data/UNITE_sh_dynamic.fasta
  mv data/UNITE*_sh_dynamic.tax data/UNITE_sh_dynamic.tax
  rm data/UNITEv*
fi

# identify ITS sequences using mothur's NBC
mothur "#classify.seqs(fasta=fasta_files/its.fasta,\
  template=data/UNITE_sh_dynamic.fasta,\
  taxonomy=data/UNITE_sh_dynamic.tax, cutoff=60, probs=T)" 
mkdir -p its_identification
mv fasta_files/its.UNITE* its_identification
mv mothur.* its_identification
removeTaxonTag \
  its_identification/its.UNITE_sh_dynamic.wang.taxonomy \
  its_identification/taxonomy_boot.csv

# create a copy of the taxonomy file without bootstrap values
cp its_identification/taxonomy_boot.csv its_identification/taxonomy.csv
sed -i 's/([^()]*)//g' its_identification/taxonomy.csv


### split ITS regions and store in independent files
mkdir -p its_regions

if command -v ITSx >/dev/null 2>&1
then
  ITSx -i fasta_files/its.fasta -o its_split -t F --summary F \
    --graphical F --preserve T --not-found F --save_regions ITS1,5.8S,ITS2
else
  scripts/ITSx/ITSx -i fasta_files/its.fasta -o its_split -t F --summary F \
    --graphical F --preserve T --not-found F --save_regions ITS1,5.8S,ITS2 \
    --reset T
fi
mv its_split* its_regions
cp its_regions/its_split.ITS1.fasta fasta_files/its_ITS1.fasta
cp its_regions/its_split.ITS2.fasta fasta_files/its_ITS2.fasta
cp its_regions/its_split.5_8S.fasta fasta_files/its_5_8S.fasta


### group isolates into ITS-based OTUs
# cluster sequences into OTUs
mkdir -p otu_clustering
blastclust -i fasta_files/its.fasta -o  \
  otu_clustering/otu_clusters.csv -p F -b T -S 97
mv error.log otu_clustering
otuList otu_clustering/otu_clusters.csv otu_clustering/otu_list.csv
otuList2mothur otu_clustering/otu_list.csv otu_clustering/otu_list_mothur


### generate average taxonomies for OTUs
mothur "#classify.otu(taxonomy=its_identification/its.UNITE_sh_dynamic.wang.taxonomy,\
  list=otu_clustering/otu_list_mothur, probs=F)"
mv mothur.* otu_clustering
removeTaxonTag otu_clustering/otu_list_mothur0.03.cons.taxonomy \
  otu_clustering/otu_taxonomy.csv

# select representative OTU sequences
sed 's/\s.*$//' otu_clustering/otu_clusters.csv > otu_clustering/otu_reps.csv
for i in fasta_files/its.fasta fasta_files/lsu.fasta  \
   fasta_files/its_ITS1.fasta fasta_files/its_ITS2.fasta  \
   fasta_files/its_5_8S.fasta
do
  out=${i::-6}
  getRepSeqs otu_clustering/otu_reps.csv $i ${out}_rep.fasta
  annotateRepSeqs otu_clustering/otu_list.csv ${out}_rep.fasta \
    ${out}_rep_otus.fasta
done


### phylogeny with LSU sequences
# with all sequences
mkdir -p lsu_phylogeny
mafft --auto fasta_files/lsu.fasta > lsu_phylogeny/lsu_align.fasta
Gblocks lsu_phylogeny/lsu_align.fasta -t=DNA -b4=5 -b5=h
raxmlHPC-SSE3 -s lsu_phylogeny/lsu_align.fasta-gb -n LSU -f a -m GTRGAMMA  \
  -x 12345 -p 12345 -# 100
mv RAxML* lsu_phylogeny

# with OTU representative sequences
mafft --auto fasta_files/lsu_rep.fasta > lsu_phylogeny/lsu_rep_align.fasta
Gblocks lsu_phylogeny/lsu_rep_align.fasta -t=DNA -b4=5 -b5=h
raxmlHPC-SSE3 -s lsu_phylogeny/lsu_rep_align.fasta-gb -n LSU_rep  \
  -f a -m GTRGAMMA -x 12345 -p 12345 -# 100
mv RAxML* lsu_phylogeny


### multilocus phylogeny with representative sequences
# create alignments for all fasta files
mkdir -p multilocus_phylogeny
for i in $( ls fasta_files/*_rep.fasta )
do
  out=multilocus_phylogeny/${i:12:-6}_align.fasta
  mafft --auto $i > $out
  Gblocks $out -t=DNA -b4=5 -b5=h    3
done

# concatenate sequences
perl scripts/BeforePhylo/BeforePhylo.pl -type=dna -fillgaps=100 -conc=raxml \
  multilocus_phylogeny/its_ITS1_rep_align.fasta-gb \
  multilocus_phylogeny/its_5_8S_rep_align.fasta-gb \
  multilocus_phylogeny/its_ITS2_rep_align.fasta-gb \
  multilocus_phylogeny/lsu_rep_align.fasta-gb
perl scripts/BeforePhylo/BeforePhylo.pl -type=dna -fillgaps=100 -conc=none \
  multilocus_phylogeny/its_ITS1_rep_align.fasta-gb \
  multilocus_phylogeny/its_5_8S_rep_align.fasta-gb \
  multilocus_phylogeny/its_ITS2_rep_align.fasta-gb \
  multilocus_phylogeny/lsu_rep_align.fasta-gb
mv output_partitions.txt multilocus_phylogeny/partitions.txt
mv output.fasta multilocus_phylogeny/concatenated_alignment.fasta
rm output.phy

# run phylogenetic analysis
raxmlHPC-SSE3 -s multilocus_phylogeny/concatenated_alignment.fasta -n ITS_LSU \
  -f a -m GTRGAMMA -q multilocus_phylogeny/partitions.txt -x 12345 -p 12345 \
  -# 100
mv RAxML_* multilocus_phylogeny/


### end
