#!/bin/bash
# This is a preprocessing for RNA-seq using RSEM for 4-read samples
# 11 paired-end sample with 2 replicates

# RSEM 1 & 2
for i in ~/Desktop/RSEM/data/4_reads/*_1.fq
do
j=`echo $i |sed "s/_1.fq/_2.fq/"`
ii=`echo $i |sed "s/_1.fq/_3.fq/"`
jj=`echo $i |sed "s/_1.fq/_4.fq/"`
k=${i:34:16}
printf "___________________________________________________________\n$i\n$j\n$ii\n$jj\n$k\n"
software/RSEM-1.3/rsem-calculate-expression -p 12 --paired-end --bowtie2 --bowtie2-path software/bowtie2-2.2 --estimate-rspd --append-names --sort-bam-by-coordinate --output-genome-bam $i $j ref/rat_ref exp4a/$k
done

# RSEM 3 & 4
for i in ~/Desktop/RSEM/data/4_reads/*_1.fq
do
date
j=`echo $i |sed "s/_1.fq/_2.fq/"`
ii=`echo $i |sed "s/_1.fq/_3.fq/"`
jj=`echo $i |sed "s/_1.fq/_4.fq/"`
k=${i:34:16}
printf "___________________________________________________________\n$i\n$j\n$ii\n$jj\n$k\n"
software/RSEM-1.3/rsem-calculate-expression -p 12 --paired-end --bowtie2 --bowtie2-path software/bowtie2-2.2 --estimate-rspd --append-names --sort-bam-by-coordinate --output-genome-bam $ii $jj ref/rat_ref exp4b/$k
done

# Concatication
for file in ~/Desktop/RSEM/data/4_reads/*_1.fq; do
i=${file##*/}
j=`echo $i |sed "s/_1.fq/_2.fq/"`
ii=`echo $i |sed "s/_1.fq/_3.fq/"`
jj=`echo $i |sed "s/_1.fq/_4.fq/"`
iicat=`echo $i |sed "s/_1.fq/_F.fq/"`
jjcat=`echo $i |sed "s/_1.fq/_R.fq/"`
k=${i:0:16}
echo "-------------------------" $k "---------------"
echo $i
echo $j
echo $ii
echo $jj
echo $iicat
echo $jjcat
cat $i $ii > $iicat
cat $j $jj > $jjcat
../../software/RSEM-1.3/rsem-calculate-expression -p 12 --paired-end --bowtie2 --bowtie2-path ../../software/bowtie2-2.2 --estimate-rspd --append-names --sort-bam-by-coordinate --output-genome-bam $iicat $jjcat ../../ref/rat_ref ../../exp4cat/$k
done



# R
data = read.table("R43_s_00077E7B8A.genes.results", header=T, stringsAsFactors=F)
idx = order(data[,"TPM"], decreasing=T)
data[idx[1:10], c("gene_id", "expected_count", "TPM")]


# PDF for sample R43
../software/RSEM-1.3/rsem-plot-model R43_s_00077E7B8A R43_s_00077E7B8A_run1_diagnostic.pdf
../software/RSEM-1.3/rsem-plot-model R43_s_00077E7B8A R43_s_00077E7B8A_run2_diagnostic.pdf
../software/RSEM-1.3/rsem-plot-transcript-wiggles --gene-list --show-unique R43_s_00077E7B8A gene_ids.txt R43_s_00077E7B8A_run1_wiggle.pdf
../software/RSEM-1.3/rsem-plot-transcript-wiggles --gene-list --show-unique R43_s_00077E7B8A gene_ids.txt R43_s_00077E7B8A_run2_wiggle.pdf






