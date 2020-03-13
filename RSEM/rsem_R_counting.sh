#!/bin/bash
# This is the R script to counting the genes based on the results of RSEM analysis


clear 
TopGenes = 10
GeneResultsFiles <- list.files(pattern="*.genes.results")
for (sample in 1:length(GeneResultsFiles)){
	print("---------------------------------------------------------------")
	print(GeneResultsFiles[sample])
	data = read.table(GeneResultsFiles[sample], header=T, stringsAsFactors=F)
	idx = order(data[,"TPM"], decreasing=T)
	print(data[idx[1:TopGenes], c("gene_id", "expected_count", "TPM")])}

