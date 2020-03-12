#!/bin/bash
# This is a sample script for producing the PDF files based on the results of RSEM analysis for sample s_00077E7B8Athe 

clear 
../software/RSEM-1.3/rsem-plot-model s_00077E7B8A s_00077E7B8A_run1_diagnostic.pdf
../software/RSEM-1.3/rsem-plot-model s_00077E7B8A s_00077E7B8A_run2_diagnostic.pdf
../software/RSEM-1.3/rsem-plot-transcript-wiggles --gene-list --show-unique s_00077E7B8A gene_ids.txt s_00077E7B8A_run1_wiggle.pdf
../software/RSEM-1.3/rsem-plot-transcript-wiggles --gene-list --show-unique s_00077E7B8A gene_ids.txt s_00077E7B8A_run2_wiggle.pdf
