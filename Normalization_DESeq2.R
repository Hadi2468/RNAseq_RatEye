#######################################################################################
### DESeq2 Analysis and Normalization of the transcriptomic data for 53 rat samples ###
#######################################################################################

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.10")
# BiocManager::install("tximport")
# BiocManager::install("tximportData")
# BiocManager::install("DESeq2")
# install.packages("tidyverse")
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("vctrs")
# install.packages("fs")

Packages <- c("tximport", "tximportData", "DESeq2", "tidyverse", "dplyr", "vctrs", "fs", "ggplot2")
lapply(Packages, library, character.only = TRUE)
# Packages <- c("vsn", "pheatmap", "RColorBrewer", "corrplot", "ggpubr")

###################################################
### Step 1: Import 53 samples' information file ###
###################################################
setwd("./Data"); getwd()
samples <- read.csv("samples53.csv", sep=",", header=TRUE)
# new_samples <- left_join(samples53, samples, by=c("sample", "SeqType"))
# row.names(samples) <- samples$sample
# for (i in (1:53)) {
#   samples$Avg_OD[i] <- as.integer(round( (samples$OD1[i] + samples$OD2[i] + samples$OD3[i]) / 3 ))
#   samples$Avg_OS[i] <- as.integer(round( (samples$OS1[i] + samples$OS2[i] + samples$OS3[i]) / 3 ))
#   samples$Avg_IOP[i] <- as.integer(round( (samples$OD1[i] + samples$OD2[i] + samples$OD3[i] + samples$OS1[i] + samples$OS2[i] + samples$OS3[i]) / 6 )) }
# for (i in (1:53)) {
#   if (samples$Avg_OD[i] == 0) {samples$Class_OD[i] <- "NoData"}
#   else if (0 < samples$Avg_OD[i] & samples$Avg_OD[i] < 15) {samples$Class_OD[i] <- "Normal"}
#   else if (15 <= samples$Avg_OD[i] & samples$Avg_OD[i] < 20) {samples$Class_OD[i] <- "Elevated"}
#   else {samples$Class_OD[i] <- "High"}
#   if (samples$Avg_OS[i] == 0) {samples$Class_OS[i] <- "NoData"}
#   else if (0 < samples$Avg_OS[i] & samples$Avg_OS[i] < 15) {samples$Class_OS[i] <- "Normal"}
#   else if (15 <= samples$Avg_OS[i] & samples$Avg_OS[i] < 20) {samples$Class_OS[i] <- "Elevated"}
#   else {samples$Class_OS[i] <- "High"}
#   if (samples$Avg_IOP[i] == 0) {samples$Class_IOP[i] <- "NoData"}
#   else if (0 < samples$Avg_IOP[i] & samples$Avg_IOP[i] < 15) {samples$Class_IOP[i] <- "Normal"}
#   else if (15 <= samples$Avg_IOP[i] & samples$Avg_IOP[i] < 20) {samples$Class_IOP[i] <- "Elevated"}
#   else {samples$Class_IOP[i] <- "High"} }
# for (i in (1:length(samples$Batch))) {
#   samples$Batch[i] <- if(samples$Batch[i] == "12") "B12" else if(samples$Batch[i] == "13") "B13" else if(samples$Batch[i] == "14") "B14" else "B15" }
# for (i in (1:length(samples$Class))) {
#   samples$Class[i] <- if(samples$Class[i] == "H") "High" else if(samples$Class[i] == "E") "Elevated" else "Normal" }
# samples <- new_samples[complete.cases(new_samples), ]         # remove 8 samples which have no IOP data ("NA")
# dim(samples)
# head(samples)
samples$Class_OD <- relevel(samples$Class_OD, "Normal")
samples$Class_OS <- relevel(samples$Class_OS, "Normal")
samples$Class_IOP <- relevel(samples$Class_IOP, "Normal")
# write.table(samples, file="samples53.csv", sep=",")
############# Plot Pie Chart for classes ###############
# par(mar=c(0.05, 0.05, 0.05, 0.05))
# slices <- c(count(samples$Class=="Normal"), count(samples$Class=="Elevated"),count(samples$Class=="High"))
# lbls <- paste(c("Normal ", "Elevated ", "High "), round(slices/sum(slices)*100), "%", sep="")
# pie(slices, labels=lbls, col=rainbow(3))

##########################################################################################
### Step 2: Import RSEM outputs(transcripts aboundance estimation) file using tximport ###
##########################################################################################
setwd("./rsem_export_dataset"); getwd()
RSEM_output_gene_files <- file.path(getwd(), paste0(samples$sample, ".genes.results"))
names(RSEM_output_gene_files) <- samples$sample
txi.rsem <- tximport(RSEM_output_gene_files, type="rsem", txIn=FALSE, txOut=FALSE)
# head(RSEM_output_gene_files)
# names(txi.rsem)
# row.names(txi.rsem$length) 
# head(txi.rsem$abundance)[,1:5]
# head(txi.rsem$counts)[,1:5]
# head(txi.rsem$length)[,1:5]
# ########### Preparing the isoform count estimation matrix##########
# RSEM_output_isoform_files <- file.path(getwd(), paste0(samples$sample, ".isoforms.results"))
# names(RSEM_output_isoform_files) <- samples$sample
# txi_i.rsem <- tximport(RSEM_output_isoform_files, type="rsem", txIn=TRUE, txOut=TRUE)
# head(txi_i.rsem$abundance)[,1:5]
# head(txi_i.rsem$counts)[,1:5]
# head(txi_i.rsem$length)[,1:5]

#########################################################################################
### Step 3: Producing the DESeq2 data frame based on the outputs of Step 1 and Step 2 ###
#########################################################################################
# colnames(txi.rsem$counts)
# row.names(samples)
txi.rsem$length[txi.rsem$length == 0] <- 1
dds <- DESeqDataSetFromTximport(txi.rsem, samples, ~Class_IOP)
# dds <- DESeqDataSetFromTximport(txi.rsem, samples, ~Class_OD)
# dds <- DESeqDataSetFromTximport(txi.rsem, samples, ~Class_OS)
# dds
# colData(dds)     # Green box
# rowData(dds)     # Blue box
# assay(dds)       # Pink box
# ############ Filtering ############## 
# cds <- dds
# dds <- cds
# nrow(dds)
# keep <- rowSums(counts(dds)) > 10   # remove genes without expression in more than 10 sample
# dds <- dds[keep,]
# nrow(dds)

#################################################
### Step 4: Normalization based on rlog & vsd ###
#################################################
setwd("../"); getwd()
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
# head(assay(rld), 3)
# ntd <- normTransform(dds)
# meanSdPlot(assay(ntd))
# meanSdPlot(assay(vsd))
# meanSdPlot(assay(rld), bins = 100)
# sub_genes_rld <- data.frame(t(assay(rld)))
# sub_genes_vst <- data.frame(t(assay(vsd)))
# dim(sub_genes_rld)
# dim(sub_genes_vst)
# sub_genes_rld <- sub_genes_rld[-c(3, 14, 42, 43, 44, 45, 46, 47), ]
# sub_genes_vst <- sub_genes_vst[-c(3, 14, 42, 43, 44, 45, 46, 47), ]
# write.table(sub_genes_rld, file="normalized_rlog_class_IOP.csv", sep=",", quote=F, row.names=TRUE, col.names=TRUE,)
# write.table(sub_genes_vst, file="normalized_vst_class_IOP.csv", sep=",", quote=F, row.names=TRUE, col.names=TRUE,)
sub_genes <- read.csv("normalized_rlog_class_IOP.csv", sep=",", header=TRUE)
sub_genes <- read.csv("normalized_vst_class_IOP.csv", sep=",", header=TRUE)
dim(sub_genes)
