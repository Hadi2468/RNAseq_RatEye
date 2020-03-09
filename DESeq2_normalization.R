############################################################################################################
### Step 1: Import RSEM outputs(transcripts aboundance estimation) and samples_name files using tximport ###
############################################################################################################
library(dplyr)
library(tximport)
library(DESeq2)

dir <- setwd("/home/h/Desktop/RSEM/export"); getwd()
# samples <- read.csv("samples.csv", sep=",", header=TRUE)
# samples_iop <- read.csv("p50_iop_all_data.csv", sep=",", header=TRUE, stringsAsFactors=TRUE, fill=TRUE)
# samples_iop[,1] <- gsub("#", "s_", samples_iop[,1])
# names(samples_iop)[1] <- "sample"
# write.table(samples_iop, file="p50_iop_all_data_modified.csv", sep=",")
# new_samples <- left_join(samples, samples_iop, by="sample")
# row.names(new_samples) <- row.names(samples)
# write.table(new_samples, file="samples_modified.csv", sep=",")
# samples <- read.csv("samples_modified.csv", sep=",", header=TRUE)
samples <- read.csv("samples_test_labled.csv", sep=",", header=TRUE)
head(samples)
samples <- na.omit(samples)
dim(samples)
head(samples)
row.names(samples)
dir <- setwd("/home/h/Desktop/RSEM/export"); getwd()
RSEM_output_gene_files <- file.path(dir, paste0(samples$sample, ".genes.results"))
names(RSEM_output_gene_files) <- samples$sample
head(RSEM_output_gene_files)
txi.rsem <- tximport(RSEM_output_gene_files, type="rsem", txIn=FALSE, txOut=FALSE)
names(txi.rsem)
row.names(txi.rsem$length) 
head(txi.rsem$abundance)[,1:5]
head(txi.rsem$counts)[,1:5]
head(txi.rsem$length)[,1:5]

# ########################### Preparing the isoform count estimation matrix#############################
# RSEM_output_isoform_files <- file.path(dir, paste0(samples$sample, ".isoforms.results"))
# names(RSEM_output_isoform_files) <- samples$sample
# txi_i.rsem <- tximport(RSEM_output_isoform_files, type="rsem", txIn=TRUE, txOut=TRUE)
# head(txi_i.rsem$abundance)[,1:5]
# head(txi_i.rsem$counts)[,1:5]
# head(txi_i.rsem$length)[,1:5]

#####################################################
### Step 2: Produsing the DESeq2-based data frame ###
#####################################################
colnames(txi.rsem$counts)
row.names(samples)
txi.rsem$length[txi.rsem$length == 0] <- 1
# dds <- DESeqDataSetFromTximport(txi.rsem, samples, ~read_type)
dds <- DESeqDataSetFromTximport(txi.rsem, samples, ~class)

dds
colData(dds)     # Green box
rowData(dds)     # Blue box
assay(dds)       # Pink box

# ############ Filtering ############## 
# cds <- dds
# dds <- cds
# nrow(dds)
# keep <- rowSums(counts(dds)) > 10       # remove genes without expression in more than 10 sample
# dds <- dds[keep,]
# nrow(dds)

###########################################
### Step 3: Normalization based on rlog ###
###########################################
# vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(rld), 3)
  # write.table(assay(rld), file="normalized_counts.csv", sep=",", quote=F, col.names=NA)

ntd <- normTransform(dds)
library(vsn)
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

###########################################
### Step 4: Heatmap of the count matrix ###
###########################################
library(pheatmap)
library(RColorBrewer)
heat.colors <- brewer.pal(6, "PuBu")

dds <- estimateSizeFactors(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("sample", "condition")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=F, cluster_cols=TRUE, annotation_col=df)
pheatmap(cor(assay(rld)))
pheatmap(cor(assay(rld)), color=heat.colors, border_color=NA, fontsize=10, fontsize_row=10, height=20, main="Heatmap of sample-to-sample distances")
hist(assay(rld))
