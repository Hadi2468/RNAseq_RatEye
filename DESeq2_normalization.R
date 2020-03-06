############################################################################################################
### Step 1: Import RSEM outputs(transcripts aboundance estimation) and samples_name files using tximport ###
############################################################################################################
dir <- setwd("/home/h/Desktop/RSEM/export"); getwd()
samples <- read.table(file.path(dir, "samples.csv"), sep=",", header=TRUE)
# samples$condition <- factor(c(rep("paired-end", each=3), rep("2_paired-ends", each=50)))
# write.table(samples, file="samples.csv", sep=",")
head(samples)
row.names(samples)

RSEM_output_gene_files <- file.path(dir, paste0(samples$sample, ".genes.results"))
names(RSEM_output_gene_files) <- samples$sample
head(RSEM_output_gene_files)
txi.rsem <- tximport(RSEM_output_gene_files, type="rsem", txIn=FALSE, txOut=FALSE)
names(txi.rsem)
row.names(txi.rsem$length) 
head(txi.rsem$abundance)[,1:5]
head(txi.rsem$counts)[,1:5]
head(txi.rsem$length)[,1:5]
geneMat <- (txi.rsem$counts)    # geneMat is the estimation of gene counts (integers)
head(geneMat)
dim(geneMat)
colnames(geneMat)
row.names(geneMat)

# ########################### Preparing the isoform count estimation matrix#############################
# RSEM_output_isoform_files <- file.path(dir, paste0(samples$sample, ".isoforms.results"))
# names(RSEM_output_isoform_files) <- samples$sample
# txi_i.rsem <- tximport(RSEM_output_isoform_files, type="rsem", txIn=TRUE, txOut=TRUE)
# head(txi_i.rsem$abundance)[,1:5]
# head(txi_i.rsem$counts)[,1:5]
# head(txi_i.rsem$length)[,1:5]
# isoformMat <- (txi_i.rsem$counts)    # isoformMat is the estimation of isoform counts (non-integers)
# head(isoformMat)
# dim(isoformMat)

#####################################################
### Step 2: Produsing the DESeq2-based data frame ###
#####################################################
library(DESeq2)
head(txi.rsem$counts)[,1:5]
colnames(txi.rsem$counts)
row.names(samples)
txi.rsem$length[txi.rsem$length == 0] <- 1
dds <- DESeqDataSetFromTximport(txi.rsem, samples, ~condition)
dds
colData(dds)     # Green box
rowData(dds)     # Blue box
assay(dds)       # Pink box

# ############ Filtering ############## 
# cds <- dds
dds <- cds
# nrow(dds)
# keep <- rowSums(counts(dds)) > 10       # remove genes without expression in more than 10 sample
# dds <- dds[keep,]
# nrow(dds)

###########################################
### Step 3: Normalization based on rlog ###
###########################################
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(rld), 3)
write.table(assay(rld), file="normalized_counts.csv", sep=",", quote=F, col.names=NA)

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
