#######################################################################################################
### DESeq2 Analysis and finding the relation between transcriptomic and IOP information for 45 Rats ###
#######################################################################################################
library(dplyr)
library(tximport)
library(DESeq2)
library(vsn)
library(pheatmap)
library(RColorBrewer)

################################################
### Step 1: Import samples' information file ###
################################################
dir <- setwd("/home/h/Desktop/RSEM/export"); getwd()
# samples <- read.csv("samples.csv", sep=",", header=TRUE)
# samples_iop <- read.csv("p50_iop_all_data.csv", sep=",", header=TRUE, stringsAsFactors=TRUE, fill=TRUE)
# samples_iop[,1] <- gsub("#", "s_", samples_iop[,1])
# names(samples_iop)[1] <- "sample"
# write.table(samples_iop, file="p50_iop_all_data_modified.csv", sep=",")
# new_samples <- left_join(samples, samples_iop, by="sample")
# row.names(new_samples) <- row.names(samples)
# samples <- new_samples[complete.cases(new_samples), ]         # remove 8 samples which have no IOP data ("NA")
# for (i in (1:53)) {
#   samples$OS[i] <- mean(samples$OS1[i], samples$OS2[i], samples$OS3[i])
#   samples$OD[i] <- mean(samples$OD1[i], samples$OD2[i], samples$OD3[i])
#   samples$IOP[i] <- mean(samples$OD1[i], samples$OD2[i], samples$OD3[i], samples$OS1[i], samples$OS2[i], samples$OS3[i])
#   samples$Class[i] <- if(samples$IOP[i] <= 15) "N" else if(15 < samples$IOP[i] & samples$IOP[i] <= 20) "E" else "H" }
# write.table(samples, file="samples_modified.csv", sep=",")
samples <- read.csv("samples_modified.csv", sep=",", header=TRUE)
dim(samples)
head(samples)
#################### Plot Pie Chart for classes ###################
# par(mar=c(0.05, 0.05, 0.05, 0.05))
# slices <- c(count(samples$Class=="N"), count(samples$Class=="E"),count(samples$Class=="H"))
# lbls <- paste(c("Normal ", "Elevated ", "High "), round(slices/sum(slices)*100), "%", sep="") 
# pie(slices, labels=lbls, col=rainbow(3))

##########################################################################################
### Step 2: Import RSEM outputs(transcripts aboundance estimation) file using tximport ###
##########################################################################################
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
dds <- DESeqDataSetFromTximport(txi.rsem, samples, ~Class)
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
# ntd <- normTransform(dds)
# meanSdPlot(assay(ntd))
# meanSdPlot(assay(vsd))
meanSdPlot(assay(rld), bins = 100)

selectedGenes <- data.frame(t(assay(rld)[c("ENSRNOG00000004712_Angptl1", "ENSRNOG00000055293_Ptprb",
                                "ENSRNOG00000016696_Angpt2",  "ENSRNOG00000020603_Angptl6",
                                "ENSRNOG00000020173_Tie1",    "ENSRNOG00000008587_Tek",
                                "ENSRNOG00000025562_Ang2",    "ENSRNOG00000043032_Ang2"),]))
selectedGenes$sample <- samples$sample
selectedGenes <- selectedGenes[c(9,1,2,3,4,5,6,7,8)]
new_samples <- left_join(samples, selectedGenes, by="sample")
row.names(new_samples) <- row.names(samples)


###########################################
### Step 4: Heatmap of the count matrix ###
###########################################
heat.colors <- brewer.pal(6, "PuBu")
dds <- estimateSizeFactors(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("sample", "class")])
# pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
# pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=F, cluster_cols=TRUE, annotation_col=df)
pheatmap(cor(assay(rld)))
pheatmap(cor(assay(rld)), color=heat.colors, border_color=NA, fontsize=10, fontsize_row=10, height=20, main="Heatmap of sample-to-sample distances")
hist(assay(rld))
