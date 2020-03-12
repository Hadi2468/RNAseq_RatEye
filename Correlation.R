##############################################################################################
### DESeq2 Analysis and finding the corrolation between transcriptomic and IOP for 45 Rats ###
##############################################################################################

# BiocManager::install("")
# install.packages("corrplot")
install.packages("RColorBrewer")

Packages <- c("tximport", "tximportData", "DESeq2", "tidyverse", "dplyr", "vctrs", "fs", "ggplot2", "corrplot", "RColorBrewer")
lapply(Packages, library, character.only = TRUE)
# Packages <- c("vsn", "pheatmap", "ggpubr")

####################################################
### Step 1: Import 45 samples' information files ###
####################################################
setwd("./datasets"); getwd()
samples <- read.csv("samples53.csv", sep=",", header=TRUE)
dim(samples)
# sub_samples <- samples[-c(3, 14, 42, 43, 44, 45, 46, 47), ]
# write.table(sub_samples, file="samples_subset.csv", sep=",", quote=F, row.names=TRUE, col.names=TRUE,)
sub_samples <- read.csv("samples_subset.csv", sep=",", header=TRUE)
dim(sub_samples)
sub_genes <- read.csv("normalized_subset_IOP_Class.csv", sep=",", header=TRUE)
dim(sub_genes)

####################################
### Step 2: Correlation Analysis ###
####################################
selGenes <- subset(sub_genes, select=c(ENSRNOG00000055293_Ptprb, ENSRNOG00000016696_Angpt2, ENSRNOG00000008587_Tek))
corrTable <- cbind(sub_samples$Avg_OD, sub_samples$Avg_OS, sub_samples$Avg_IOP, selGenes)    # Correlation tables for three genes 
names(corrTable) <- c("OD", "OS", "IOP", "PTPRB", "AngPt2", "TEK")
corrplot(cor(corrTable), method="color", type="lower", order="hclust", col=brewer.pal(n=7, name="PRGn"),
         addCoef.col = "black", tl.col="black", tl.cex=1)
ggscatter(corrTable, x="IOP", y=c("PTPRB", "AngPt2", "TEK"), size = 1, color = "Blue", cor.method="pearson", 
          combine = TRUE, add="reg.line", conf.int=TRUE, cor.coef=TRUE, xlab="IOP", ylab="Normalized gene counts") 
ggscatter(corrTable, x="PTPRB", y="TEK", size = 1, color = "Blue", cor.method="pearson", 
          combine = TRUE, add="reg.line", conf.int=TRUE, cor.coef=TRUE, xlab="PTPRB", ylab="TEK") 
ggscatter(corrTable, x="PTPRB", y="AngPt2", size = 1, color = "Red", cor.method="pearson", 
          combine = TRUE, add="reg.line", conf.int=TRUE, cor.coef=TRUE, xlab="PTPRB", ylab="AngPt2") 

###########################################
### Step 6: Heatmap of the count matrix ###
###########################################
gene_subset <- c("ENSRNOG00000004712_Angptl1", "ENSRNOG00000055293_Ptprb",
                 "ENSRNOG00000016696_Angpt2",  "ENSRNOG00000020603_Angptl6",
                 "ENSRNOG00000020173_Tie1",    "ENSRNOG00000008587_Tek",
                 "ENSRNOG00000025562_Ang2",    "ENSRNOG00000043032_Ang2")
txi.rsem$abundance <- txi.rsem$abundance[gene_subset,]
txi.rsem$counts <- txi.rsem$counts[gene_subset,]
txi.rsem$length <- txi.rsem$length[gene_subset,]

heat.colors <- brewer.pal(6, "PuBu")
dds <- estimateSizeFactors(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("sample", "Class")])
# pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
# pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=F, cluster_cols=TRUE, annotation_col=df)
pheatmap(cor(assay(rld)))
pheatmap(cor(assay(rld)), color=heat.colors, border_color=NA, fontsize=10, fontsize_row=10, height=20, main="Heatmap of sample-to-sample distances")
hist(assay(rld))





samples$OS1[3] <- 0L
samples$OS1[14] <- 0L
samples$OS1[42] <- 0L
samples$OS1[43] <- 0L
samples$OS1[44] <- 0L
samples$OS1[45] <- 0L
samples$OS1[46] <- 0L
samples$OS1[47] <- 0L