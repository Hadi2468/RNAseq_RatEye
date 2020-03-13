##############################################################################################
### DESeq2 Analysis and finding the corrolation between transcriptomic and IOP for 45 Rats ###
##############################################################################################

# BiocManager::install("")
# install.packages("corrplot")
# install.packages("RColorBrewer")
# install.packages("ggpubr")
# install.packages("pheatmap")
# install.packages("ppcor")

Packages <- c("tximport", "tximportData", "DESeq2", "tidyverse", "dplyr", "vctrs", "fs", "ggplot2", 
              "corrplot", "RColorBrewer", "ggpubr", "pheatmap", "ppcor")
lapply(Packages, library, character.only = TRUE)

####################################################
### Step 1: Import 45 samples' information files ###
####################################################

getwd(); setwd("./Data"); getwd()
samples <- read.csv("samples53.csv", sep=",", header=TRUE)
dim(samples)
# sub_samples <- samples[-c(3, 14, 42, 43, 44, 45, 46, 47), ]
# write.table(sub_samples, file="samples45.csv", sep=",", quote=F, row.names=TRUE, col.names=TRUE,)
sub_samples <- read.csv("samples45.csv", sep=",", header=TRUE)
dim(sub_samples)
sub_genes <- read.csv("normalized_rlog_class_IOP.csv", sep=",", header=TRUE)
sub_genes <- read.csv("normalized_vst_class_IOP.csv", sep=",", header=TRUE)
dim(sub_genes)

####################################
### Step 2: Correlation Analysis ###
####################################

selGenes <- subset(sub_genes, select=c(ENSRNOG00000016696_Angpt2, ENSRNOG00000055293_Ptprb, ENSRNOG00000008587_Tek))
# corrTable <- cbind(sub_samples$Avg_OD, sub_samples$Avg_OS, sub_samples$Avg_IOP, selGenes)    # Correlation tables for three genes 
# names(corrTable) <- c("_OD", "_OS", "_IOP", "ANGPT2", "PTPRB", "TEK")
corrTable <- cbind(sub_samples$Avg_IOP, selGenes)    # Correlation tables for three genes 
names(corrTable) <- c("IOP", "ANGPT2", "PTPRB", "TEK")

corrplot(cor(corrTable, method="pearson"), method="color", type="upper", order="hclust", 
         col=colorRampPalette(c("dodgerblue", "aliceblue", "brown1"))(7), addCoef.col="black", tl.col="black", tl.cex=1, addrect=3)
corrplot(cor(corrTable, method="spearman"), method="color", type="upper", order="hclust", 
         col=colorRampPalette(c("dodgerblue", "aliceblue", "brown1"))(7), addCoef.col="black", tl.col="black", tl.cex=1, addrect=3)
# corrplot(cor(corrTable, method="kendall"), method="color", type="upper", order="hclust", 
#          col=colorRampPalette(c("dodgerblue", "aliceblue", "brown1"))(7), addCoef.col="black", tl.col="black", tl.cex=1, addrect=3)
pheatmap(cor(corrTable))

####### "Sex" based scatter plot ############
corrTable_new <- cbind(sub_samples$Sex, corrTable)    # Correlation tables for three genes 
# names(corrTable_new) <- c("Sex", "OD", "OS", "IOP", "ANGPT2", "PTPRB", "TEK")
names(corrTable_new) <- c("Sex", "IOP", "ANGPT2", "PTPRB", "TEK")

ggscatter(corrTable_new, x="IOP", y=c("ANGPT2", "PTPRB", "TEK"), size=3, shape=19, color="Sex", cor.method="pearson", title="Correlation: Pearson,    Normalization: rlog",
          combine = TRUE, add="reg.line", conf.int=F, cor.coef=TRUE, xlab="IOP", ylab="Normalized gene counts", palette=c("red", "pink", "blue"))

# ggscatter(corrTable_new, x="OD", y=c("PTPRB", "ANGPT2", "TEK"), size=3, shape=19, color="Sex", cor.method="pearson", title="Correlation: Pearson, Normalization: rlog",
#           combine = TRUE, add="reg.line", conf.int=F, cor.coef=TRUE, xlab="OD", ylab="Normalized gene counts", palette=c("red", "pink", "blue")) 
# ggscatter(corrTable_new, x="OS", y=c("PTPRB", "ANGPT2", "TEK"), size=3, shape=19, color="Sex", cor.method="pearson", title="Correlation: Pearson, Normalization: rlog",
#           combine = TRUE, add="reg.line", conf.int=F, cor.coef=TRUE, xlab="OS", ylab="Normalized gene counts", palette=c("red", "pink", "blue")) 

# ggscatter(corrTable_new, x="OD", y=c("PTPRB", "ANGPT2", "TEK"), size=3, shape=19, color="Sex", cor.method="spearman", title="Correlation: Spearman, Normalization: rlog",
#           combine = TRUE, add="reg.line", conf.int=F, cor.coef=TRUE, xlab="OD", ylab="Normalized gene counts", palette=c("red", "pink", "blue")) 
# ggscatter(corrTable_new, x="OS", y=c("PTPRB", "ANGPT2", "TEK"), size=3, shape=19, color="Sex", cor.method="spearman", title="Correlation: Spearman, Normalization: rlog",
#           combine = TRUE, add="reg.line", conf.int=F, cor.coef=TRUE, xlab="OS", ylab="Normalized gene counts", palette=c("red", "pink", "blue")) 

# ggscatter(corrTable_new, x="OD", y=c("PTPRB", "ANGPT2", "TEK"), size=3, shape=19, color="Sex", cor.method="kendall", title="Correlation: Kendall, Normalization: rlog",
#           combine = TRUE, add="reg.line", conf.int=F, cor.coef=TRUE, xlab="OD", ylab="Normalized gene counts", palette=c("red", "pink", "blue")) 
# ggscatter(corrTable_new, x="OS", y=c("PTPRB", "ANGPT2", "TEK"), size=3, shape=19, color="Sex", cor.method="kendall", title="Correlation: Kendall, Normalization: rlog",
#           combine = TRUE, add="reg.line", conf.int=F, cor.coef=TRUE, xlab="OS", ylab="Normalized gene counts", palette=c("red", "pink", "blue")) 
# ggscatter(corrTable_new, x="IOP", y=c("PTPRB", "ANGPT2", "TEK"), size=3, shape=19, color="Sex", cor.method="kendall", title="Correlation: Kendall, Normalization: rlog",
#           combine = TRUE, add="reg.line", conf.int=F, cor.coef=TRUE, xlab="IOP", ylab="Normalized gene counts", palette=c("red", "pink", "blue"))

ggscatter(corrTable_new, x="ANGPT2", y=c("PTPRB", "TEK"), size=3, shape=19, color="Sex", cor.method="pearson", title="Correlation: Pearson,    Normalization: rlog",
          combine = TRUE, add="reg.line", conf.int=FALSE, cor.coef=TRUE, xlab="ANGPT2", ylab="rlog-normalized gene counts", palette=c("red", "pink", "blue"))
ggscatter(corrTable_new, x="PTPRB", y="TEK", size=3, shape=19, color="Sex", cor.method="pearson", title="Correlation: Pearson,    Normalization: rlog",
          combine = TRUE, add="reg.line", conf.int=FALSE, cor.coef=TRUE, xlab="PTPRB", ylab="TEK", palette=c("red", "pink", "blue"))

ggscatter(corrTable_new, x="ANGPT2", y=c("PTPRB", "TEK"), size=3, shape=19, color="Sex", cor.method="spearman", title="Correlation: Spearman,    Normalization: rlog",
          combine = TRUE, add="reg.line", conf.int=FALSE, cor.coef=TRUE, xlab="ANGPT2", ylab="rlog-normalized gene counts", palette=c("red", "pink", "blue"))
ggscatter(corrTable_new, x="PTPRB", y="TEK", size=3, shape=19, color="Sex", cor.method="spearman", title="Correlation: Spearman,    Normalization: rlog",
          combine = TRUE, add="reg.line", conf.int=FALSE, cor.coef=TRUE, xlab="PTPRB", ylab="TEK", palette=c("red", "pink", "blue"))

####### "Age" based scatter plot ############
corrTable_new <- cbind(sub_samples$Class_Age, corrTable)    # Correlation tables for three genes 
names(corrTable_new)[1] <- "Age"

ggscatter(corrTable_new, x="IOP", y=c("ANGPT2", "PTPRB", "TEK"), size=3, shape=19, color="Age", cor.method="pearson", title="Correlation: Pearson,    Normalization: rlog",
          combine = TRUE, add="reg.line", conf.int=F, cor.coef=TRUE, xlab="IOP", ylab="rlog-normalized gene counts", palette=c("green", "red", "blue", "black"))
ggscatter(corrTable_new, x="ANGPT2", y=c("PTPRB", "TEK"), size=3, shape=19, color="Age", cor.method="pearson", title="Correlation: Pearson,    Normalization: rlog",
          combine = TRUE, add="reg.line", conf.int=FALSE, cor.coef=TRUE, xlab="ANGPT2", ylab="rlog-normalized gene counts", palette=c("green", "red", "blue", "black")) 
ggscatter(corrTable_new, x="PTPRB", y="TEK", size=3, shape=19, color="Age", cor.method="pearson", title="Correlation: Pearson,    Normalization: rlog",
          combine = TRUE, add="reg.line", conf.int=FALSE, cor.coef=TRUE, xlab="PTPRB", ylab="TEK", palette=c("green", "red", "blue", "black"))

ggscatter(corrTable_new, x="IOP", y=c("ANGPT2", "PTPRB", "TEK"), size=3, shape=19, color="Age", cor.method="spearman", title="Correlation: Spearman,    Normalization: rlog",
          combine = TRUE, add="reg.line", conf.int=F, cor.coef=TRUE, xlab="IOP", ylab="rlog-normalized gene counts", palette=c("green", "red", "blue", "black"))
ggscatter(corrTable_new, x="ANGPT2", y=c("PTPRB", "TEK"), size=3, shape=19, color="Age", cor.method="spearman", title="Correlation: Spearman,    Normalization: rlog",
          combine = TRUE, add="reg.line", conf.int=FALSE, cor.coef=TRUE, xlab="ANGPT2", ylab="rlog-normalized gene counts", palette=c("green", "red", "blue", "black")) 
ggscatter(corrTable_new, x="PTPRB", y="TEK", size=3, shape=19, color="Age", cor.method="spearman", title="Correlation: Spearman,    Normalization: rlog",
          combine = TRUE, add="reg.line", conf.int=FALSE, cor.coef=TRUE, xlab="PTPRB", ylab="TEK", palette=c("green", "red", "blue", "black"))

####### "Batch" based scatter plot ############
corrTable_new <- cbind(sub_samples$Batch, corrTable)    # Correlation tables for three genes 
names(corrTable_new)[1] <- "Batch"

ggscatter(corrTable_new, x="IOP", y=c("ANGPT2", "PTPRB", "TEK"), size=3, shape=19, color="Batch", cor.method="pearson", title="Correlation: Pearson,    Normalization: rlog",
          combine = TRUE, add="reg.line", conf.int=F, cor.coef=TRUE, xlab="IOP", ylab="rlog-normalized gene counts", palette=c("green", "red", "blue", "black"))
ggscatter(corrTable_new, x="ANGPT2", y=c("PTPRB", "TEK"), size=3, shape=19, color="Batch", cor.method="pearson", title="Correlation: Pearson,    Normalization: rlog",
          combine = TRUE, add="reg.line", conf.int=FALSE, cor.coef=TRUE, xlab="ANGPT2", ylab="rlog-normalized gene counts", palette=c("green", "red", "blue", "black")) 
ggscatter(corrTable_new, x="PTPRB", y="TEK", size=3, shape=19, color="Batch", cor.method="pearson", title="Correlation: Pearson,    Normalization: rlog",
          combine = TRUE, add="reg.line", conf.int=FALSE, cor.coef=TRUE, xlab="PTPRB", ylab="TEK", palette=c("green", "red", "blue", "black"))

ggscatter(corrTable_new, x="IOP", y=c("ANGPT2", "PTPRB", "TEK"), size=3, shape=19, color="Batch", cor.method="spearman", title="Correlation: Spearman,    Normalization: rlog",
          combine = TRUE, add="reg.line", conf.int=F, cor.coef=TRUE, xlab="IOP", ylab="rlog-normalized gene counts", palette=c("green", "red", "blue", "black"))
ggscatter(corrTable_new, x="ANGPT2", y=c("PTPRB", "TEK"), size=3, shape=19, color="Batch", cor.method="spearman", title="Correlation: Spearman,    Normalization: rlog",
          combine = TRUE, add="reg.line", conf.int=FALSE, cor.coef=TRUE, xlab="ANGPT2", ylab="rlog-normalized gene counts", palette=c("green", "red", "blue", "black")) 
ggscatter(corrTable_new, x="PTPRB", y="TEK", size=3, shape=19, color="Batch", cor.method="spearman", title="Correlation: Spearman,    Normalization: rlog",
          combine = TRUE, add="reg.line", conf.int=FALSE, cor.coef=TRUE, xlab="PTPRB", ylab="TEK", palette=c("green", "red", "blue", "black"))
           
####### Partial correlation test ####### 
corrTable <- cbind(sub_samples$Avg_IOP, sub_samples$AgeInDays, sub_samples$Batch, selGenes)    # Correlation tables for three genes 
names(corrTable) <- c("IOP", "Age", "Batch", "ANGPT2", "PTPRB", "TEK")

pcor.test(corrTable$ANGPT2, corrTable$IOP, corrTable$Age, method="pearson")
summary(lm(corrTable$ANGPT2~corrTable$IOP+corrTable$Age))

pcor.test(corrTable$ANGPT2, corrTable$IOP, corrTable$Batch, method="pearson")
summary(lm(corrTable$ANGPT2~corrTable$IOP+corrTable$Age+corrTable$Batch))

summary(lm(corrTable$ANGPT2~corrTable$IOP))

########################################
### Others helpful codes for memorize ###
########################################
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

samples$AgeInDays[3] <- 0L
samples$AgeInDays[14] <- 0L
samples$AgeInDays[42] <- 0L
samples$AgeInDays[43] <- 0L
samples$AgeInDays[44] <- 0L
samples$AgeInDays[45] <- 0L
samples$AgeInDays[46] <- 0L
samples$AgeInDays[47] <- 0L