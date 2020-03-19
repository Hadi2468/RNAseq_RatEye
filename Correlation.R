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

getwd()
setwd("./Data"); getwd()
samples <- read.csv("samples53.csv", sep=",", header=TRUE)
dim(samples)
# sub_samples <- samples[-c(3, 14, 42, 43, 44, 45, 46, 47), ]
# write.table(sub_samples, file="samples45.csv", sep=",", quote=F, row.names=TRUE, col.names=TRUE,)
sub_samples <- read.csv("samples45.csv", sep=",", header=TRUE)
dim(sub_samples)
# sub_genes_2 <- read.csv("normalized_log2.csv", sep=",", header=TRUE)
sub_genes_r <- read.csv("normalized_rlog.csv", sep=",", header=TRUE)
# sub_genes_v <- read.csv("normalized_vst.csv", sep=",", header=TRUE)
# sub_genes_c <- read.csv("real_counts.csv", sep=",", header=TRUE)
# dim(sub_genes_2)
dim(sub_genes_r)
# dim(sub_genes_v)
# dim(sub_genes_c)

###########################################################################
### Step 2: Correlation Analysis: Heatmap, Clustering, and Distribution ###
###########################################################################

# selGenes <- subset(sub_genes_2, select=c(ENSRNOG00000016696_Angpt2, ENSRNOG00000055293_Ptprb, ENSRNOG00000008587_Tek))
selGenes <- subset(sub_genes_r, select=c(ENSRNOG00000016696_Angpt2, ENSRNOG00000055293_Ptprb, ENSRNOG00000008587_Tek))
# selGenes <- subset(sub_genes_v, select=c(ENSRNOG00000016696_Angpt2, ENSRNOG00000055293_Ptprb, ENSRNOG00000008587_Tek))
# selGenes <- subset(sub_genes_c, select=c(ENSRNOG00000016696_Angpt2, ENSRNOG00000055293_Ptprb, ENSRNOG00000008587_Tek))

corrTable <- cbind(sub_samples$Avg_IOP, selGenes)    # Correlation tables for three genes 
names(corrTable) <- c("IOP", "ANGPT2", "PTPRB", "TEK")
summary(corrTable)     # Basic statistical analysis

pheatmap(cor(corrTable))

######## Pearson Correlatoin ########

corrplot(cor(corrTable, method="pearson"), method="color", type="upper", order="hclust", 
         col=colorRampPalette(c("dodgerblue", "aliceblue", "brown1"))(7), 
         addCoef.col="black", tl.col="black", tl.cex=1, addrect=3)

######## Spearman Correlatoin ########

corrplot(cor(corrTable, method="spearman"), method="color", type="upper", order="hclust", 
         col=colorRampPalette(c("dodgerblue", "aliceblue", "brown1"))(7), 
         addCoef.col="black", tl.col="black", tl.cex=1, addrect=3)

######## Probability Density Function ########

ThreeGenes <- read.csv("ThreeGenes.csv", sep=",", header=TRUE)
ThreeGenes$Expression[1:45] <- corrTable$ANGPT2
ThreeGenes$Expression[46:90] <- corrTable$PTPRB
ThreeGenes$Expression[91:135] <- corrTable$TEK
ggplot(ThreeGenes, aes(Expression, fill=Gene)) + geom_density(alpha=0.6) +   scale_fill_manual(values=c("red", "blue", "yellow"))
ggplot(corrTable, aes(IOP)) + geom_density(fill="green") + scale_x_continuous(limits=c(10, 25))

###################################################
### Step 3: Correlation Analysis: Scatter Plots ###
###################################################

corrTable <- cbind(sub_samples$Class_IOP, corrTable)    # Correlation tables for three genes and three groups of IOP 
names(corrTable)[1] <- "Class_IOP"
levels(corrTable$Class_IOP)
corrTable$Class_IOP <- relevel(corrTable$Class_IOP, "Normal")
levels(corrTable$Class_IOP)

corrTable_IOP2 <- corrTable           
for (i in (1:45)) {if (corrTable_IOP2$Class_IOP[i] == "Elevated") {corrTable_IOP2$Class_IOP[i] <- "High"}}
corrTable_IOP2$Class_IOP <- factor(corrTable_IOP2$Class_IOP)    # Correlation tables for just two groups of IOP
levels(corrTable_IOP2$Class_IOP)

ggscatter(corrTable, x="IOP", y=c("ANGPT2", "PTPRB", "TEK"), size=3, shape=19, color="Class_IOP", 
          cor.method="spearman", title="Correlation: Spearman,    Normalization: rlog", combine=TRUE, 
          add="reg.line", conf.int=TRUE, cor.coef=TRUE, xlab="IOP", ylab="Gene expression") + 
  stat_cor(aes(color=Class_IOP), label.x=3)

ggscatter(corrTable_IOP2, x="IOP", y=c("ANGPT2", "PTPRB", "TEK"), size=3, shape=19, color="Class_IOP", 
          cor.method="spearman", title="Correlation: Spearman,    Normalization: rlog", combine=TRUE, 
          add="reg.line", conf.int=TRUE, cor.coef=FALSE, xlab="IOP", ylab="Gene expression") + 
  stat_cor(aes(color=Class_IOP), label.x=3) 

######### Scatter plot using ggplot2 ##########

ggANG <- ggplot(corrTable, aes(x=IOP, y=ANGPT2, color=Class_IOP)) + scale_size_manual(values=c(2,2,2)) +
  geom_point(aes(size=Class_IOP)) + geom_smooth(method=lm, aes(fill=Class_IOP), se=FALSE, fullrange=TRUE) + 
  theme(legend.position="top") + labs(y="Expression")
ggPT <- ggplot(corrTable, aes(x=IOP, y=PTPRB, color=Class_IOP)) + scale_size_manual(values=c(2,2,2)) +
  geom_point(aes(size=Class_IOP)) + geom_smooth(method=lm, aes(fill=Class_IOP), se=FALSE, fullrange=TRUE) + 
  theme(legend.position="top", axis.title.y = element_blank())
ggTEK <- ggplot(corrTable, aes(x=IOP, y=TEK, color=Class_IOP)) + scale_size_manual(values=c(2,2,2)) +
  geom_point(aes(size=Class_IOP)) + geom_smooth(method=lm, aes(fill=Class_IOP), se=FALSE, fullrange=TRUE) + 
  theme(legend.position="top", axis.title.y = element_blank())
ggarrange(ggANG, ggPT, ggTEK + rremove("x.text"), labels=c("ANGPT2", "PTPRB", "TEK"), ncol=3, nrow=1)

#########################################
### Step 4: t-test and ANOVA Analysis ###
#########################################

############ t-test Analysis ############

# H0: Mean of ANGPT2 expression for samples with Normal_IOP = Mean of ANGPT2 expression for samples with High_IOP 
ggplot(corrTable_IOP2, aes(x=Class_IOP, y=ANGPT2)) + geom_boxplot(color="blue", fill="pink", size=1.2)
t.test(ANGPT2 ~ Class_IOP, data=corrTable_IOP2, mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)

############ ANOVA Analysis #############

# H0: Mean of ANGPT2 expression for samples with all IOP gropus are equal 
ggplot(corrTable, aes(x=Class_IOP, y=ANGPT2)) + geom_boxplot(color="blue", fill="pink", size=1)
Anova_results <- aov(ANGPT2 ~ Class_IOP, data=corrTable)
summary(Anova_results)

############################################
### Step 5: Partial Correlation Analysis ###
############################################

corrTable <- cbind(sub_samples$Class_Age, sub_samples$Sex, sub_samples$Batch, corrTable)
names(corrTable)[1:3] <- c("Age", "Sex", "Batch")
levels(corrTable$Age)
levels(corrTable$Sex)
levels(corrTable$Batch)

######## "Sex" based scatter plot ######## 

ggscatter(corrTable, x="IOP", y=c("ANGPT2", "PTPRB", "TEK"), size=3, shape=19, color="blue", cor.method="spearman", 
          title="Correlation: Spearman,    Normalization: rlog", combine = TRUE, add="reg.line", conf.int=TRUE, 
          cor.coef=TRUE, xlab="IOP", ylab="Expression")
ggscatter(corrTable, x="IOP", y=c("ANGPT2", "PTPRB", "TEK"), size=3, shape=19, color="Sex", cor.method="spearman", 
          title="Correlation: Spearman,    Normalization: rlog", combine = TRUE, add="reg.line", conf.int=FALSE, 
          cor.coef=FALSE, xlab="IOP", ylab="Expression", palette=c("red", "pink", "blue"))
ggscatter(corrTable, x="IOP", y=c("ANGPT2", "PTPRB", "TEK"), size=3, shape=19, color="Sex", cor.method="spearman", 
          title="Correlation: Spearman,    Normalization: rlog", combine=TRUE, add="reg.line", conf.int=FALSE, 
          cor.coef=FALSE, xlab="IOP", ylab="Expression") + stat_cor(aes(color=Sex), label.x=3)
ggscatter(corrTable, x="ANGPT2", y=c("PTPRB", "TEK"), size=3, shape=19, color="Sex", cor.method="spearman", 
          title="Correlation: Spearman,    Normalization: rlog", combine=TRUE, add="reg.line", conf.int=FALSE, 
          cor.coef=FALSE, xlab="ANGPT2", ylab="Expression") + stat_cor(aes(color=Sex))
# ggscatter(corrTable_new, x="PTPRB", y="TEK", size=3, shape=19, color="Sex", cor.method="spearman",
#           title="Correlation: Spearman,    Normalization: rlog", combine = TRUE, add="reg.line", conf.int=FALSE, 
#           cor.coef=FALSE, xlab="PTPRB", ylab="TEK", palette=c("red", "pink", "blue"))

######## "Age" based scatter plot ########

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

######## "Batch" based scatter plot ########

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
           
######## Partial Correlation Tests(Dr. Chen) ######## 
corrTable <- cbind(sub_samples$Avg_IOP, sub_samples$AgeInDays, sub_samples$Batch, selGenes)    # Correlation tables for three genes 
names(corrTable) <- c("IOP", "Age", "Batch", "ANGPT2", "PTPRB", "TEK")

pcor.test(corrTable$ANGPT2, corrTable$IOP, corrTable$Age, method="pearson")
summary(lm(corrTable$ANGPT2~corrTable$IOP+corrTable$Age))

pcor.test(corrTable$ANGPT2, corrTable$IOP, corrTable$Batch, method="pearson")
summary(lm(corrTable$ANGPT2~corrTable$IOP+corrTable$Age+corrTable$Batch))

summary(lm(corrTable$ANGPT2~corrTable$IOP))









# write.table(ThreeGenes, file="ThreeGenes.csv", sep=",", quote=F, row.names=TRUE, col.names=TRUE,)
