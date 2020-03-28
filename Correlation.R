##############################################################################################
### DESeq2 Analysis and finding the corrolation between transcriptomic and IOP for 45 Rats ###
##############################################################################################

# BiocManager::install("")
# install.packages("corrplot")
# install.packages("RColorBrewer")
# install.packages("ggpubr")
# install.packages("pheatmap")
# install.packages("ppcor")
# install.packages("rcompanion")
# install.packages("caret")
# install.packages("kader")
# install.packages("moments")
# install.packages("remotes")
# remotes::install_github("ProcessMiner/nlcor")
# remove.packages("nlcor")
# install.packages("readxl")
# install.packages('NNS')
# install.packages("mgcv")
# install.packages("gam")
# install.packages("nlme")

Packages <- c("tximport", "tximportData", "DESeq2", "tidyverse", "dplyr", "vctrs", "fs", "ggplot2", 
              "kader", "remotes", "nlcor", "corrplot", "RColorBrewer", "ggpubr", "pheatmap", "NNS", 
              "ppcor", "BBmisc", "rcompanion", "caret", "moments", "devtools", "readxl", "mgcv")
lapply(Packages, library, character.only=TRUE)

####################################################
### Step 1: Import 45 samples' information files ###
####################################################
getwd()
setwd("./Data"); getwd()
# samples <- read.csv("samples53.csv", sep=",", header=TRUE)
sub_samples <- read.csv("samples45.csv", sep=",", header=TRUE)
dim(sub_samples)

####### Looking for a best partion of Age ########
for (i in (1:45)) {
    if (113 <= sub_samples$AgeInDays[i] & sub_samples$AgeInDays[i] < 127) {sub_samples$Class_Age3[i] <- "Adolescent"}
    else if (127 <= sub_samples$AgeInDays[i] & sub_samples$AgeInDays[i] < 167) {sub_samples$Class_Age3[i] <- "Adult"}
    else {sub_samples$Class_Age3[i] <- "Aged"} }

corrTable <- cbind(sub_samples$Class_Age, sub_samples$Class_Age3, sub_samples$Sex, sub_samples$Batch, sub_samples$Avg_IOP, selGenes)
names(corrTable) <- c("Age", "Age3", "Sex", "Batch", "IOP", "ANGPT2", "PTPRB", "TEK")
ggscatter(corrTable, x="IOP", y=c("ANGPT2", "PTPRB", "TEK"), size=3, shape=19, color="Age3", 
          cor.method="spearman", title="Correlation: Spearman,    Normalization: rlog", combine=TRUE, 
          add="reg.line", conf.int=FALSE, cor.coef=FALSE, xlab="IOP", ylab="Expression") + 
  stat_cor(aes(color=Age3), label.x=5)
sort(sub_samples$AgeInDays)
ggscatter(corrTable, x="ANGPT2", y=c("PTPRB", "TEK"), size=3, shape=19, color="Age3", 
          cor.method="spearman", title="Correlation: Spearman,    Normalization: rlog", combine=TRUE, 
          add="reg.line", conf.int=FALSE, cor.coef=FALSE, xlab="ANGPT2", ylab="Expression") + 
  stat_cor(aes(color=Age3), label.x=6.3)
write.table(sub_samples, file="samples45.csv", sep=",", quote=F, row.names=TRUE, col.names=TRUE,)
####################################################################################################

# sub_genes_2 <- read.csv("normalized_log2.csv", sep=",", header=TRUE)
# dim(sub_genes_2)
# sub_genes_r <- read.csv("normalized_rlog.csv", sep=",", header=TRUE)
# sub_genes_r <- read.csv("normalized_rlog_IOP.csv", sep=",", header=TRUE)     # Genes + IOP
sub_genes_r <- read.csv("normalized_rlog_Michael.csv", sep=",", header=TRUE)   # Solution by Michael Love
dim(sub_genes_r)
# sub_genes_v <- read.csv("normalized_vst.csv", sep=",", header=TRUE)
# dim(sub_genes_v)
# sub_genes_c <- read.csv("real_counts.csv", sep=",", header=TRUE)
# dim(sub_genes_c)

###########################################################################
### Step 2: Correlation Analysis: Heatmap, Clustering, and Distribution ###
###########################################################################

# selGenes <- subset(sub_genes_2, select=c(ENSRNOG00000016696_Angpt2, ENSRNOG00000055293_Ptprb, ENSRNOG00000008587_Tek))
selGenes <- subset(sub_genes_r, select=c(ENSRNOG00000016696_Angpt2, ENSRNOG00000055293_Ptprb, ENSRNOG00000008587_Tek))
# names(sub_genes_r)[32884] <- "IOP_norm"
# selGenes <- subset(sub_genes_r, select=c(IOP_norm, ENSRNOG00000016696_Angpt2, ENSRNOG00000055293_Ptprb, ENSRNOG00000008587_Tek))
# selGenes <- subset(sub_genes_v, select=c(ENSRNOG00000016696_Angpt2, ENSRNOG00000055293_Ptprb, ENSRNOG00000008587_Tek))
# selGenes <- subset(sub_genes_c, select=c(ENSRNOG00000016696_Angpt2, ENSRNOG00000055293_Ptprb, ENSRNOG00000008587_Tek))

corrTable <- cbind(sub_samples$Avg_IOP, selGenes)    # Correlation tables for three genes 
names(corrTable) <- c("IOP", "ANGPT2", "PTPRB", "TEK")
# names(corrTable) <- c("IOP", "IOP_norm", "ANGPT2", "PTPRB", "TEK")
summary(corrTable)     # Basic statistical analysis

# Dr_Chen <- cbind(samples[,c(1, 22, 20, 21, 19, 8, 23, 4)], 
#                  subset(sub_genes_r, select=c(ENSRNOG00000016696_Angpt2, ENSRNOG00000055293_Ptprb, ENSRNOG00000008587_Tek)),
#                  subset(sub_genes_v, select=c(ENSRNOG00000016696_Angpt2, ENSRNOG00000055293_Ptprb, ENSRNOG00000008587_Tek)),
#                  subset(sub_genes_c, select=c(ENSRNOG00000016696_Angpt2, ENSRNOG00000055293_Ptprb, ENSRNOG00000008587_Tek)))
# names(Dr_Chen) <- c("sample", "IOP", "OD", "OS", "IOP_class", "Age", "Age_class", "Sex", "ANGPT2_rlog", "PTPRB_rlog", "TEK_rlog",
#                     "ANGPT2_vst", "PTPRB_vst", "TEK_vst", "ANGPT2_raw", "PTPRB_raw", "TEK_raw")
# write.table(Dr_Chen, file="Dr_Chen.csv", sep=",", quote=F, row.names=TRUE, col.names=TRUE,)

######## Check Normality ########
# ggqqplot(corrTable$ANGPT2)
# ggqqplot(corrTable$PTPRB)
# ggqqplot(corrTable$TEK)
# ggqqplot(corrTable$IOP)
# shapiro.test(corrTable$ANGPT2)
# shapiro.test(corrTable$PTPRB)
# shapiro.test(corrTable$TEK)
# shapiro.test(corrTable$IOP)

######## IOP Normalization Methods ########
##### Method 2: standard transformation #####
# corrTable$IOP <- (corrTable$IOP - mean(corrTable$IOP)) / sd(corrTable$IOP)
# corrTable$IOP <- normalize(corrTable$IOP, method="standardize", range=c(0, 1), margin=1L, on.constant="quiet")
##### Method 3: Box-Cox power transformation #####
# Lambda = preProcess(corrTable, method=c("BoxCox"))[2]$bc$IOP$lambda
# corrTable$IOP = (corrTable$IOP ^ Lambda - 1) / Lambda
##### Method 4: Log transformation #####
# corrTable$IOP <- log(corrTable$IOP)
# corrTable$IOP <- log2(corrTable$IOP)
# corrTable$IOP <- log10(corrTable$IOP)
##### Method 5: Square root transformation #####
# corrTable$IOP <- sqrt(corrTable$IOP)
##### Method 6: Cube root transformation #####
# corrTable$IOP <- kader:::cuberoot(corrTable$IOP)

# skewness(corrTable$IOP)

##### Method 7: Dr. Chen Adjusting #####
# df0 <- read.csv(file="../Dr.Chen/Dr_Chen.csv", header=T)
# meta<-read.csv(file="../Dr.Chen/p50_iop_all_data.csv", header=T)
# df0 <- subset(df0, 0 < IOP)
# df0$sample <- gsub("s_", "", df0$sample)
# df1 <- merge(df0, meta, by.x="sample", by.y="RatID")
# df1$iop <- apply(df1[,c("OD1", "OD2", "OD3", "OS1", "OS2", "OS3")] , 1, mean)
# df1$iopage <- df1$iop / df1$AgeInDays
# par(mfcol=c(1,2))
# gg1 <- ggplot(df1, aes(IOP)) + geom_density(fill="green") + ggtitle("Density plot for IOP")
# gg2 <- ggplot(df1, aes(iopage)) + geom_density(fill="green") + ggtitle("Density plot for age-ajusted IOP")
# ggarrange(gg1, gg2, ncol=2)
# shapiro.test(df1$iopage)
# skewness(df1$iopage)
# ## Correlation ##
# corrTable <- cbind(sub_samples$Avg_IOP, selGenes, df1$iopage)    # Correlation tables for three genes 
# names(corrTable) <- c("IOP", "ANGPT2", "PTPRB", "TEK", "IOP_Age")
# cor.test(df1$iopage, df1$ANGPT2_rlog, method="spearman")
# plot(df1$iopage, df1$ANGPT2_rlog, main="correlation") 
# mtext(side=3, "rho=0.29, p=0.057")
# abline(lm(df1$ANGPT2_rlog ~ df1$iopage))
# dev.off()

######## Pearson Correlatoin ########
# pheatmap(cor(corrTable))
# 
# corrplot(cor(corrTable, method="pearson"), method="color", type="upper", order="hclust", 
#          col=colorRampPalette(c("dodgerblue", "aliceblue", "brown1"))(7), 
#          addCoef.col="black", tl.col="black", tl.cex=1, addrect=3)

######## Spearman Correlatoin ########
plot(corrTable$IOP, corrTable$ANGPT2)
corrplot(cor(corrTable, method="spearman"), method="color", type="upper", order="hclust", 
         col=colorRampPalette(c("dodgerblue", "aliceblue", "brown1"))(7), 
         addCoef.col="black", tl.col="black", tl.cex=1, addrect=3)

#################### NonLinear Correlatoin (nlcor) ####################
nlcor(corrTable$IOP, corrTable$ANGPT2, refine = 0.7, plt=T)

nlcor(corrTable$IOP, corrTable$IOP)
nlcor(corrTable$IOP, corrTable$ANGPT2)
nlcor(corrTable$IOP, corrTable$PTPRB)
nlcor(corrTable$IOP, corrTable$TEK)

nlcor(corrTable$ANGPT2, corrTable$IOP)
nlcor(corrTable$ANGPT2, corrTable$ANGPT2)
nlcor(corrTable$ANGPT2, corrTable$PTPRB)
nlcor(corrTable$ANGPT2, corrTable$TEK)

nlcor(corrTable$PTPRB, corrTable$IOP)
plot(corrTable$PTPRB, corrTable$ANGPT2)
nlcor(corrTable$PTPRB, corrTable$ANGPT2, plt=T)       # R=0.3664747    p=0.006690839
nlcor(corrTable$PTPRB, corrTable$PTPRB)
nlcor(corrTable$PTPRB, corrTable$TEK)

nlcor(corrTable$TEK, corrTable$IOP, plt=T)            # R=0.4438115     p=0.01823332
nlcor(corrTable$TEK, corrTable$ANGPT2)
nlcor(corrTable$TEK, corrTable$PTPRB, plt=T)          # R=0.3986897     p=0.01190667
nlcor(corrTable$TEK, corrTable$TEK)

##################### NonLinear Correlatoin (NNS) ##################################################
## Nonlinear Nonparametric Statistics (NNS) using partial moments.                                ##
## Partial moments are the elements of variance and asymptotically approximate the area of f(x).  ##
## NNS offers Nonlinear Correlation & Dependence analysis.                                        ##
####################################################################################################
NNS.cor(corrTable$IOP, corrTable$ANGPT2)
NNS.cor(corrTable$IOP, corrTable$PTPRB)
plot(corrTable$IOP, corrTable$TEK)
NNS.cor(corrTable$IOP, corrTable$TEK)                  # R (IOP_TEK) = 0.8501074
NNS.cor(corrTable$ANGPT2, corrTable$PTPRB)
NNS.cor(corrTable$ANGPT2, corrTable$TEK)               # R (ANGPT2_TEK) = 0.5067027
plot(corrTable$ANGPT2, corrTable$TEK)
NNS.cor(corrTable$PTPRB, corrTable$TEK)

######## Probability Density Function ########
ThreeGenes <- read.csv("ThreeGenes.csv", sep=",", header=TRUE)
ThreeGenes$Expression[1:45] <- corrTable$ANGPT2
ThreeGenes$Expression[46:90] <- corrTable$PTPRB
ThreeGenes$Expression[91:135] <- corrTable$TEK
ggplot(ThreeGenes, aes(Expression, fill=Gene)) + geom_density(alpha=0.6) +   scale_fill_manual(values=c("red", "blue", "yellow"))
ggplot(corrTable, aes(IOP)) + geom_density(fill="green") + scale_x_continuous(limits=c(10, 25))
ggplot(corrTable, aes(IOP)) + geom_density(fill="green") 

#######################################################
### Step 3: IOP Correlation Analysis: Scatter Plots ###
#######################################################

corrTable <- cbind(sub_samples$Class_IOP, corrTable)    # Correlation tables for three genes and three groups of IOP 
names(corrTable)[1] <- "Class_IOP"
levels(corrTable$Class_IOP)
corrTable$Class_IOP <- relevel(corrTable$Class_IOP, "Normal")
levels(corrTable$Class_IOP)

corrTable_IOP2 <- corrTable    # Two subgroups of IOP
for (i in (1:45)) {if (corrTable_IOP2$Class_IOP[i] == "Elevated") {corrTable_IOP2$Class_IOP[i] <- "High"}}
corrTable_IOP2$Class_IOP <- factor(corrTable_IOP2$Class_IOP)    # Correlation tables for just two groups of IOP
levels(corrTable_IOP2$Class_IOP)

ggscatter(corrTable, x="IOP", y=c("ANGPT2", "PTPRB", "TEK"), size=3, shape=19, color="blue", 
          cor.method="spearman", title="Correlation: Spearman,    Normalization: rlog", combine=TRUE, 
          add="reg.line", conf.int=TRUE, cor.coef=TRUE, xlab="IOP", ylab="Expression")
ggscatter(corrTable, x="IOP", y=c("ANGPT2", "PTPRB", "TEK"), size=3, shape=19, color="Class_IOP", 
          cor.method="spearman", title="Correlation: Spearman,    Normalization: rlog", combine=TRUE, 
          add="reg.line", conf.int=TRUE, cor.coef=FALSE, xlab="IOP", ylab="Expression") + 
  stat_cor(aes(color=Class_IOP), label.x=12)
ggscatter(corrTable_IOP2, x="IOP", y=c("ANGPT2", "PTPRB", "TEK"), size=3, shape=19, color="Class_IOP", 
          cor.method="spearman", title="Correlation: Spearman,    Normalization: rlog", combine=TRUE, 
          add="reg.line", conf.int=TRUE, cor.coef=FALSE, xlab="IOP", ylab="Expression") + 
  stat_cor(aes(color=Class_IOP), label.x=-3) 

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

#### IOP_Age (Dr. Chen) #####
ggscatter(corrTable, x="IOP_Age", y=c("ANGPT2", "PTPRB", "TEK"), size=3, shape=19, color="Class_IOP", 
          cor.method="spearman", title="Correlation: Spearman,    Normalization: rlog", combine=TRUE, 
          add="reg.line", conf.int=TRUE, cor.coef=FALSE, xlab="IOP_Age", ylab="Expression") + 
  stat_cor(aes(color=Class_IOP), label.x=-0.15)
ggscatter(corrTable_IOP2, x="IOP_Age", y=c("ANGPT2", "PTPRB", "TEK"), size=3, shape=19, color="Class_IOP", 
          cor.method="spearman", title="Correlation: Spearman,    Normalization: rlog", combine=TRUE, 
          add="reg.line", conf.int=TRUE, cor.coef=FALSE, xlab="IOP_Age", ylab="Expression") + 
  stat_cor(aes(color=Class_IOP), label.x=0.05) 
ggANG <- ggplot(corrTable, aes(x=IOP_Age, y=ANGPT2, color=Class_IOP)) + scale_size_manual(values=c(2,2,2)) +
  geom_point(aes(size=Class_IOP)) + geom_smooth(method=lm, aes(fill=Class_IOP), se=FALSE, fullrange=TRUE) + 
  theme(legend.position="top") + labs(y="Expression")
ggPT <- ggplot(corrTable, aes(x=IOP_Age, y=PTPRB, color=Class_IOP)) + scale_size_manual(values=c(2,2,2)) +
  geom_point(aes(size=Class_IOP)) + geom_smooth(method=lm, aes(fill=Class_IOP), se=FALSE, fullrange=TRUE) + 
  theme(legend.position="top", axis.title.y = element_blank())
ggTEK <- ggplot(corrTable, aes(x=IOP_Age, y=TEK, color=Class_IOP)) + scale_size_manual(values=c(2,2,2)) +
  geom_point(aes(size=Class_IOP)) + geom_smooth(method=lm, aes(fill=Class_IOP), se=FALSE, fullrange=TRUE) + 
  theme(legend.position="top", axis.title.y = element_blank())
ggarrange(ggANG, ggPT, ggTEK + rremove("x.text"), labels=c("ANGPT2", "PTPRB", "TEK"), ncol=3, nrow=1)

#################################################
### Step 4: t-test and ANOVA Analysis for IOP ###
#################################################

############ t-test Analysis ############

# H0: Mean of ANGPT2 expression for samples with Normal_IOP = Mean of ANGPT2 expression for samples with High_IOP 
ggplot(corrTable_IOP2, aes(x=Class_IOP, y=ANGPT2)) + geom_boxplot(color="blue", fill="pink", size=1.2)
t.test(ANGPT2 ~ Class_IOP, data=corrTable_IOP2, mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)

############ ANOVA Analysis #############

# H0: Mean of ANGPT2 expression for samples with all IOP gropus are equal 
ggplot(corrTable, aes(x=Class_IOP, y=ANGPT2)) + geom_boxplot(color="blue", fill="pink", size=1)
Anova_results <- aov(ANGPT2 ~ Class_IOP, data=corrTable)
summary(Anova_results)

##############################################################
### Step 5: Partial Correlation Analysis for IOP subgroups ###
##############################################################

corrTable <- cbind(sub_samples$Avg_IOP, selGenes)    # Correlation tables for three genes 
names(corrTable) <- c("IOP", "ANGPT2", "PTPRB", "TEK")

corrTable <- cbind(sub_samples$Class_Age, sub_samples$Class_Age3, sub_samples$Sex, sub_samples$Batch, corrTable)
names(corrTable)[1:3] <- c("Age", "Age3", "Sex", "Batch")
levels(corrTable$Age)
levels(corrTable$Age3)
levels(corrTable$Sex)
levels(corrTable$Batch)

############################################
### Step 6: Partial Correlation Analysis ###
###            "Age" based               ###
############################################

ggscatter(corrTable, x="IOP", y=c("ANGPT2", "PTPRB", "TEK"), size=3, shape=19, color="Age", 
          cor.method="spearman", title="Correlation: Spearman,    Normalization: rlog", combine=TRUE, 
          add="reg.line", conf.int=FALSE, cor.coef=FALSE, xlab="IOP", ylab="Expression") + 
  stat_cor(aes(color=Age), label.x=12)
ggscatter(corrTable, x="ANGPT2", y=c("PTPRB", "TEK"), size=3, shape=19, color="Age", 
          cor.method="spearman", title="Correlation: Spearman,    Normalization: rlog", combine=TRUE, 
          add="reg.line", conf.int=FALSE, cor.coef=FALSE, xlab="ANGPT2", ylab="Expression") + 
  stat_cor(aes(color=Age), label.x=6.3)
ggscatter(corrTable, x="IOP", y=c("ANGPT2", "PTPRB", "TEK"), size=3, shape=19, color="Age3", 
          cor.method="spearman", title="Correlation: Spearman,    Normalization: rlog", combine=TRUE, 
          add="reg.line", conf.int=FALSE, cor.coef=FALSE, xlab="IOP", ylab="Expression") + 
  stat_cor(aes(color=Age3), label.x=12)
ggscatter(corrTable, x="ANGPT2", y=c("PTPRB", "TEK"), size=3, shape=19, color="Age3", 
          cor.method="spearman", title="Correlation: Spearman,    Normalization: rlog", combine=TRUE, 
          add="reg.line", conf.int=FALSE, cor.coef=FALSE, xlab="ANGPT2", ylab="Expression") + 
  stat_cor(aes(color=Age3), label.x=6.3)

############ ANOVA Analysis #############

# H0: Mean of ANGPT2 expression for samples with all Age gropus are equal 
ggplot(corrTable, aes(x=Age, y=ANGPT2)) + geom_boxplot(color="blue", fill="pink", size=1)
Anova_results <- aov(ANGPT2 ~ Age, data=corrTable)
summary(Anova_results)

Anova_results <- aov(ANGPT2 ~ IOP + Age, data=corrTable)
summary(Anova_results)

######## Division 1  ########
corrTable_Age2 <- corrTable
levels(corrTable_Age2$Age) <- c("Adolescent", "Adolescent_&_Adult", "Middle_Aged", "Middle_Aged_&_Aged")
for (i in (1:45)) {if (corrTable_Age2$Age[i] == "Adolescent") {corrTable_Age2$Age[i] <- "Adolescent_&_Adult"}}
for (i in (1:45)) {if (corrTable_Age2$Age[i] == "Middle_Aged") {corrTable_Age2$Age[i] <- "Middle_Aged_&_Aged"}} 
corrTable_Age2$Age <- factor(corrTable_Age2$Age)    # Correlation tables for just two groups of Age
levels(corrTable_Age2$Age)

ggplot(corrTable_Age2, aes(x=Age, y=ANGPT2)) + geom_boxplot(color="blue", fill="pink", size=1)
Anova_results <- aov(ANGPT2 ~ Age, corrTable_Age2); summary(Anova_results)
Anova_results <- aov(ANGPT2 ~ IOP + Age, corrTable_Age2); summary(Anova_results)
t.test(ANGPT2 ~ Age, corrTable_Age2, mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)

######## Division 2 ########
corrTable_Age2 <- corrTable
levels(corrTable_Age2$Age) <- c("Others", "Adult", "Middle_Aged", "Aged")
for (i in (1:45)) {if (corrTable_Age2$Age[i] == "Middle_Aged") {corrTable_Age2$Age[i] <- "Others"}} 
for (i in (1:45)) {if (corrTable_Age2$Age[i] == "Aged") {corrTable_Age2$Age[i] <- "Others"}}
corrTable_Age2$Age <- factor(corrTable_Age2$Age)    # Correlation tables for just two groups of Age
levels(corrTable_Age2$Age)

ggplot(corrTable_Age2, aes(x=Age, y=ANGPT2)) + geom_boxplot(color="blue", fill="pink", size=1)
Anova_results <- aov(ANGPT2 ~ Age, corrTable_Age2); summary(Anova_results)
Anova_results <- aov(ANGPT2 ~ IOP + Age, corrTable_Age2); summary(Anova_results)
t.test(ANGPT2 ~ Age, corrTable_Age2, mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)

######## Division 3 ########
corrTable_Age2 <- corrTable
levels(corrTable_Age2$Age) <- c("Adolescent", "Others", "Middle_Aged", "Aged")
for (i in (1:45)) {if (corrTable_Age2$Age[i] == "Adolescent") {corrTable_Age2$Age[i] <- "Others"}}
for (i in (1:45)) {if (corrTable_Age2$Age[i] == "Middle_Aged") {corrTable_Age2$Age[i] <- "Others"}} 
corrTable_Age2$Age <- factor(corrTable_Age2$Age)    # Correlation tables for just two groups of Age
levels(corrTable_Age2$Age)

ggplot(corrTable_Age2, aes(x=Age, y=ANGPT2)) + geom_boxplot(color="blue", fill="pink", size=1)
Anova_results <- aov(ANGPT2 ~ Age, corrTable_Age2); summary(Anova_results)
Anova_results <- aov(ANGPT2 ~ IOP + Age, corrTable_Age2); summary(Anova_results)
t.test(ANGPT2 ~ Age, corrTable_Age2, mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)

######## Division 4 ########
corrTable_Age2 <- corrTable
levels(corrTable_Age2$Age) <- c("Adolescent", "Others", "Middle_Aged", "Aged")
for (i in (1:45)) {if (corrTable_Age2$Age[i] == "Middle_Aged") {corrTable_Age2$Age[i] <- "Others"}} 
for (i in (1:45)) {if (corrTable_Age2$Age[i] == "Aged") {corrTable_Age2$Age[i] <- "Others"}}
corrTable_Age2$Age <- factor(corrTable_Age2$Age)    # Correlation tables for just two groups of Age
levels(corrTable_Age2$Age)

ggplot(corrTable_Age2, aes(x=Age, y=ANGPT2)) + geom_boxplot(color="blue", fill="pink", size=1)
Anova_results <- aov(ANGPT2 ~ Age, corrTable_Age2); summary(Anova_results)
Anova_results <- aov(ANGPT2 ~ IOP + Age, corrTable_Age2); summary(Anova_results)
t.test(ANGPT2 ~ Age, corrTable_Age2, mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)

############################################
### Step 7: Partial Correlation Analysis ###
###            "Sex" based               ###
############################################

summary(corrTable$Sex)
ggscatter(corrTable, x="IOP", y=c("ANGPT2", "PTPRB", "TEK"), size=3, shape=19, color="Sex", 
          cor.method="spearman", title="Correlation: Spearman,    Normalization: rlog", combine=TRUE, 
          add="reg.line", conf.int=FALSE, cor.coef=FALSE, xlab="IOP", ylab="Expression") + 
  stat_cor(aes(color=Sex), label.x=3)
ggscatter(corrTable, x="ANGPT2", y=c("PTPRB", "TEK"), size=3, shape=19, color="Sex", 
          cor.method="spearman", title="Correlation: Spearman,    Normalization: rlog", combine=TRUE, 
          add="reg.line", conf.int=FALSE, cor.coef=FALSE, xlab="ANGPT2", ylab="Expression") + 
  stat_cor(aes(color=Sex), label.x=6)
# ggscatter(corrTable_new, x="PTPRB", y="TEK", size=3, shape=19, color="Sex", cor.method="spearman",
#           title="Correlation: Spearman,    Normalization: rlog", combine = TRUE, add="reg.line", conf.int=FALSE, 
#           cor.coef=FALSE, xlab="PTPRB", ylab="TEK", palette=c("red", "pink", "blue"))

############ ANOVA Analysis #############

ggplot(corrTable, aes(x=Sex, y=ANGPT2)) + geom_boxplot(color="blue", fill="pink", size=1)
Anova_results <- aov(ANGPT2 ~ Sex, corrTable); summary(Anova_results)
Anova_results <- aov(ANGPT2 ~ IOP + Sex, corrTable); summary(Anova_results)
Linear_model <- lm(ANGPT2 ~ IOP + Sex, corrTable); summary(Linear_model)

######## Division 1  ########
corrTable_Sex2 <- corrTable
levels(corrTable_Sex2$Sex) <- c("F", "FP", "M")
for (i in (1:45)) {if (corrTable_Sex2$Sex[i] == "FP") {corrTable_Sex2$Sex[i] <- "F"}}
corrTable_Sex2$Sex <- factor(corrTable_Sex2$Sex)    # Correlation tables for just two groups of Sex
levels(corrTable_Sex2$Sex)
summary(corrTable_Sex2$Sex)

ggplot(corrTable_Sex2, aes(x=Sex, y=ANGPT2)) + geom_boxplot(color="blue", fill="pink", size=1)
Anova_results <- aov(ANGPT2 ~ Sex, corrTable_Sex2); summary(Anova_results)
Anova_results <- aov(ANGPT2 ~ IOP + Sex, corrTable_Sex2); summary(Anova_results)
Anova_results <- aov(ANGPT2 ~ IOP + Sex + Age, corrTable_Sex2); summary(Anova_results)
t.test(ANGPT2 ~ Sex, corrTable_Sex2, mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)

######## Division 2  ########
corrTable_Sex2 <- corrTable
levels(corrTable_Sex2$Sex) <- c("NonP", "P", "M")
for (i in (1:45)) {if (corrTable_Sex2$Sex[i] == "M") {corrTable_Sex2$Sex[i] <- "NonP"}}
corrTable_Sex2$Sex <- factor(corrTable_Sex2$Sex)    # Correlation tables for just two groups of Sex
levels(corrTable_Sex2$Sex)
summary(corrTable_Sex2$Sex)

ggplot(corrTable_Sex2, aes(x=Sex, y=ANGPT2)) + geom_boxplot(color="blue", fill="pink", size=1)
Anova_results <- aov(ANGPT2 ~ Sex, corrTable_Sex2); summary(Anova_results)
Anova_results <- aov(ANGPT2 ~ IOP + Sex, corrTable_Sex2); summary(Anova_results)
Anova_results <- aov(ANGPT2 ~ IOP + Sex + Age, corrTable_Sex2); summary(Anova_results)
t.test(ANGPT2 ~ Sex, corrTable_Sex2, mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)


############################################
### Step 8: Partial Correlation Analysis ###
###            "Batch" based             ###
############################################

summary(corrTable$Batch)
ggscatter(corrTable, x="IOP", y=c("ANGPT2", "PTPRB", "TEK"), size=3, shape=19, color="Batch", 
          cor.method="spearman", title="Correlation: Spearman,    Normalization: rlog", combine=TRUE, 
          add="reg.line", conf.int=FALSE, cor.coef=FALSE, xlab="IOP", ylab="Expression") + 
  stat_cor(aes(color=Batch), label.x=3)
ggscatter(corrTable, x="ANGPT2", y=c("PTPRB", "TEK"), size=3, shape=19, color="Batch", 
          cor.method="spearman", title="Correlation: Spearman,    Normalization: rlog", combine=TRUE, 
          add="reg.line", conf.int=FALSE, cor.coef=FALSE, xlab="ANGPT2", ylab="Expression") + 
  stat_cor(aes(color=Batch), label.x=6)

############ ANOVA Analysis #############

ggplot(corrTable, aes(x=Batch, y=ANGPT2)) + geom_boxplot(color="blue", fill="pink", size=1)
Anova_results <- aov(ANGPT2 ~ Batch, corrTable); summary(Anova_results)
Anova_results <- aov(ANGPT2 ~ IOP + Batch, corrTable); summary(Anova_results)
Anova_results <- aov(ANGPT2 ~ IOP + Batch + Sex + Age, corrTable); summary(Anova_results)
Linear_model <- lm(ANGPT2 ~ IOP + Batch, corrTable); summary(Linear_model)

######## Division 1  ########
corrTable_Batch2 <- corrTable
levels(corrTable_Batch2$Batch) <- c("B_1", "B13", "B_2", "B15")
for (i in (1:45)) {if (corrTable_Batch2$Batch[i] == "B13") {corrTable_Batch2$Batch[i] <- "B_1"}}
for (i in (1:45)) {if (corrTable_Batch2$Batch[i] == "B15") {corrTable_Batch2$Batch[i] <- "B_2"}}
corrTable_Batch2$Batch <- factor(corrTable_Batch2$Batch)    # Correlation tables for just two groups of Sex
levels(corrTable_Batch2$Batch)
summary(corrTable_Batch2$Batch)

ggplot(corrTable_Batch2, aes(x=Batch, y=ANGPT2)) + geom_boxplot(color="blue", fill="pink", size=1)
Anova_results <- aov(ANGPT2 ~ Batch, corrTable_Batch2); summary(Anova_results)
Anova_results <- aov(ANGPT2 ~ IOP + Batch, corrTable_Batch2); summary(Anova_results)
Anova_results <- aov(ANGPT2 ~ IOP + Batch + Age + Sex, corrTable_Batch2); summary(Anova_results)
t.test(ANGPT2 ~ Batch, corrTable_Batch2, mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)

######## Division 2  ########
corrTable_Batch2 <- corrTable
levels(corrTable_Batch2$Batch) <- c("B12", "B13", "B14", "Others")
for (i in (1:45)) {if (corrTable_Batch2$Batch[i] == "B12") {corrTable_Batch2$Batch[i] <- "Others"}}
for (i in (1:45)) {if (corrTable_Batch2$Batch[i] == "B13") {corrTable_Batch2$Batch[i] <- "Others"}}
corrTable_Batch2$Batch <- factor(corrTable_Batch2$Batch)    # Correlation tables for just two groups of Sex
levels(corrTable_Batch2$Batch)
summary(corrTable_Batch2$Batch)

ggplot(corrTable_Batch2, aes(x=Batch, y=ANGPT2)) + geom_boxplot(color="blue", fill="pink", size=1)
Anova_results <- aov(ANGPT2 ~ Batch, corrTable_Batch2); summary(Anova_results)
Anova_results <- aov(ANGPT2 ~ IOP + Batch, corrTable_Batch2); summary(Anova_results)
Anova_results <- aov(ANGPT2 ~ IOP + Batch + Age + Sex, corrTable_Batch2); summary(Anova_results)
t.test(ANGPT2 ~ Batch, corrTable_Batch2, mu=0, alt="two.sided", conf=0.95, var.eq=F, paired=F)

#######################################################
### Step 9: Partial Correlation Analysis (Dr. Chen) ###
#######################################################

corrTable <- cbind(sub_samples$Avg_IOP, sub_samples$AgeInDays, sub_samples$Sex, selGenes)    # Correlation tables for three genes 
names(corrTable) <- c("IOP", "Age", "Sex", "ANGPT2", "PTPRB", "TEK")

######## Partial Correlation (Dr. Chen) ##########
pcor.test(corrTable$ANGPT2, corrTable$IOP, corrTable$Age, method="spearman")    # p-value: 0.091
pcor.test(corrTable$ANGPT2, corrTable$IOP, corrTable$Sex, method="spearman")

######## Regression (Dr. Chen) ##########
fit1 <- lm(ANGPT2 ~ IOP, corrTable); summary(fit1)                         # p-value: 0.132
fit2 <- lm(ANGPT2 ~ IOP + Age, corrTable); summary(fit2)                   # p-value: 0.2225
fit3 <- lm(ANGPT2 ~ IOP + Sex, corrTable); summary(fit3)                   # p-value: 0.3689
fit4 <- lm(ANGPT2 ~ IOP + Age + Sex, corrTable); summary(fit4)             # p-value: 0.4302

########### Coefficients with labled data ###############
fit4 <- lm(ANGPT2 ~ IOP + Age + Sex, corrTable); summary(fit4)
#                    Estimate  std. Error t value  Pr(>|t|)
# (Intercept)        6.922953   0.119444  57.960   <2e-16
# IOP                0.010553   0.006894   1.531    0.134
# AgeA2_Adult       -0.050848   0.076139  -0.668    0.508
# AgeA3_Middle_Aged  0.043947   0.073989   0.594    0.556
# AgeA4_Aged        -0.085289   0.071434  -1.194    0.240
# SexFP             -0.083323   0.072568  -1.148    0.258
# SexM              -0.083855   0.063257  -1.326    0.193

######## ANOVA ##########
anova(fit1, fit2, fit3, fit4)

anova(lm(ANGPT2 ~ IOP, corrTable))                           # p-value: 0.132
anova(lm(ANGPT2 ~ IOP + Age, corrTable))                     # p-value: 0.1332                   
anova(lm(ANGPT2 ~ IOP + Sex, corrTable))                     # p-value: 0.1371
anova(lm(ANGPT2 ~ IOP + Age + Sex, corrTable))               # p-value: 0.1387

######## Omitting Intercept ##########
summary(lm(ANGPT2 ~ IOP - 1, corrTable))                     # p-value: < 2.2e-16
summary(lm(ANGPT2 ~ IOP + Age - 1, corrTable))               # p-value: < 2.2e-16
summary(lm(ANGPT2 ~ IOP + Sex - 1, corrTable))               # p-value: < 2.2e-16
summary(lm(ANGPT2 ~ IOP + Age + Sex - 1, corrTable))         # p-value: < 2.2e-16



