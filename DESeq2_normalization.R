#########################################################################################
### Step 1: Transfer gene counts from RSEM outputs to a .csv file, for all 53 sampels ###
#########################################################################################

setwd("/home/h/Desktop/RSEM/export"); getwd()
TopGenes = 32883
GeneResultsFiles <- list.files(pattern = "*.genes.results")
for (sample in 1:length(GeneResultsFiles)){
  sample_name <- substr(GeneResultsFiles[sample], 1, 12)
  cat ("\n---------------", sample_name , "-----------------------\n")
  data = read.table(GeneResultsFiles[sample], header = T, stringsAsFactors = F)
  idx = order(data[,"gene_id"], decreasing = F)
  # print(data[idx[1:TopGenes], c("gene_id", "TPM")])
  if (sample == 1) all_data = data[idx[1:TopGenes], c("gene_id", "TPM")]
  all_data[, sample+1] = data[idx[1:TopGenes], c("TPM")]
  names(all_data)[sample+1] <- sample_name}
setwd("/home/h/Desktop/RSEM/R_files"); getwd()
raw = write.table(all_data, file = "raw_gene_counts.csv", sep = ",", row.names = F)

###################################################################
### Step 2: Converting raw counts file to the DESeq2 data frame ###
###################################################################

# RawCountTable <- read.csv("raw_gene_counts.csv", header = TRUE, sep = ",", row.names = 1)
RawCountTable <- read.csv("raw_gene_counts_abstract.csv", header = TRUE, sep = ",", row.names = 1)
head(RawCountTable)
RoundCountTable = round(RawCountTable[,1:5],0)
class(RoundCountTable)
head(RoundCountTable)

library("DESeq2")
DataStructure <- data.frame(row.names = colnames(RoundCountTable), 
                      sample_name = colnames(RoundCountTable),
                      library_type = factor(c(rep("two_reads", 2), rep("four_reads", 3))))
DataStructure
colnames(RoundCountTable)
rownames(DataStructure)
dds <- DESeqDataSetFromMatrix(countData = RoundCountTable, 
                              colData = DataStructure, 
                              design = ~library_type)   # sample_name, library_type 
dds
colData(dds)     # Green box
rowData(dds)     # Blue box
assay(dds)       # Pink box

#########################################
### Step 3: Pre-filtering the dataset ###
#########################################

nrow(dds)
keep <- rowSums(counts(dds)) > 0       # remove genes without expression in any sample
# keep <- rowSums(counts(dds)) > 1       # remove genes without expression in more than one sample
dds <- dds[keep,]
nrow(dds)

##############################  ####################################
### Step 4: The Variance stabilizing transformation using rlog ###
##################################################################
