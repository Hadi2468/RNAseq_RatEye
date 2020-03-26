## Read in the data and metadata
df0 <- read.csv(file="./Dr_Chen.csv", header=T)
meta<-read.csv(file="./p50_iop_all_data.csv", header=T)
head(df0)
head(meta)

## remove IOP missing
df0 <- subset(df0, 0 < IOP)

## harmonize the names
df0$sample <- gsub("s_", "", df0$sample)

## merge data with metadata
df1 <- merge(df0, meta, by.x="sample", by.y="RatID")

## calculate mean iop
df1$iop <- apply(df1[,c("OD1", "OD2", "OD3", "OS1", "OS2", "OS3")] , 1, mean)

pdf(file="IOP_angpt2_correlation_hadi.pdf", width=11, height=8.5)
par(mfcol=c(1,3))

## check IOP distribution
plot(density(df1$IOP), main="Density plot for IOP")
plot(density(df1$iop), main="Density plot for iop")

## adjust IOP by age
df1$iopage <- df1$iop / df1$AgeInDays

## check density
plot(density(df1$iopage), main="Density plot for age-ajusted IOP")

cor.test(df1$iopage, df1$ANGPT2_rlog, method="spearman")
plot(df1$iopage, df1$ANGPT2_rlog, main="correlation")
mtext(side=3, "rho=0.29, p=0.057")
abline(lm(df1$ANGPT2_rlog ~ df1$iopage))
dev.off()

