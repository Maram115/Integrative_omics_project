
##################################### 
#load Libraries whyspecific lib?
#####################################

#install.packages("BiocManager")
#BiocManager::install("DESeq2")
#install_github("vqv/ggbiplot")
library(DESeq2)
library(ggplot2)
library(devtools) # installing and managing R packages hosted on GitHub.
install_github("vqv/ggbiplot")
library(ggbiplot)
library(ggsci) #provides a collection of themes and color palettes for ggplot2
library(RColorBrewer)#provides color palettes for creating visually appealing 
library(NMF)# Non-negative Matrix Factorization:used for analyzing and visualizing patterns in non-negative data matrices.
library(pheatmap) #for creating heatmaps
library(devtools)
library(ggfortify) #extends ggplot2 for visualizing statistical results from various modeling techniques.
library(gtools)#working with vectors, sorting data, and manipulating strings.


#####################################################
# Read Data
#####################################################

#loading the MIRNA Data as matrix to use it in deseq2...why?
mirna_data = as.matrix(read.table("mirna", quote="\"", comment.char="" , row.names=1))
#loading the melanoma Data (Metadata)
melanoma_sample <- read.delim("melanoma", header=T ,row.names = 1 )
#########################################################################################################################################################################
#show the dimension of the Data as how many columns that present samples and how many rows that present genes in this data
dim(mirna_data)
dim(melanoma_sample)
###########################################################################################################
########Quality control steps####
#check if is there any missing value of expression 
sum(is.na(mirna_data))
sum(is.na(melanoma_sample))

# Check for duplicates
any_duplicates <- any(duplicated(melanoma_sample))
any_duplicates_mirna <- any(duplicated(mirna_data))
####remove duplicates
#mirna_data_duplicaties_removing <- mirna_data[!duplicated(mirna_data), ]
#dim(mirna_data_duplicaties_removing )
#summary(mirna_data_duplicaties_removing )
############################################################################################################
# We interesting about vital status because it gives us the information about who dead and who is living !!!!!!!and i know why ? 
table(melanoma_sample$vital_status)
table(melanoma_sample$gender)
#remove NA from condition column (vital_status)
meta_data <- melanoma_sample[!is.na(melanoma_sample$vital_status), ] 

################################################################################################################33

#patients (samples) name are dotted in mirna data and dashed in meta data (melanoma) to make them matched in dashed 
colnames(mirna_data) <- gsub (".", "-", colnames(mirna_data), fixed = TRUE) 
colnames(mirna_data)


#intersect mirna data and meta data to make cols of mirna data equal rows of meta data to apply deseq2...why?
mirna_data1=mirna_data[,intersect(colnames(mirna_data),rownames(meta_data))]
#any_duplicates_mirna1<- any(duplicated(mirna_data1))
dim(mirna_data1)
dim(meta_data)

#plot histogram to visualize the distribution of the data 
hist(log2(mirna_data1+1), col = "orange", main="Histogram") 
# The histogram result shows that the data is right skewed 

# The scale the data by log2 transformation for better visualization, the +1 at the end of command is to avoid the infinity at log the values equal to zero 
boxplot(log(mirna_data1[1:5,]+1))

#draw QQ plot for checking the normality of data
qqnorm(mirna_data1[1,])
qqline(mirna_data1[1,])

# the qq plot results show that the data is not normally distributed as not all the points are on the line of qq norm 
#in deseq2, it is critical that the columns of the data and the rows of the metadata (melanoma) that gives information about samples, to be in the same order as DESeq2 will not know which column of the count matrix belongs to which row of the column of the data, these must be provided to deseq2 in a consistent order so the following code will do that

meta_data2=meta_data[colnames(mirna_data1),]
dim(mirna_data1)
dim(meta_data2)

#check the order to know if they are in a consistent order or not
all(colnames(mirna_data1) == rownames(meta_data2))
all(colnames(mirna_data1) %in% rownames(meta_data2))
identical(rownames(meta_data2), colnames(mirna_data1))
###################################################################################################

#The deseq2 package require the count data values to be integers 
#save the gene names in a variable
genes=row.names(mirna_data1)
#convert the data values to integers
mirna_data1=apply(round(mirna_data1),2,as.integer)
#rename the rows of the data
row.names(mirna_data1)=genes

#################### DO the differential EXP analysis using DeSeq2 ####################################
#creat a deseq dataset object
dds= DESeqDataSetFromMatrix( countData = mirna_data1 , colData = meta_data2, design = ~ vital_status )

# Run pipeline for differential expression steps
dds.run = DESeq(dds)
plotDispEsts(dds.run)
#specify how many conditions do you want to compare according to the phenotypic table
cond1="DECEASED"
cond2="LIVING"
str(meta_data2$vital_status)
#specifying the contrast (to make a res object based on two specific conditions)
res=results(dds.run, contrast = c("vital_status",cond1,cond2))

# remove nulls
res=as.data.frame(res[complete.cases(res), ])

#choose the statstical significant differentaily expressed transcripts of mirnas (DETs) based on the p adjusted value less than 0.05 and biological significance  based on the fold change more than 2
deseq.deg=res[res$padj < 0.05 & abs(res$log2FoldChange)> 0.5,]
#plotDispEsts(deseq.deg)
deseq.deg1=res[res$padj < 0.05 & abs(res$log2FoldChange)> 1.2,]
#export the Dets into your current folder for further analysthis
#write.csv(as.matrix(deseq.deg1),file="deseq.detf.csv", quote=F,row.names=T)


#draw DETs volcano plot for logfold>1.2
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(padj), pch=20, main="DECEASED vs LIVING"))
with(subset(res, padj<.05 & (log2FoldChange)> 1.2 ), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(res, padj<.05 & (log2FoldChange)< -1.2), points(log2FoldChange, -log10(padj), pch=20, col="red"))
legend(x=-0.5,y=4,c("Upregulated","Downgulated"), cex=.8, bty="n", col=c("blue","red"),pch=20)

#draw DETs volcano plot for logfold>0.5
par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(padj), pch=20, main="DECEASED vs LIVING"))
with(subset(res, padj<.05 & (log2FoldChange)> 0.5 ), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(res, padj<.05 & (log2FoldChange)< -0.5), points(log2FoldChange, -log10(padj), pch=20, col="red"))
legend(x=-2,y=10,c("Upregulated","Downgulated"), cex=.8, bty="n", col=c("blue","red"),pch=20)

####draw heatmap####
#normalize the data
dds2 <- estimateSizeFactors(dds)
normalized_counts <- as.data.frame(counts(dds2, normalized=TRUE))
write.csv(as.matrix(normalized_counts),file="count_table.csv", quote=F,row.names=T)
#plotDispEsts(normalized_counts)
#extract counts values of DETs(diffrentially expressed transcripts(micrornas)) only for each stage for logfold>1.2
exp.degs=as.matrix(normalized_counts[rownames(normalized_counts) %in% rownames(deseq.deg), ])
aheatmap(log2(exp.degs+1), annCol =meta_data2$vital_status, col = rev(brewer.pal(9,"RdBu")))

#extract counts values of DETs(diffrentially expressed transcripts(micrornas)) only for each stage for logfold>0.5
exp.degs=as.matrix(normalized_counts[rownames(normalized_counts) %in% rownames(deseq.deg1), ])
aheatmap(log2(exp.degs+1), annCol =meta_data2$vital_status, col = rev(brewer.pal(9,"RdBu")))


############################################################################################################################
#PCA 
#transpose rows with columns 
mirna_data2 <- t(mirna_data1)
dim(mirna_data2)
dim(mirna_data1)
PCA <- prcomp(mirna_data2)
plot(PCA$x[,1], PCA$x[,2])
plot(PCA)
plot(PCA ,type="l")
biplot(PCA , scale = 0)
# Extract PC scores
pcaData <- as.data.frame(PCA$x[, 1:2])
pcaData <- cbind(pcaData, meta_data2$vital_status) 
colnames(pcaData) <- c("PC1", "PC2", "Species")
# plot with ggplot 
library(ggplot2)
ggplot(pcaData) +
  aes(PC1, PC2, color = Species, shape = Species) + # define plot area
  geom_point(size = 2) 
#adding data points
percentVar <- round(100 * summary(PCA)$importance[2, 1:2], 0)
ggplot(pcaData, aes(PC1, PC2, color = Species, shape = Species)) + # starting ggplot2
  geom_point(size = 2) + # add data points
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + # x label
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + # y label
  ggtitle("Principal component analysis (PCA)") + # title
  theme(aspect.ratio = 1) # width and height ratio

ggplot(pcaData, aes(PC1, PC2, col = Species, fill = Species)) +
  stat_ellipse(geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, col = "black")
#####################################################################
#PCA (gender)
#transpose rows with columns 
mirna_data2 <- t(mirna_data1)
dim(mirna_data2)
dim(mirna_data1)
PCA <- prcomp(mirna_data2)
plot(PCA$x[,1], PCA$x[,2])
plot(PCA)
plot(PCA ,type="l")
biplot(PCA , scale = 0)
# Extract PC scores
pcaData <- as.data.frame(PCA$x[, 1:2])
pcaData <- cbind(pcaData, meta_data2$gender) 
colnames(pcaData) <- c("PC1", "PC2", "Gender")
# plot with ggplot 
library(ggplot2)
ggplot(pcaData) +
  aes(PC1, PC2, color = Gender, shape = Gender) + # define plot area
  geom_point(size = 2) 
#adding data points
percentVar <- round(100 * summary(PCA)$importance[2, 1:2], 0)
ggplot(pcaData, aes(PC1, PC2, color = Gender, shape = Gender)) + # starting ggplot2
  geom_point(size = 2) + # add data points
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + # x label
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + # y label
  ggtitle("Principal component analysis (PCA)") + # title
  theme(aspect.ratio = 1) # width and height ratio

ggplot(pcaData, aes(PC1, PC2, col = Gender, fill = Gender)) +
  stat_ellipse(geom = "polygon", col = "black", alpha = 0.5) +
  geom_point(shape = 21, col = "black")

