#Skin Cutaneous Melanoma
#BiocManager::install("DESeq2", force = TRUE)#one time
library(BiocManager)
library(DESeq2)#each time you open your R
library(ggplot2)
library(RColorBrewer)
library(NMF)

### read all files##############
count=as.matrix(read.table('exp', sep = " ", row.names = 1 , header = T)) 
methyl=read.delim('methy', sep =" ",row.names = 1 , header = T )
mirna=read.delim('mirna', sep =" ", row.names = 1 , header = T)
survival=read.delim('survival', sep ='\t', row.names = 1)
metadata=read.table('melanoma', sep = '\t', header = T, row.names = 1)

#####################################################################################################################
dim(count)  #it contains 20531 genes (rows)  && 473 samples (columns)
dim(metadata) #it contains 481 samples (rows) && 102 columns
### number of samples are not the same in count and metadata file
all(colnames(count) %in% rownames(metadata)) #FALSE

## i need to remove samples from count file that not present in metadata
# Step 1: Check which samples are present in the metadata
#samples_to_keep <- rownames(metadata)
# Step 2: Subset the count data to keep only the samples present in the metadata
#count_filtered <- count[, colnames(count) %in% samples_to_keep]

#all(colnames(count_filtered) %in% rownames(metadata))
# Now, 'count_filtered' contains only columns (samples) present in the 'metadata'

#check for NAs 
sum(is.na(count))
sum(is.na(metadata))

# I will use vital_status column as condition while creating Deseq object
table(metadata$vital_status)

#remove NA from vital_status column
meta_data_new <- metadata[!is.na(metadata$vital_status), ] 

#patients (samples) name are dotted in count file and dashed in metadata file so make the same dashed in both
colnames(count) <- gsub (".", "-", colnames(count), fixed = TRUE) 
colnames(count)

### NOW number of samples are the same in count and metadata file
all(colnames(count) %in% rownames(meta_data_new)) #TRUE

#############################################################################################
#normality check
#explore the data distribution using the histogram plot
hist(count, col = "orange", main="Histogram expression counts")

#scaling the data using log2 transformation to better visulization# we use (+1) to avoid the infinity character when we log zero valus 
hist(log2(count+1), col = "orange", main="Histogram of expression counts") #negative binomial distribution 
#zero inflated negative binomial distribution data

#draw QQ plot for checking the normality of data
qqnorm(log2(count[1000,]+1))
qqline(log2(count[1000,]+1))
# the qq plot results show that the data are NOT normally distributed

shapiro.test(count[20531,])
#p-value < 2.2e-16 means data are not normally distributed 

#variance check 
boxplot(log2(count[,1:50]+1)) # for samples to check normality
boxplot(log2(count[1:3,]+1)) # for first 3 genes

#explore if is there any missing expression value (empty cell)
sum(is.na(count))    #no missing values 
is.na(count)
is.null(count)
is.nan(count)

#It is absolutely critical that the columns of the count and the rows of 
#the meta_data_new (information about samples) are in the same order.DESeq2 will
#not make guesses as to which column of the count matrix belongsto which row
#of the column data, these must be provided to DESeq2 already in consistent order
meta_data_new2=meta_data_new[colnames(count),]
dim(count)
dim(meta_data_new2)
dim(meta_data_new)

all(colnames(count) == rownames(meta_data_new2)) #TRUE
identical(rownames(meta_data_new2), colnames(count)) #TRUE

#The deseq2 package require the count data values to be integers 
#save the gene names in a variable
genes=row.names(count)
#convert the data values to integers
count=apply(round(count),2,as.integer)
#view the data
head(count)
#rename the rows of the data
row.names(count)=genes
#view the data
#######################################################################################################
###### DO the differential EXP analysis using DeSeq2

#specify how many conditions do you want to compare according to metadata table
cond1="DECEASED"
cond2="LIVING"
str(meta_data_new2$vital_status) #Character
#convert vitual status column to factor
meta_data_new2$vital_status <- factor(meta_data_new2$vital_status, levels = c("LIVING", "DECEASED"))
#create a deseq dataset object as a matrix 
dds= DESeqDataSetFromMatrix( countData = count , colData = meta_data_new2, design = ~ vital_status)

#show design 
model.matrix(~0+meta_data_new2$vital_status)

#run the deseq2 workflow
dds.run = DESeq(dds)   #refit outliers
#specifying the contrast (to make a res object based on two specific conditions)
res=results(dds.run, contrast = c("vital_status",cond1 ,cond2))

####normalization after DESeq2
#dds2 <- estimateSizeFactors(dds)
#normalized_counts <- as.data.frame(counts(dds2, normalized=TRUE))
#write.csv(as.matrix(normalized_counts),file="normalized_counts_mRNA.csv", quote=F,row.names=T)

# remove nulls
res=as.data.frame(res[complete.cases(res), ])
###########################################################################################
####################visualization ###################################
#chose the statstical significant differentaily expressed genes (DEGs) based
#on the p adjusted value less than 0.05 and biological significance  based
#on the fold change more than 1
deseq.deg=res[res$padj < 0.05 & abs(res$log2FoldChange)>0.5,]

#export the Degs into your current folder for further analysthis
write.csv(as.matrix(deseq.deg),file="deseq.deg.csv", quote=F,row.names=T)

#drow DEGs volcano plot

par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(padj), pch=20, main="DECEASED vs LIVINGDEG padj<0.5 &FC0.5"))
with(subset(res, padj<0.05 & (log2FoldChange)> 0.5 ), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(res, padj<0.05 & (log2FoldChange)< -0.5), points(log2FoldChange, -log10(padj), pch=20, col="red"))
legend(x=3,y=20,c("Upregulated","Downgulated"), cex=.9, bty="n", col=c("blue","red"),pch=19)

par(mfrow=c(1,1))
with(res, plot(log2FoldChange, -log10(padj), pch=20, main="DECEASED vs LIVINGDEG padj<1.2 &FC1.2"))
with(subset(res, padj<0.05 & (log2FoldChange)> 1.2 ), points(log2FoldChange, -log10(padj), pch=20, col="blue"))
with(subset(res, padj<0.05 & (log2FoldChange)< -1.2), points(log2FoldChange, -log10(padj), pch=20, col="red"))
legend(x=3,y=20,c("Upregulated","Downgulated"), cex=.9, bty="n", col=c("blue","red"),pch=19)
###########################################################################################################################
####drow heatmap####
#normalize the data
####normalization after DESeq2
dds2 <- estimateSizeFactors(dds)
normalized_counts <- as.data.frame(counts(dds2, normalized=TRUE))
write.csv(as.matrix(normalized_counts),file="normalized_counts_mRNA.csv", quote=F,row.names=T)

#extract counts values of DEGs only for each stage for log FC >0.5
deseq.deg=res[res$padj < 0.05 & abs(res$log2FoldChange)>0.5,]
exp.degs=as.matrix(normalized_counts[rownames(normalized_counts) %in% rownames(deseq.deg), ])
aheatmap(log2(exp.degs+1), annCol =meta_data_new2$vital_status, col = rev(brewer.pal(9,"RdBu")), main=" DECEASED vs LIVINGDEG all DEGs FC>0.5")
dim(exp.degs) #943 473
#subset the exp.degs (first 10 genes across all samples) for better visualization of clusters 
new_exp.degs = as.matrix(exp.degs[1:10,])
dim(new_exp.degs) #10 473
aheatmap(log2(new_exp.degs+1), annCol =meta_data_new2$vital_status, col = rev(brewer.pal(9,"RdBu")), main=" DECEASED vs LIVINGDEG of top 10 genes")

#View(exp.degs)
#heatmap(log2(new_exp.degs+1))

######################################################
#extract counts values of DEGs only for each stage for log FC >1.2

deseq.deg=res[res$padj < 0.05 & abs(res$log2FoldChange)>1.2,]
exp.degs=as.matrix(normalized_counts[rownames(normalized_counts) %in% rownames(deseq.deg), ])
aheatmap(log2(exp.degs+1), annCol =meta_data_new2$vital_status, col = rev(brewer.pal(9,"RdBu")) , main=" DECEASED vs LIVINGDEG all DEGs FC>1.2")
dim(exp.degs) #118 473
#subset the exp.degs (first 10 genes across all samples) for better visualization of clusters 
new_exp.degs = as.matrix(exp.degs[1:10,])
dim(new_exp.degs) #10 473
aheatmap(log2(new_exp.degs+1), annCol =meta_data_new2$vital_status, col = rev(brewer.pal(9,"RdBu")), main=" DECEASED vs LIVINGDEG of top 10 genes")

#####cann't see any clusters of data so we need the integration analysis to see the hidden factors.

###################################################################################################
#Draw PCA 
#answer1
vsd = vst(dds, blind = FALSE)
se <- SummarizedExperiment(log2(counts(dds, normalized=FALSE) + 1),
                           colData=colData(dds))
plotPCA( DESeqTransform( se ), intgroup = "vital_status" )

vsd = vst(dds, blind = FALSE)
se <- SummarizedExperiment(log2(counts(dds, normalized=FALSE) + 1),
                           colData=colData(dds))
plotPCA( DESeqTransform( se ), intgroup = "gender" )

####no clustring of the data

#answer2
#PCA (vital_status)
#transpose rows with columns 
#count2 <- t(count)
#dim(count) #20531   473
#dim(count2) #473 20531
#PCA <- prcomp(count2)
#plot(PCA$x[,1], PCA$x[,2])
#plot(PCA)
#plot(PCA ,type="l")
#biplot(PCA , scale = 0)
# Extract PC scores
#pcaData <- as.data.frame(PCA$x[, 1:2])
#pcaData <- cbind(pcaData, meta_data_new2$vital_status) 
#colnames(pcaData) <- c("PC1", "PC2", "Species")
#warnings()
# plot with ggplot 
#library(ggplot2)
#ggplot(pcaData) +
#  aes(PC1, PC2, color = Species, shape = Species) + # define plot area
#  geom_point(size = 2) 
#adding data points
#percentVar <- round(100 * summary(PCA)$importance[2, 1:2], 0)
#ggplot(pcaData, aes(PC1, PC2, color = Species, shape = Species)) + # starting ggplot2
#  geom_point(size = 2) + # add data points
#  xlab(paste0("PC1: ", percentVar[1], "% variance")) + # x label
#  ylab(paste0("PC2: ", percentVar[2], "% variance")) + # y label
#  ggtitle("Principal component analysis (PCA)") + # title
#  theme(aspect.ratio = 1) # width and height ratio

#ggplot(pcaData, aes(PC1, PC2, col = Species, fill = Species)) +
#  stat_ellipse(geom = "polygon", col = "black", alpha = 0.5) +
#  geom_point(shape = 21, col = "black")

#PCA (gender)
# Extract PC scores
#pcaData <- as.data.frame(PCA$x[, 1:2])
#pcaData <- cbind(pcaData, meta_data_new2$gender) 
#colnames(pcaData) <- c("PC1", "PC2", "Gender")
# plot with ggplot 
#library(ggplot2)
#ggplot(pcaData) +
#  aes(PC1, PC2, color = Gender, shape = Gender) + # define plot area
#  geom_point(size =2) 
#adding data points
#percentVar <- round(100 * summary(PCA)$importance[2, 1:4], 0)
#ggplot(pcaData, aes(PC1, PC2, color = Gender, shape = Gender)) + # starting ggplot2
#  geom_point(size = 2) + # add data points
#  xlab(paste0("PC1: ", percentVar[1], "% variance")) + # x label
#  ylab(paste0("PC2: ", percentVar[2], "% variance")) + # y label
#  ggtitle("Principal component analysis (PCA)") + # title
#  theme(aspect.ratio = 1) # width and height ratio

#ggplot(pcaData, aes(PC1, PC2, col = Gender, fill = Gender)) +
#  stat_ellipse(geom = "polygon", col = "black", alpha = 0.5) +
#  geom_point(shape = 21, col = "black")

