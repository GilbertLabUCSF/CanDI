#!/usr/bin/env Rscript
library(DESeq2)

args = commandArgs(trailingOnly=TRUE)

#read data
counts.mat <- read.csv(args[1])
coldata <- read.csv(args[2])

#name metadata columns
colnames(coldata) <- c("sample.id", "condition")
#convert datatype to factor
coldata[,-1] <- as.factor(coldata[, -1])
#Match sample ids to sample columns in counts matrix
coldata$sample.id <- sub("-", ".", coldata$sample.id)

#init dds object
dds <- DESeqDataSetFromMatrix(countData=counts.mat, 
                              colData=coldata, 
                              design= ~condition,
                              tidy = TRUE)

dds <- estimateSizeFactors(dds)
#Filter lowly expressed genes
idx <- rowSums(counts(dds, normalized=TRUE) >= 5) >= 3
dds <- dds[idx,]

dds <- DESeq(dds) #run deseq
res <- results(dds) #get results
#Show results
print(summary(res))
print(head(res))

write.csv(res, args[3]) #save results

