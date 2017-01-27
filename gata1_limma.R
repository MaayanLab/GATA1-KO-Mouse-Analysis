
library(affy)
library(Biobase)
library(edgeR)

setwd("~/Documents/GATA1-Analysis")

marrowdf <-  read.csv('GATA1_marrow_profile.txt',sep='\t',header=TRUE,row.names=1)
marrow <- ExpressionSet(assayData=as.matrix(marrowdf))

spleendf <-  read.csv('GATA1_spleen_profile.txt',sep='\t',header=TRUE,row.names=1)
spleen <- ExpressionSet(assayData=as.matrix(spleendf))
# m <- as.matrix(marrowdf)
# marrow <- new("ExpressionSet",expers=marrowdf)

# 
# spleendf <-  read.csv('GATA1_spleen_profile.txt',sep='\t',header=TRUE,row.names=1)
# eset <- AnnotatedDataFrame(spleendf)
# files <- targets[,c("GATA1_marrow_profile.txt","GATA1_spleen_profile.txt")]
# RG <- read.maimages(files, source="imagene")

# marrow <- new("ExpressionSet",expers=marrow,)