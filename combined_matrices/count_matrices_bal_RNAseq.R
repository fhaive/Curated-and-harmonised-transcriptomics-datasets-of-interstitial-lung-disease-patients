#########make the consensus count matrices############
#These paths in this bunch of code doesnt exist in this folder. They are psudopaths
#symbol_norm folder contains all the RNAseq bal expression matrices and all the other 
#folders in this script are under this folder


setwd("/Lung_data/BAL/RNA-seq/symbol_norm")


GSE166036 <- read.table("Normalized_expression_matrix_BAL_symbol_GSE166036.txt", sep="\t", header = T)
rownames(GSE166036)<-GSE166036$X
GSE166036<-GSE166036[,2:13] #Exclude the SSc samples
library(xlsx)
metadata_GSE166036<-read.xlsx("/Lung_data/BAL/RNA-seq/GSE166036/phenodata/GSE166036_curated.xlsx", sheetIndex = 1)
metadata_GSE166036<-metadata_GSE166036[metadata_GSE166036$disease%in%c("IPF", "healthy"),]
metadata_GSE166036<-metadata_GSE166036[metadata_GSE166036$tissue_source%in%c("BAL"),]

####Disease and healthy

disease_GSE166036_subset<-GSE166036[,metadata_GSE166036$disease=="IPF"]

healthy_GSE166036_subset<-GSE166036[,metadata_GSE166036$disease=="healthy"]

setwd("/Lung_data/BAL/RNA-seq/symbol_norm/disease")
write.table(disease_GSE166036_subset, file="GSE166036_symbol_expression_matrix_disease_subset.txt", sep="\t")

setwd("/Lung_data/BAL/RNA-seq/symbol_norm/healthy")
write.table(healthy_GSE166036_subset, file="GSE166036_symbol_expression_matrix_healthy_subset.txt", sep="\t")

