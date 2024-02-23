#########make the consensus count matrices############

#These paths in this bunch of code doesnt exist in this folder. They are psudopaths
#This folder contains all the RNAseq fibroblast expression matrices and
#all the other folders in this script are under this folder


setwd("/Lung_data/Fibroblast/RNA-seq/symbol_norm")


GSE185492 <- read.table("Normalized_expression_matrix_ALL_symbol_GSE185492.txt", sep="\t", header = T)
rownames(GSE185492)<-GSE185492$X
GSE185492<-GSE185492[,2:length(colnames(GSE185492))]

GSE180415<- read.table("Normalized_expression_matrix_symbol_GSE180415.txt", sep="\t", header = T)
rownames(GSE180415)<-GSE180415$X
GSE180415<-GSE180415[,2:length(colnames(GSE180415))]

GSE97038 <- read.table("Normalized_expression_matrix_ALL_symbol_GSE97038.txt", sep="\t", header = T)
rownames(GSE97038)<-GSE97038$X
GSE97038<-GSE97038[,2:length(colnames(GSE97038))]


library(purrr)
common_row_names <- Reduce(intersect, list(rownames(GSE185492), rownames(GSE180415), rownames(GSE97038)))

GSE185492_subset <- GSE185492[rownames(GSE185492) %in% common_row_names, ]
GSE185492_subset<- GSE185492_subset[order(row.names(GSE185492_subset)), ]

GSE180415_subset <- GSE180415[rownames(GSE180415) %in% common_row_names, ]
GSE180415_subset<-GSE180415_subset[order(row.names(GSE180415_subset)), ]

GSE97038_subset <-GSE97038[rownames(GSE97038) %in% common_row_names, ]
GSE97038_subset<-GSE97038_subset[order(row.names(GSE97038_subset)), ]

setwd("/Lung_data/Fibroblast/RNA-seq/symbol_norm/subsets")

write.table(GSE185492_subset, file="GSE185492_symbol_expression_matrix_subset.txt", sep="\t")
write.table(GSE180415_subset, file="GSE180415_symbol_expression_matrix_subset.txt", sep="\t")
write.table(GSE97038_subset, file="GSE97038_symbol_expression_matrix_subset.txt", sep="\t")

####Disease and healthy

library(xlsx)

##################################################################
meta_GSE185492<-read.xlsx("/Lung_data/Fibroblast/RNA-seq/GSE185492/phenodata/GSE185492_curated.xlsx", sheetIndex = 1)
disease_GSE185492<-meta_GSE185492$disease
disease_GSE185492_subset<-GSE185492_subset[,c(grep(pattern = "IPF", disease_GSE185492))]

healthy_GSE185492_subset<-GSE185492_subset[,c(grep(pattern = "healthy", disease_GSE185492))]

###################################################################
meta_GSE180415<-read.xlsx("/Lung_data/Fibroblast/RNA-seq/GSE180415/phenodata/GSE180415_curated.xlsx", sheetIndex = 1)
meta_GSE180415<-meta_GSE180415[meta_GSE180415$disease%in%c("IPF", "healthy"),]
disease_GSE180415<-meta_GSE180415$disease

disease_GSE180415_subset<- GSE180415_subset[,c(grep(pattern = "IPF", disease_GSE180415))]

healthy_GSE180415_subset<-GSE180415_subset[,c(grep(pattern = "healthy", disease_GSE180415))]

##############################################################################
meta_GSE97038<-read.xlsx("/Lung_data/Fibroblast/RNA-seq/GSE97038/phenodata/GSE97038_curated.xlsx", sheetIndex = 1)
disease_GSE97038<-meta_GSE97038$disease


disease_GSE97038_subset<- GSE97038_subset[,c(grep(pattern = "IPF", disease_GSE97038))]

healthy_GSE97038_subset<-GSE97038_subset[,c(grep(pattern = "healthy", disease_GSE97038))]
#################################################################################

setwd("/Lung_data/Fibroblast/RNA-seq/symbol_norm/subsets/disease")

write.table(disease_GSE185492_subset, file="GSE185492_symbol_expression_matrix_disease_subset.txt", sep="\t")
write.table(disease_GSE180415_subset, file="GSE180415_symbol_expression_matrix_disease_subset.txt", sep="\t")
write.table(disease_GSE97038_subset, file="GSE97038_symbol_expression_matrix_disease_subset.txt", sep="\t")

setwd("/Lung_data/Fibroblast/RNA-seq/symbol_norm/subsets/healthy")

write.table(healthy_GSE185492_subset, file="GSE185492_symbol_expression_matrix_healthy_subset.txt", sep="\t")
write.table(healthy_GSE180415_subset, file="GSE180415_symbol_expression_matrix_healthy_subset.txt", sep="\t")
write.table(healthy_GSE97038_subset, file="GSE97038_symbol_expression_matrix_healthy_subset.txt", sep="\t")


setwd("/Lung_data/Fibroblast/RNA-seq/symbol_norm/subsets/disease")

combined_table_disease <- cbind(disease_GSE185492_subset, disease_GSE180415_subset,  disease_GSE97038_subset)
samples_disease <-colnames(combined_table_disease)
batch<-c(rep("GSE185492", length(colnames(disease_GSE185492_subset))), rep("GSE180415", length(colnames(disease_GSE180415_subset))),  rep("GSE97038", length(colnames(disease_GSE97038_subset)))) 

write.table(combined_table_disease, file="combined_symbol_expression_matrix_disease_fibroblast_rnaseq_subset.txt", sep="\t")
write.table(samples_disease, file="samples_biopsy_rnaseq_disease.txt", sep="\t")
write.table(batch, file="batch_biopsy_rnaseq_disease.txt", sep="\t")

setwd("/Lung_data/Fibroblast/RNA-seq/symbol_norm/subsets/healthy")

combined_table_healthy <- cbind(healthy_GSE185492_subset, healthy_GSE180415_subset, healthy_GSE97038_subset)
samples_healthy <-colnames(combined_table_healthy)
batch_healthy<-c((rep("GSE185492", length(colnames(healthy_GSE185492_subset)))), rep("GSE180415", length(colnames(healthy_GSE180415_subset))),rep("GSE97038", length(colnames(healthy_GSE97038_subset)))) 

write.table(combined_table_healthy, file="combined_symbol_expression_matrix_healthy_fibroblast_rnaseq_subset.txt", sep="\t")
write.table(samples_healthy, file="samples_biopsy_rnaseq_healthy.txt", sep="\t")
write.table(batch_healthy, file="batch_biopsy_rnaseq_healthy.txt", sep="\t")


