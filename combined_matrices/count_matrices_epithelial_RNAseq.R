#########make the consensus count matrices############
#These paths in this bunch of code doesnt exist in this folder. They are psudopaths
#symbol_norm folder contains all the RNAseq epithelial expression matrices and all 
#the other folders in this script are under this folder



setwd("/Lung_data/Epithelial/RNA-seq/symbol_norm")


GSE151673<- read.table("Normalized_expression_matrix_symbol_GSE151673.txt", sep="\t", header = T)
rownames(GSE151673)<-GSE151673$X
GSE151673<-GSE151673[,2:length(colnames(GSE151673))]

GSE205525<- read.table("Normalized_expression_matrix_symbol_GSE205525.txt", sep="\t", header = T)
rownames(GSE205525)<-GSE205525$X
GSE205525<-GSE205525[,2:length(colnames(GSE205525))]


library(purrr)
common_row_names <- Reduce(intersect, list(rownames(GSE151673), rownames(GSE205525)))

GSE151673_subset <- GSE151673[rownames(GSE151673) %in% common_row_names, ]
GSE151673_subset<- GSE151673_subset[order(row.names(GSE151673_subset)), ]

GSE205525_subset <- GSE205525[rownames(GSE205525) %in% common_row_names, ]
GSE205525_subset<-GSE205525_subset[order(row.names(GSE205525_subset)), ]

setwd("/Lung_data/Epithelial/RNA-seq/symbol_norm/subsets")
write.table(GSE151673_subset, file="GSE151673_symbol_expression_matrix_subset.txt", sep="\t")
write.table(GSE205525_subset, file="GSE205525_symbol_expression_matrix_subset.txt", sep="\t")


####Disease and healthy

library(xlsx)

##################################################################
meta_GSE151673<-read.xlsx("/Lung_data/Epithelial/RNA-seq/GSE151673/phenodata/GSE151673_curated.xlsx", sheetIndex = 1)
disease_GSE151673<-meta_GSE151673$disease
disease_GSE151673_subset<-GSE151673_subset[,c(grep(pattern = "IPF", disease_GSE151673))]

healthy_GSE151673_subset<-GSE151673_subset[,c(grep(pattern = "healthy", disease_GSE151673))]

###################################################################
meta_GSE205525<-read.xlsx("/Lung_data/Epithelial/RNA-seq/GSE205525/phenodata/GSE205525_curated.xlsx", sheetIndex = 1)

disease_GSE205525_subset<- GSE205525_subset[,c(which(meta_GSE205525$sample_description=="IPF_no_fibroblast"))]

healthy_GSE205525_subset<-GSE205525_subset[,c(which(meta_GSE205525$sample_description=="healthy_no_fibroblast"))]


##################################################################################################
setwd("/Lung_data/Epithelial/RNA-seq/symbol_norm/subsets/disease")

write.table(disease_GSE151673_subset, file="GSE151673_symbol_expression_matrix_disease_subset.txt", sep="\t")
write.table(disease_GSE205525_subset, file="GSE205525_symbol_expression_matrix_disease_subset.txt", sep="\t")

setwd("/Lung_data/Epithelial/RNA-seq/symbol_norm/subsets/healthy")

write.table(healthy_GSE151673_subset, file="GSE151673_symbol_expression_matrix_healthy_subset.txt", sep="\t")
write.table(healthy_GSE205525_subset, file="GSE205525_symbol_expression_matrix_healthy_subset.txt", sep="\t")



setwd("/Lung_data/Epithelial/RNA-seq/symbol_norm/subsets/disease")

combined_table_disease <- cbind(disease_GSE151673_subset, disease_GSE205525_subset)
samples_disease <-colnames(combined_table_disease)
batch<-c(rep("GSE151673", length(colnames(disease_GSE151673_subset))), rep("GSE205525", length(colnames(disease_GSE205525_subset)))) 

write.table(combined_table_disease, file="combined_symbol_expression_matrix_disease_epithelial_rnaseq_subset.txt", sep="\t")
write.table(samples_disease, file="samples_epithelial_rnaseq_disease.txt", sep="\t")
write.table(batch, file="batch_epithelial_rnaseq_disease.txt", sep="\t")

setwd("/Lung_data/Epithelial/RNA-seq/symbol_norm/subsets/healthy")

combined_table_healthy <- cbind(healthy_GSE151673_subset, healthy_GSE205525_subset)
samples_healthy <-colnames(combined_table_healthy)
batch_healthy<-c((rep("GSE151673", length(colnames(healthy_GSE151673_subset)))), rep("GSE205525", length(colnames(healthy_GSE205525_subset)))) 

write.table(combined_table_healthy, file="combined_symbol_expression_matrix_healthy_epithelial_rnaseq_subset.txt", sep="\t")
write.table(samples_healthy, file="samples_epithelial_rnaseq_healthy.txt", sep="\t")
write.table(batch_healthy, file="batch_epithelial_rnaseq_healthy.txt", sep="\t")


