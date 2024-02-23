#########make the consensus count matrices############
#These paths in this bunch of code doesnt exist in this folder. They are psudopaths
#This folder contains all the microarray macrophage expression matrices and
#all the other folders in this script are under this folder

setwd("/Lung_data/Macrophage/Microarray/symbol_norm")


GSE49072 <- read.table("expression_matrix_symbol_GSE49072.csv", sep="\t")
GSE90010<- read.table("expression_matrix_GSE90010_symbol.csv", sep="\t")

meta_GSE49072<-read.xlsx("/Lung_data/Macrophage/Microarray/GSE49072/phenodata/GSE49072_curated.xlsx", sheetIndex = 1)
meta_GSE90010<-read.xlsx("/Lung_data/Macrophage/Microarray/GSE90010/phenodata/GSE90010_curated.xlsx", sheetIndex = 1)


library(purrr)
common_row_names <- Reduce(intersect, list(rownames(GSE49072),rownames(GSE90010)))

GSE49072_subset <- GSE49072[rownames(GSE49072) %in% common_row_names, ]
GSE90010_subset <- GSE90010[rownames(GSE90010) %in% common_row_names, ]

setwd("/Lung_data/Macrophage/Microarray/symbol_norm/subsets")

write.table(GSE49072_subset, file="GSE49072_symbol_expression_matrix_subset.txt", sep="\t")
write.table(GSE90010_subset, file="GSE90010_symbol_expression_matrix_subset.txt", sep="\t")


####Disease_samples#######

disease_GSE49072_subset<-GSE49072_subset[,meta_GSE49072$disease_state %in% c("spontaneous_IPF")]

disease_GSE90010_subset<-GSE90010_subset[,meta_GSE90010$cell_line_treatment %in% c("AM_from_IPF")]

setwd("/Lung_data/Macrophage/Microarray/symbol_norm/subsets/disease")

write.table(disease_GSE49072_subset, file="GSE49072_symbol_expression_matrix_disease_subset.txt", sep="\t")
write.table(disease_GSE90010_subset, file="GSE90010_symbol_expression_matrix_disease_subset.txt", sep="\t")

##################Healthy_samples####################
healthy_GSE49072_subset<-GSE49072_subset[,meta_GSE49072$disease_state %in% c("healthy_volunteer")]

healthy_GSE90010_subset<-GSE90010_subset[,meta_GSE90010$cell_line_treatment %in% c("MDM_no_treatment")]

setwd("/Lung_data/Macrophage/Microarray/symbol_norm/subsets/healthy")

write.table(healthy_GSE49072_subset, file="GSE49072_symbol_expression_matrix_healthy_subset.txt", sep="\t")
write.table(healthy_GSE90010_subset, file="GSE90010_symbol_expression_matrix_healthy_subset.txt", sep="\t")


################Count_matrices########################################

setwd("/Lung_data/Macrophage/Microarray/symbol_norm/subsets/disease")

combined_table_disease <- cbind(disease_GSE49072_subset, disease_GSE90010_subset)
samples_disease <-colnames(combined_table_disease)
batch<-c(rep("GSE49072", length(colnames(disease_GSE49072_subset))), rep("GSE90010", length(colnames(disease_GSE90010_subset)))) 

write.table(combined_table_disease, file="combined_symbol_expression_matrix_disease_macrophage_microarray_subset.txt", sep="\t")
write.table(samples_disease, file="samples_macrophage_microarray_disease.txt", sep="\t")
write.table(batch, file="batch_macrophage_microarray_disease.txt", sep="\t")

setwd("/Lung_data/Macrophage/Microarray/symbol_norm/subsets/healthy")

combined_table_healthy <- cbind(healthy_GSE49072_subset, healthy_GSE90010_subset)
samples_healthy <-colnames(combined_table_healthy)
batch_healthy<-c(rep("GSE49072", length(colnames(healthy_GSE49072_subset))), rep("GSE90010", length(colnames(healthy_GSE90010_subset)))) 

write.table(combined_table_healthy, file="combined_symbol_expression_matrix_healthy_macrophage_microarray_subset.txt", sep="\t")
write.table(samples_healthy, file="samples_macrophage_microarray_healthy.txt", sep="\t")
write.table(batch_healthy, file="batch_macrophage_microarray_healthy.txt", sep="\t")
