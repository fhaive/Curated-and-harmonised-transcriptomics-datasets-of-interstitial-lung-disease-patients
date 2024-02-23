#########make the consensus count matrices############

#These paths in this bunch of code doesnt exist in this folder. They are psudopaths
#symbol_norm folder contains all the microarray biopsy expression matrices and all the other folders in this script are under this folder

setwd("/Lung_data/Biopsy/Microarray/symbol_norm")

files <- list.files(pattern = ".csv")
#files
#[1] "expression_matrix_symbol_GSE10667.csv"  "expression_matrix_symbol_GSE110147.csv" "expression_matrix_symbol_GSE21369.csv"  "expression_matrix_symbol_GSE24206.csv" 
#[5] "expression_matrix_symbol_GSE53845.csv"  "expression_matrix_symbol_GSE72073.csv" 

matrix_list <- list()

for (i in 1:length(files)) {
  expr_mat <- read.table(file = files[i], sep = "\t", header = TRUE, row.names = 1)
  
  # Extract GSE ID from the filename using regular expression
  gse_id <- sub(".*GSE([0-9]+).*\\.csv", "\\1", files[i])
  
  # Use GSE ID as the name for the matrix in the list
  matrix_list[[paste("GSE", gse_id, sep = "")]] <- expr_mat
}

library(purrr)
common_row_names <- Reduce(intersect, lapply(matrix_list, function(mat) rownames(mat)))



subset_matrices<-list()
for (i in 1:length(matrix_list)) {
  
  subset<-matrix_list[[i]][rownames(matrix_list[[i]]) %in% common_row_names,]
  subset_matrices[[paste0(names(matrix_list[i]), "_symbol_expression_matrix_subset")]] <- subset
  filename<-paste0(names(matrix_list[i]), "_symbol_expression_matrix_subset.txt")
  write.table(subset, file=paste0("subsets/", filename), sep="\t")
  
}



####Disease and healthy

####From now on I did not want to loop through the files because I wanted to make sure that just the samples that are wanted 
##for the downstream analysis are included. I had to subset the samples slightly differently based on the metadata to get just
#the samples that are IPF or UIP, biopsy samples. 

###Base_path for metadata

base_path <- "/Lung_data/Biopsy/Microarray"

# Construct file paths for each matrix in matrix_list
metadata_paths <- lapply(names(matrix_list), function(matrix_name) {
  path <- file.path(base_path, matrix_name, "phenodata", paste0(matrix_name, "_curated.xlsx"))
  return(path)
})




####THE ORDER IN metadata_paths and subset_matrices from indexes 1 to 11 is 
#"GSE10667"  "GSE110147" "GSE21369"  "GSE24206"  "GSE53845"  "GSE72073" 

#!!!!!!!!!Check the order before proceeding!!!!!!

#The reason of the order of the datasets below is just that I had done this previously in this order so 
#They happen to be now in this order


####Disease and healthy

library(xlsx)
meta_GSE110147<-read.xlsx(metadata_paths[[2]], sheetIndex=1)
meta_GSE21369<-read.xlsx(metadata_paths[[3]], sheetIndex=1)
meta_GSE24206<-read.xlsx(metadata_paths[[4]], sheetIndex=1)
meta_GSE53845<-read.xlsx(metadata_paths[[5]], sheetIndex=1)
meta_GSE72073<-read.xlsx(metadata_paths[[6]], sheetIndex=1)
meta_GSE10667<-read.xlsx(metadata_paths[[1]], sheetIndex=1)



disease_GSE110147<-meta_GSE110147$disease
disease_GSE21369<-meta_GSE21369$disease
disease_GSE24206<-meta_GSE24206$disease
#disease_GSE53845<-meta_GSE53845$disease #Comes on later
disease_GSE72073<-meta_GSE72073$disease
disease_GSE10667<-meta_GSE10667$disease

disease_subsets<-list()

healthy_subsets<-list()

###########################################################################################################
disease_GSE110147_subset<-subset_matrices[[2]][,1:22] # First 22 samples IPF, samples 33-37 mixed_IPF_NSIP
disease_subsets[["GSE110147_symbol_expression_matrix_disease_subset"]] <- disease_GSE110147_subset
healthy_GSE110147_subset<-subset_matrices[[2]][,c(grep(pattern = "healthy", disease_GSE110147))]
healthy_subsets[["GSE110147_symbol_expression_matrix_healthy_subset"]] <- healthy_GSE110147_subset
#############################################################################################3
disease_GSE21369_subset<- subset_matrices[[3]][,c(grep(pattern = "UIP", disease_GSE21369))]
disease_subsets[["GSE21369_symbol_expression_matrix_disease_subset"]] <- disease_GSE21369_subset
healthy_GSE21369_subset<- subset_matrices[[3]][,c(grep(pattern = "healthy", disease_GSE21369))]
healthy_subsets[["GSE21369_symbol_expression_matrix_healthy_subset"]] <- healthy_GSE21369_subset
#######################################################################

disease_GSE24206_subset<- subset_matrices[[4]][,c(grep(pattern = "IPF", disease_GSE24206))]
disease_subsets[["GSE24206_symbol_expression_matrix_disease_subset"]] <- disease_GSE24206_subset
healthy_GSE24206_subset<- subset_matrices[[4]][,c(grep(pattern = "healthy", disease_GSE24206))]
healthy_subsets[["GSE24206_symbol_expression_matrix_healthy_subset"]] <- healthy_GSE24206_subset
###########################################################################################

#GSE53845 metadata and actual data in different order. Reorder

subset_column_order <- colnames(subset_matrices[[5]])

# Match the values in meta_GSE53845$title with subset_column_order
matched_indexes <- match(subset_column_order,meta_GSE53845$title)

meta_GSE53845<-meta_GSE53845[matched_indexes,]

disease_GSE53845<-meta_GSE53845$disease

disease_GSE53845_subset<- subset_matrices[[5]][,c(grep(pattern = "IPF", disease_GSE53845))]
disease_subsets[["GSE53845_symbol_expression_matrix_disease_subset"]] <- disease_GSE53845_subset
healthy_GSE53845_subset<- subset_matrices[[5]][,c(grep(pattern = "healthy", disease_GSE53845))]
healthy_subsets[["GSE53845_symbol_expression_matrix_healthy_subset"]] <- healthy_GSE53845_subset
#########################################################################################################
disease_GSE72073_subset<- subset_matrices[[6]][,c(grep(pattern = "IPF", disease_GSE72073))]
disease_subsets[["GSE72073_symbol_expression_matrix_disease_subset"]] <- disease_GSE72073_subset
healthy_GSE72073_subset<- subset_matrices[[6]][,c(grep(pattern = "primary_spontaneous_pneumothorax", disease_GSE72073))]
healthy_subsets[["GSE72073_symbol_expression_matrix_healthy_subset"]] <- healthy_GSE72073_subset

############################################################################################################
disease_GSE10667_subset<- subset_matrices[[1]][,c(grep(pattern = "UIP", disease_GSE10667))]
disease_subsets[["GSE10667_symbol_expression_matrix_disease_subset"]] <- disease_GSE10667_subset
healthy_GSE10667_subset<- subset_matrices[[1]][,c(grep(pattern = "healthy", disease_GSE10667))]
healthy_subsets[["GSE10667_symbol_expression_matrix_healthy_subset"]] <- healthy_GSE10667_subset

##########################################################################################

######################################################################################
###############################################################################################

disease_subsets_ordered<-list()
for (i in 1:length(disease_subsets)) {
  subset_ordered<-disease_subsets[[i]][order(rownames(disease_subsets[[i]])),]
  disease_subsets_ordered[[names(disease_subsets[i])]] <- subset_ordered
  
}

healthy_subsets_ordered<-list()
for (i in 1:length(healthy_subsets)) {
  subset_ordered<-healthy_subsets[[i]][order(rownames(healthy_subsets[[i]])),]
  healthy_subsets_ordered[[names(healthy_subsets[i])]] <- subset_ordered
  
}

setwd("/Lung_data/Biopsy/Microarray/symbol_norm/subsets/disease")

for (i in 1:length(disease_subsets)) {
  write.table(disease_subsets_ordered[[i]], file=paste0(names(disease_subsets_ordered[i]),".txt"), sep="\t")
}


setwd("/Lung_data/Biopsy/Microarray/symbol_norm/subsets/healthy")

for (i in 1:length(healthy_subsets_ordered)) {
  write.table(healthy_subsets_ordered[[i]], file=paste0(names(healthy_subsets_ordered[i]),".txt"), sep="\t")
}


setwd("/Lung_data/Biopsy/Microarray/symbol_norm/subsets/disease")

combined_table_disease <- do.call(cbind, disease_subsets_ordered)
names(combined_table_disease) <- sub("^[^.]+\\.", "", names(combined_table_disease))
samples_disease<-colnames(combined_table_disease)

batch<-c()
for (i in 1:length(names(disease_subsets_ordered))) {
  GSE<-sub("^GSE([0-9]+)_.*", "GSE\\1", names(disease_subsets)[i])
  batch_individual<-rep(GSE, length(colnames(disease_subsets_ordered[[i]])))
  batch<-c(batch, batch_individual)
}


write.table(combined_table_disease, file="combined_symbol_expression_matrix_disease_biopsy_rnaseq_subset.txt", sep="\t")
write.table(samples_disease, file="samples_biopsy_rnaseq_disease.txt", sep="\t")
write.table(batch, file="batch_biopsy_rnaseq_disease.txt", sep="\t")

setwd("/Lung_data/Biopsy/Microarray/symbol_norm/subsets/healthy")

combined_table_healthy <- do.call(cbind, healthy_subsets_ordered)
names(combined_table_healthy) <- sub("^[^.]+\\.", "", names(combined_table_healthy))
samples_healthy<-colnames(combined_table_healthy)

batch_healthy<-c()
for (i in 1:length(names(healthy_subsets_ordered))) {
  GSE<-sub("^GSE([0-9]+)_.*", "GSE\\1", names(healthy_subsets)[i])
  batch_individual<-rep(GSE, length(colnames(healthy_subsets_ordered[[i]])))
  batch_healthy<-c(batch_healthy, batch_individual)
}


write.table(combined_table_healthy, file="combined_symbol_expression_matrix_healthy_biopsy_rnaseq_subset.txt", sep="\t")
write.table(samples_healthy, file="samples_biopsy_rnaseq_healthy.txt", sep="\t")
write.table(batch_healthy, file="batch_biopsy_rnaseq_healthy.txt", sep="\t")
