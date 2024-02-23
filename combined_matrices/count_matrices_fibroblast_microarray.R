#########make the consensus count matrices############
#These paths in this bunch of code doesnt exist in this folder. They are psudopaths
#symbol_norm contains all the microarray fibroblast expression matrices and all 
#the other folders in this script are under this folder


setwd("/Lung_data/Fibroblast/Microarray/symbol_norm")

files <- list.files(pattern = ".csv")
#files
#[1] "expression_matrix_symbol_GSE11196.csv"  "expression_matrix_symbol_GSE129164.csv" "expression_matrix_symbol_GSE144338.csv" "expression_matrix_symbol_GSE40839.csv" 
#[5] "expression_matrix_symbol_GSE44723.csv" 

matrix_list <- list()

for (i in 1:length(files)) {
  expr_mat <- read.table(file = files[i], sep = "\t", header = TRUE, row.names = 1)
  
  # Extract GSE ID from the filename using regular expression
  gse_id <- sub(".*GSE([0-9]+).*\\.csv", "\\1", files[i])
  
  # Use GSE ID as the name for the matrix in the list
  matrix_list[[paste("GSE", gse_id, sep = "")]] <- expr_mat
}

##########Metadatas########################

####Disease and healthy

####From now on I did not want to loop through the files because I wanted to make sure that just the samples that are wanted 
##for the downstream analysis are included. I had to subset the samples slightly differently based on the metadata to get just
#the samples that are IPF or UIP, biopsy samples. 


####Metadatas
###Base_path for metadata

base_path <- "/Lung_data/Fibroblast/Microarray"

# Construct file paths for each matrix in matrix_list
metadata_paths <- lapply(names(matrix_list), function(matrix_name) {
  path <- file.path(base_path, matrix_name, "phenodata", paste0(matrix_name, "_curated.xlsx"))
  return(path)
})

####The order in metadatas and matrix_list names(matrix_list)
#[1] "GSE11196"  "GSE129164" "GSE144338" "GSE40839"  "GSE44723" 

library(xlsx)
#####################################################################
meta_GSE11196<-read.xlsx(metadata_paths[[1]], sheetIndex = 1)
#Exctract the "total rna" samples from metadata GSE11196, already extracted in actual data
is_total_rna <- grepl("Total RNA", meta_GSE11196$title)

# Subset the dataframe to only include rows where is_total_rna is TRUE
meta_GSE11196 <- meta_GSE11196[is_total_rna, ]
#####################################################################################

meta_GSE44723<-read.xlsx(metadata_paths[[5]], sheetIndex = 1)
meta_GSE40839<-read.xlsx(metadata_paths[[4]], sheetIndex = 1)



#######subset the non-stimulated samples#######
meta_GSE129164<-read.xlsx(metadata_paths[[2]], sheetIndex = 1)
matrix_list[[2]]<-matrix_list[[2]][, meta_GSE129164$treatment%in%c("non_stimulated_24_hours_healthy", "non_stimulated_24_hours_IPF")]
meta_GSE129164<-meta_GSE129164[meta_GSE129164$treatment%in%c("non_stimulated_24_hours_healthy", "non_stimulated_24_hours_IPF"),]
###########################################################################################################
meta_GSE144338<-read.xlsx(metadata_paths[[3]], sheetIndex = 1)


library(purrr)
common_row_names <- Reduce(intersect, lapply(matrix_list, function(mat) rownames(mat)))



subset_matrices<-list()
for (i in 1:length(matrix_list)) {
  
  subset<-matrix_list[[i]][rownames(matrix_list[[i]]) %in% common_row_names,]
  subset_matrices[[paste0(names(matrix_list[i]), "_symbol_expression_matrix_subset")]] <- subset
  filename<-paste0(names(matrix_list[i]), "_symbol_expression_matrix_subset.txt")
  write.table(subset, file=paste0("subsets/", filename), sep="\t")
  
}


####Disease_samples#######

disease_GSE11196_subset<-subset_matrices[[1]][,meta_GSE11196$disease %in% c("IPF")]

disease_GSE40839_subset<-subset_matrices[[4]][,meta_GSE40839$disease %in% c("UIP")]

disease_GSE44723_subset<-subset_matrices[[5]][,meta_GSE44723$disease %in% c("IPF")]

disease_GSE129164_subset<-subset_matrices[[2]][,meta_GSE129164$disease %in% c("IPF")]

disease_GSE144338_subset<-subset_matrices[[3]][,meta_GSE144338$disease %in% c("IPF")]

setwd("/Lung_data/Fibroblast/Microarray/symbol_norm/subsets/disease")

write.table(disease_GSE11196_subset, file="GSE11196_symbol_expression_matrix_disease_subset.txt", sep="\t")
write.table(disease_GSE44723_subset, file="GSE44723_symbol_expression_matrix_disease_subset.txt", sep="\t")
write.table(disease_GSE40839_subset, file="GSE40839_symbol_expression_matrix_disease_subset.txt", sep="\t")
write.table(disease_GSE129164_subset, file="GSE129164_symbol_expression_matrix_disease_subset.txt", sep="\t")
write.table(disease_GSE144338_subset, file="GSE144338_symbol_expression_matrix_disease_subset.txt", sep="\t")


##################Healthy_samples####################
healthy_GSE11196_subset<-subset_matrices[[1]][,meta_GSE11196$disease %in% c("healthy")]

healthy_GSE40839_subset<-subset_matrices[[4]][,meta_GSE40839$disease %in% c("healthy")]

healthy_GSE44723_subset<-subset_matrices[[5]][,meta_GSE44723$disease %in% c("healthy")]

healthy_GSE129164_subset<-subset_matrices[[2]][,meta_GSE129164$disease %in% c("healthy")]

healthy_GSE144338_subset<-subset_matrices[[3]][,meta_GSE144338$disease %in% c("healthy")]

setwd("/Lung_data/Fibroblast/Microarray/symbol_norm/subsets/healthy")

write.table(healthy_GSE11196_subset, file="GSE11196_symbol_expression_matrix_healthy_subset.txt", sep="\t")
write.table(healthy_GSE44723_subset, file="GSE44723_symbol_expression_matrix_healthy_subset.txt", sep="\t")
write.table(healthy_GSE40839_subset, file="GSE40839_symbol_expression_matrix_healthy_subset.txt", sep="\t")
write.table(healthy_GSE129164_subset, file="GSE129164_symbol_expression_matrix_healthy_subset.txt", sep="\t")
write.table(healthy_GSE144338_subset, file="GSE144338_symbol_expression_matrix_healthy_subset.txt", sep="\t")


################Count_matrices########################################

setwd("/Lung_data/Fibroblast/Microarray/symbol_norm/subsets/disease")

combined_table_disease <- cbind(disease_GSE11196_subset, disease_GSE129164_subset,  disease_GSE144338_subset, disease_GSE40839_subset, disease_GSE44723_subset)
samples_disease <-colnames(combined_table_disease)
batch<-c(rep("GSE11196", length(colnames(disease_GSE11196_subset))), rep("GSE129164", length(colnames(disease_GSE129164_subset))),  rep("GSE144338", length(colnames(disease_GSE144338_subset))), rep("GSE40839", length(colnames(disease_GSE40839_subset))), rep("GSE44723", length(colnames(disease_GSE44723_subset)))) 

write.table(combined_table_disease, file="combined_symbol_expression_matrix_disease_fibroblast_microarray_subset.txt", sep="\t")
write.table(samples_disease, file="samples_fibroblast_microarray_disease.txt", sep="\t")
write.table(batch, file="batch_fibroblast_microarray_disease.txt", sep="\t")

setwd("/Lung_data/Fibroblast/Microarray/symbol_norm/subsets/healthy")

combined_table_healthy <- cbind(healthy_GSE11196_subset, healthy_GSE129164_subset,  healthy_GSE144338_subset, healthy_GSE40839_subset, healthy_GSE44723_subset)
samples_healthy <-colnames(combined_table_healthy)
batch_healthy<-c(rep("GSE11196", length(colnames(healthy_GSE11196_subset))), rep("GSE129164", length(colnames(healthy_GSE129164_subset))),  rep("GSE144338", length(colnames(healthy_GSE144338_subset))), rep("GSE40839", length(colnames(healthy_GSE40839_subset))), rep("GSE44723", length(colnames(healthy_GSE44723_subset)))) 

write.table(combined_table_healthy, file="combined_symbol_expression_matrix_healthy_fibroblast_microarray_subset.txt", sep="\t")
write.table(samples_healthy, file="samples_fibroblast_microarray_healthy.txt", sep="\t")
write.table(batch_healthy, file="batch_fibroblast_microarray_healthy.txt", sep="\t")
