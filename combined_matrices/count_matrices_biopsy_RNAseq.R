#########make the consensus count matrices############

#These paths in this bunch of code doesnt exist in this folder. They are psudopaths
#symbol_norm folder contains all the rnaseq biopsy expression matrices and all the other folders 
#in this script are under this folder

setwd("/Lung_data/Biopsy/RNA-seq/symbol_norm")

files <- list.files(pattern = ".txt")

#files
#[1] "Normalized_expression_matrix_symbol_GSE199949.txt"    "Normalized_expression_matrix_symbol_GSE213001.txt"    "Normalized_expression_matrix_biopsy_symbol_GSE166036.txt"
#[4] "Normalized_expression_matrix_symbol_GSE124685.txt"        "Normalized_expression_matrix_symbol_GSE138283.txt"        "Normalized_expression_matrix_symbol_GSE150910.txt"       
#[7] "Normalized_expression_matrix_symbol_GSE169500.txt"        "Normalized_expression_matrix_symbol_GSE184316.txt"        "Normalized_expression_matrix_symbol_GSE199152.txt"       
#[10] "Normalized_expression_matrix_symbol_GSE92592.txt"         "Normalized_expression_matrix_symbol_GSE99621.txt"        


matrix_list <- list()

for (i in 1:length(files)) {
  expr_mat <- read.table(file = files[i], sep = "\t", header = TRUE, row.names = 1)
  
  # Extract GSE ID from the filename using regular expression
  gse_id <- sub(".*GSE([0-9]+).*\\.txt", "\\1", files[i])
  
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
#the samples that are IPF or UIP, biopsy samples. For example in GSE166036 the data is also subsetted based on tissue source
#unlike other datasets
#It could have been done with loops but I wanted to go through all the 
#dataset manually in this phase to make it correct

###Base_path for metadata

base_path <- "/nasdata/sinkala/Lung_data/Biopsy/RNA-seq"

# Construct file paths for each matrix in matrix_list
metadata_paths <- lapply(names(matrix_list), function(matrix_name) {
  path <- file.path(base_path, matrix_name, "phenodata", paste0(matrix_name, "_curated.xlsx"))
  return(path)
})


disease_subsets<-list()

healthy_subsets<-list()

####THE ORDER IN metadata_paths and subset_matrices from indexes 1 to 11 is 
#GSE199949, GSE213001, GSE166036, GSE124685, GSE138283, GSE150910, GSE169500, GSE184316,
#GSE199152, GSE92592, GSE99621. 

#!!!!!!!!!Check the order before proceeding!!!!!!

#The reason of the order of the datasets below is just that I had done this previously in this order so 
#They happen to be now in this order

##################################################################
meta_GSE124685<-read.xlsx(metadata_paths[[4]], sheetIndex = 1)
disease_GSE124685<-meta_GSE124685$disease
disease_GSE124685_subset<-subset_matrices[[4]][,c(grep(pattern = "IPF", disease_GSE124685))]
disease_subsets[["GSE124685_symbol_expression_matrix_disease_subset"]] <- disease_GSE124685_subset
healthy_GSE124685_subset<-subset_matrices[[4]][,c(grep(pattern = "healthy", disease_GSE124685))]
healthy_subsets[["GSE124685_symbol_expression_matrix_healthy_subset"]] <- healthy_GSE124685_subset


#####################################################################################################
meta_GSE150910<-read.xlsx(metadata_paths[[6]], sheetIndex = 1)
meta_GSE150910<-meta_GSE150910[meta_GSE150910$disease%in%c("IPF", "healthy"),]
disease_GSE150910<-meta_GSE150910$disease

disease_GSE150910_subset<- subset_matrices[[6]][,c(grep(pattern = "ipf", colnames(subset_matrices[[6]])))]
disease_subsets[["GSE150910_symbol_expression_matrix_disease_subset"]] <- disease_GSE150910_subset
healthy_GSE150910_subset<- subset_matrices[[6]][,c(grep(pattern = "control", colnames(subset_matrices[[6]])))]
healthy_subsets[["GSE150910_symbol_expression_matrix_healthy_subset"]] <- healthy_GSE150910_subset
#################################################################################################

meta_GSE166036<-read.xlsx(metadata_paths[[3]], sheetIndex = 1)
meta_GSE166036<-meta_GSE166036[meta_GSE166036$disease%in%c("IPF", "healthy"),]
meta_GSE166036<-meta_GSE166036[meta_GSE166036$tissue_source%in%c("biopsy"),]
disease_GSE166036<-meta_GSE166036$disease

disease_GSE166036_subset<- subset_matrices[[3]][,c(grep(pattern = "IPF_Biopsy", colnames(subset_matrices[[3]])))]
disease_subsets[["GSE166036_symbol_expression_matrix_disease_subset"]] <- disease_GSE166036_subset
healthy_GSE166036_subset<- subset_matrices[[3]][,c(grep(pattern = "Normal_Biopsy", colnames(subset_matrices[[3]])))]
healthy_subsets[["GSE166036_symbol_expression_matrix_healthy_subset"]] <- healthy_GSE166036_subset

###################################################################

meta_GSE169500<-read.xlsx(metadata_paths[[7]], sheetIndex = 1)
disease_GSE169500<-meta_GSE169500$disease

disease_GSE169500_subset<- subset_matrices[[7]][,c(grep(pattern = "IPF", disease_GSE169500))]
disease_subsets[["GSE169500_symbol_expression_matrix_disease_subset"]] <- disease_GSE169500_subset
healthy_GSE169500_subset<-subset_matrices[[7]][,c(grep(pattern = "healthy", disease_GSE169500))]
healthy_subsets[["GSE169500_symbol_expression_matrix_healthy_subset"]] <- healthy_GSE169500_subset
######################################################################################

meta_GSE184316<-read.xlsx(metadata_paths[[8]], sheetIndex = 1)
meta_GSE184316<-meta_GSE184316[meta_GSE184316$disease%in%c("IPF", "healthy"),]
disease_GSE184316<-meta_GSE184316$disease

disease_GSE184316_subset<- subset_matrices[[8]][,c(grep(pattern = "IPF", colnames(subset_matrices[[8]])))]
disease_subsets[["GSE184316_symbol_expression_matrix_disease_subset"]] <- disease_GSE184316_subset
healthy_GSE184316_subset<- subset_matrices[[8]][,c(grep(pattern = "CTRL", colnames(subset_matrices[[8]])))]
healthy_subsets[["GSE184316_symbol_expression_matrix_healthy_subset"]] <- healthy_GSE184316_subset

###########################################################################################

meta_GSE199152<-read.xlsx(metadata_paths[[9]], sheetIndex = 1)
meta_GSE199152<-meta_GSE199152[meta_GSE199152$title%in%colnames(subset_matrices[[9]]),]
meta_GSE199152<-meta_GSE199152[meta_GSE199152$disease%in%c("UIP", "healthy"),]
disease_GSE199152<-meta_GSE199152$disease


disease_GSE199152_subset<- subset_matrices[[9]][,c(grep(pattern = "UIP", disease_GSE199152))]
disease_subsets[["GSE199152_symbol_expression_matrix_disease_subset"]] <- disease_GSE199152_subset
healthy_GSE199152_subset<- subset_matrices[[9]][,c(grep(pattern = "healthy",  disease_GSE199152))]
healthy_subsets[["GSE199152_symbol_expression_matrix_healthy_subset"]] <- healthy_GSE199152_subset




##############################################################################

meta_GSE199949<-read.xlsx(metadata_paths[[1]], sheetIndex = 1)
disease_GSE199949<-meta_GSE199949$disease
disease_GSE199949_subset<- subset_matrices[[1]][,c(grep(pattern = "IPF", disease_GSE199949))]
disease_subsets[["GSE199949_symbol_expression_matrix_disease_subset"]] <- disease_GSE199949_subset
healthy_GSE199949_subset<- subset_matrices[[1]][,c(grep(pattern = "healthy", disease_GSE199949))]
healthy_subsets[["GSE199949_symbol_expression_matrix_healthy_subset"]] <- healthy_GSE199949_subset

#############################################################################################################################
meta_GSE213001<-read.xlsx(metadata_paths[[2]], sheetIndex = 1)
meta_GSE213001<-meta_GSE213001[meta_GSE213001$disease%in%c("IPF", "healthy"),]
disease_GSE213001<-meta_GSE213001$disease

disease_GSE213001_subset<- subset_matrices[[2]][c(grep(pattern = "IPF", disease_GSE213001))]
disease_subsets[["GSE213001_symbol_expression_matrix_disease_subset"]] <- disease_GSE213001_subset
healthy_GSE213001_subset<-  subset_matrices[[2]][,c(grep(pattern = "healthy", disease_GSE213001))]
healthy_subsets[["GSE213001_symbol_expression_matrix_healthy_subset"]] <- healthy_GSE213001_subset

#################################################################################################

meta_GSE92592<-read.xlsx(metadata_paths[[10]], sheetIndex = 1)
disease_GSE92592<-meta_GSE92592$disease
disease_GSE92592_subset<- subset_matrices[[10]][,c(grep(pattern = "IPF", disease_GSE92592))]
disease_subsets[["GSE92592_symbol_expression_matrix_disease_subset"]] <- disease_GSE92592_subset
healthy_GSE92592_subset<-  subset_matrices[[10]][,c(grep(pattern = "healthy", disease_GSE92592))]
healthy_subsets[["GSE92592_symbol_expression_matrix_healthy_subset"]] <- healthy_GSE92592_subset

###################################################################################################33
meta_GSE99621<-read.xlsx(metadata_paths[[11]], sheetIndex = 1)
disease_GSE99621<-meta_GSE99621$disease

disease_GSE99621_subset<-subset_matrices[[11]][c(grep(pattern = "IPF", colnames(subset_matrices[[11]])))]
disease_subsets[["GSE99621_symbol_expression_matrix_disease_subset"]] <- disease_GSE99621_subset
healthy_GSE99621_subset<- subset_matrices[[11]][,c(grep(pattern = "HC", colnames(subset_matrices[[11]])))]
healthy_subsets[["GSE99621_symbol_expression_matrix_healthy_subset"]] <- healthy_GSE99621_subset

############################################################################################

meta_GSE138283<-read.xlsx(metadata_paths[[5]], sheetIndex = 1)
disease_GSE138283<-meta_GSE138283$disease

disease_GSE138283_subset<-subset_matrices[[5]][c(grep(pattern = "IPF", disease_GSE138283))]
disease_subsets[["GSE138283_symbol_expression_matrix_disease_subset"]] <- disease_GSE138283_subset
healthy_GSE138283_subset<- subset_matrices[[5]][,c(grep(pattern = "healthy", disease_GSE138283))]
healthy_subsets[["GSE138283_symbol_expression_matrix_healthy_subset"]] <- healthy_GSE138283_subset




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

setwd("/Lung_data/Biopsy/RNA-seq/symbol_norm/subsets/disease")

for (i in 1:length(disease_subsets)) {
  write.table(disease_subsets_ordered[[i]], file=paste0(names(disease_subsets_ordered[i]),".txt"), sep="\t")
}


setwd("/Lung_data/Biopsy/RNA-seq/symbol_norm/subsets/healthy")

for (i in 1:length(healthy_subsets_ordered)) {
  write.table(healthy_subsets_ordered[[i]], file=paste0(names(healthy_subsets_ordered[i]),".txt"), sep="\t")
}


setwd("/Lung_data/Biopsy/RNA-seq/symbol_norm/subsets/disease")

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

setwd("/nasdata/sinkala/Lung_data/Biopsy/RNA-seq/symbol_norm/subsets/healthy")

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


