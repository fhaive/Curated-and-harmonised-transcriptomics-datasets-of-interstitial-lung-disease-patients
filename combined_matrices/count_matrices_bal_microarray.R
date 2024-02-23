#########make the consensus count matrices############
#These paths in this bunch of code doesnt exist in this folder. They are psudopaths
#symbol_norm folder contains all the microarray bal expression matrices and all the other 
#folders in this script are under this folder

setwd("/Lung_data/BAL/Microarray/symbol_norm")


GSE70866 <- read.table("expression_matrix_symbol_GSE70866.csv", sep="\t")
library(xlsx)
metadata_GSE70866<-read.xlsx("/Lung_data/BAL/Microarray/GSE70866/phenodata/GSE70866_curated.xlsx", sheetIndex = 1)
metadata_GSE70866<-metadata_GSE70866[metadata_GSE70866$cohort==c("Freiburg"),]

#####In expression matrix already just Freiburg samples included 
####Disease and healthy

disease_GSE70866_subset<-GSE70866[,metadata_GSE70866$disease=="IPF"]

healthy_GSE70866_subset<-GSE70866[,metadata_GSE70866$disease=="healthy"]

setwd("/nasdata/sinkala/Lung_data/BAL/Microarray/symbol_norm/disease")
write.table(disease_GSE70866_subset, file="GSE70866_symbol_expression_matrix_disease_subset.txt", sep="\t")

setwd("/nasdata/sinkala/Lung_data/BAL/Microarray/symbol_norm/healthy")
write.table(healthy_GSE70866_subset, file="GSE70866_symbol_expression_matrix_healthy_subset.txt", sep="\t")

