###########Microarraydatasets_Ensembl to geneid################

setwd("/Lung_data/Fibroblast/Microarray/GSE144338/expression_data/expression_matrices")

library(xlsx)


expression_matrix<-read.table("RAW_Expression_Matrix_Normalized_2023-10-27.txt", sep="\t", header=T, row.names = NULL)

###########Expression_matrix ensg#######################

rownames_table<-c()
for (i in 1:(length(expression_matrix$row.names))) {
  rowname<-strsplit(as.character(expression_matrix$row.names[i]), "_")[[1]][1]
  rownames_table<-c(rownames_table, rowname)
}

expression_matrix$ensg<-rownames_table

expression_matrix<-expression_matrix[-grep("AFFX", expression_matrix$ensg),]

rownames(expression_matrix)<-expression_matrix$ensg

expression_matrix<-expression_matrix[,2:(length(colnames(expression_matrix))-1)]

write.table(expression_matrix, file = "expression_matrix_ensembl_GSE144338.csv", sep = "\t", col.names =TRUE)

###########Expression_matrix symbol#######################

require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

annotLookup <- getBM(
  mart=mart,
  attributes=c(
    "ensembl_gene_id",
    "gene_biotype",
    "external_gene_name"),
  filter = "ensembl_gene_id",
  values = rownames(expression_matrix),
  uniqueRows = TRUE)

expression_matrix$ensg<-rownames(expression_matrix)

table_annot <- merge(expression_matrix, annotLookup, by.x="ensg", by.y="ensembl_gene_id")


table_annot_1<-table_annot[,c(2:13,15)]
aggr_table<- table_annot_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(expression_matrix)[c(1:(length(colnames(expression_matrix))-1))]

write.table(aggr_table, file = "expression_matrix_symbol_GSE144338.csv", sep = "\t", col.names =TRUE)
###############################################################################################################

##########################DEG_analysis###########################################
#########################Ensembl_IDS##############################################
############################IPF####################################################


setwd("/Lung_data/Fibroblast/Microarray/GSE144338/expression_data/DEG_results")

dif_table_ipf<- as.data.frame(readxl::read_excel("ALL_Differential_Expression_Tables_2023-11-15.xlsx", sheet = 2))

dif_table_adeno_healthy<- as.data.frame(readxl::read_excel("ALL_Differential_Expression_Tables_2023-11-15.xlsx", sheet = 3))

dif_table_ipf_adeno<-as.data.frame(readxl::read_excel("ALL_Differential_Expression_Tables_2023-11-15.xlsx", sheet = 4))


rownames_table<-c()
for (i in 1:(length(dif_table_ipf$ID))) {
  rowname<-strsplit(as.character(dif_table_ipf$ID[i]), "_")[[1]][1]
  rownames_table<-c(rownames_table, rowname)
}

dif_table_ipf$ensg<-rownames_table

dif_table_ipf<-dif_table_ipf[-grep("AFFX", dif_table_ipf$ensg),]

rownames(dif_table_ipf)<-dif_table_ipf$ensg

dif_table_ipf<-dif_table_ipf[,1:(length(colnames(dif_table_ipf))-2)]

dif_table_ipf_filtered<-dif_table_ipf[which(dif_table_ipf$adj.P.Val<=0.01 & abs(dif_table_ipf$logFC)>=0.58),]

write.table(dif_table_ipf_filtered, file = "DEG_results_GSE144338_ipf_vs_healthy_filtered_ensembl.csv", sep = "\t", col.names =TRUE)

write.table(dif_table_ipf, file = "DEG_results_GSE144338_ipf_vs_healthy_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)


#########Gene symbols##############################

require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

annotLookup <- getBM(
  mart=mart,
  attributes=c(
    "ensembl_gene_id",
    "gene_biotype",
    "external_gene_name"),
  filter = "ensembl_gene_id",
  values = rownames(dif_table_ipf),
  uniqueRows = TRUE)

dif_table_ipf$ensg<-rownames(dif_table_ipf)

table_annot<- merge(dif_table_ipf, annotLookup, by.x="ensg", by.y="ensembl_gene_id")


table_annot_1<-table_annot[,c(2:8,10)]
aggr_table<- table_annot_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(dif_table_ipf)[c(1:(length(colnames(dif_table_ipf))-1))]

dif_table_ipf_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_ipf_filtered, file = "DEG_results_GSE144338_ipf_vs_healthy_filtered_symbol.csv", sep = "\t", col.names =TRUE)


write.table(aggr_table, file = "DEG_results_GSE144338_ipf_vs_healthy_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

##########################DEG_analysis###########################################
#########################Ensembl_IDS##############################################
############################ADENOCARCINOMA#########################################

rownames_table<-c()
for (i in 1:(length(dif_table_adeno_healthy$ID))) {
  rowname<-strsplit(as.character(dif_table_adeno_healthy$ID[i]), "_")[[1]][1]
  rownames_table<-c(rownames_table, rowname)
}

dif_table_adeno_healthy$ensg<-rownames_table

dif_table_adeno_healthy<-dif_table_adeno_healthy[-grep("AFFX", dif_table_adeno_healthy$ensg),]

rownames(dif_table_adeno_healthy)<-dif_table_adeno_healthy$ensg

dif_table_adeno_healthy<-dif_table_adeno_healthy[,1:(length(colnames(dif_table_adeno_healthy))-2)]

dif_table_adeno_healthy_filtered<-dif_table_adeno_healthy[which(dif_table_adeno_healthy$adj.P.Val<=0.01 & abs(dif_table_adeno_healthy$logFC)>=0.58),]

write.table(dif_table_adeno_healthy_filtered, file = "DEG_results_GSE144338_adeno_vs_healthy_filtered_ensembl.csv", sep = "\t", col.names =TRUE)

write.table(dif_table_adeno_healthy, file = "DEG_results_GSE144338_adeno_vs_healthy_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)


#########Gene symbols##############################

dif_table_adeno_healthy$ensg<-rownames(dif_table_adeno_healthy)

table_annot<- merge(dif_table_adeno_healthy, annotLookup, by.x="ensg", by.y="ensembl_gene_id")


table_annot_1<-table_annot[,c(2:8,10)]
aggr_table<- table_annot_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(dif_table_adeno_healthy)[c(1:(length(colnames(dif_table_adeno_healthy))-1))]

dif_table_adeno_healthy_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_adeno_healthy_filtered, file = "DEG_results_GSE144338_adeno_vs_healthy_filtered_symbol.csv", sep = "\t", col.names =TRUE)

write.table(aggr_table, file = "DEG_results_GSE144338_adeno_vs_healthy_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

##########################DEG_analysis###########################################
#########################Ensembl_IDS##############################################
############################IPF_VS_ADENOCARCINOMA#########################################

rownames_table<-c()
for (i in 1:(length(dif_table_ipf_adeno$ID))) {
  rowname<-strsplit(as.character(dif_table_ipf_adeno$ID[i]), "_")[[1]][1]
  rownames_table<-c(rownames_table, rowname)
}

dif_table_ipf_adeno$ensg<-rownames_table

dif_table_ipf_adeno<-dif_table_ipf_adeno[-grep("AFFX", dif_table_ipf_adeno$ensg),]

rownames(dif_table_ipf_adeno)<-dif_table_ipf_adeno$ensg

dif_table_ipf_adeno<-dif_table_ipf_adeno[,1:(length(colnames(dif_table_ipf_adeno))-2)]

dif_table_ipf_adeno_filtered<-dif_table_ipf_adeno[which(dif_table_ipf_adeno$adj.P.Val<=0.01 & abs(dif_table_ipf_adeno$logFC)>=0.58),]

write.table(dif_table_ipf_adeno_filtered, file = "DEG_results_GSE144338_ipf_vs_adeno_filtered_ensembl.csv", sep = "\t", col.names =TRUE)

write.table(dif_table_ipf_adeno, file = "DEG_results_GSE144338_ipf_vs_adeno_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)


#########Gene symbols##############################

dif_table_ipf_adeno$ensg<-rownames(dif_table_ipf_adeno)

table_annot<- merge(dif_table_ipf_adeno, annotLookup, by.x="ensg", by.y="ensembl_gene_id")


table_annot_1<-table_annot[,c(2:8,10)]
aggr_table<- table_annot_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(dif_table_ipf_adeno)[c(1:(length(colnames(dif_table_ipf_adeno))-1))]

dif_table_ipf_adeno_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_ipf_adeno_filtered, file = "DEG_results_GSE144338_ipf_vs_adeno_filtered_symbol.csv", sep = "\t", col.names =TRUE)

write.table(aggr_table, file = "DEG_results_GSE144338_ipf_vs_adeno_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

