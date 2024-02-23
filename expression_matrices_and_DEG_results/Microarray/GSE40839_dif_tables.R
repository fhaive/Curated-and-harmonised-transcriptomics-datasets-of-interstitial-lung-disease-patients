###########Microarraydatasets_Ensembl to geneid################

setwd("/Lung_data/Fibroblast/Microarray/GSE40839/expression_data/DEG_results")

library(xlsx)


dif_table_scild <- read.xlsx("ALL_diseases_Differential_Expression_Tables_2023-10-27.xlsx", sheetIndex = 3)

dif_table_uip <- read.xlsx("ALL_diseases_Differential_Expression_Tables_2023-10-27.xlsx", sheetIndex = 2)

dif_table_uip<-dif_table_uip[order(dif_table_uip$ID),]

dif_table_scild<-dif_table_scild[order(dif_table_scild$ID),]

#identical(dif_table_scild$ID, dif_table_uip$ID)
#[1] TRUE

require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

annotLookup <- getBM(
  mart=mart,
  attributes=c(
    "affy_hg_u133a_2",
    "ensembl_gene_id",
    "gene_biotype",
    "external_gene_name"),
  filter = "affy_hg_u133a_2",
  values = dif_table_scild$ID,
  uniqueRows = TRUE)

############Scild_dif_tables###############################
table_with_gene_id <- merge(dif_table_scild, annotLookup, by.x="ID", by.y="affy_hg_u133a_2")

table_with_gene_id_1<-table_with_gene_id[,c(2:(length(colnames(table_with_gene_id))-2))] # REMOVE excess columns
aggr_table<- table_with_gene_id_1 %>% group_by(ensembl_gene_id) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
#aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$ensembl_gene_id
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(dif_table_scild)[c(1:(length(colnames(dif_table_scild))-1))]

dif_table_scild_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_scild_filtered, file = "DEG_results_GSE40839_scild_vs_healthy_filtered_ensembl.csv", sep = "\t", col.names =TRUE)

write.table(aggr_table, file = "DEG_results_GSE40839_scild_vs_healthy_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)

table_with_gene_id_1<-table_with_gene_id[,c(2:(length(colnames(table_with_gene_id))-3),11)] # REMOVE excess columns
aggr_table<- table_with_gene_id_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(dif_table_scild)[c(1:(length(colnames(dif_table_scild))-1))]

dif_table_scild_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_scild_filtered, file = "DEG_results_GSE40839_scild_vs_healthy_filtered_symbol.csv", sep = "\t", col.names =TRUE)


write.table(aggr_table, file = "DEG_results_GSE40839_scild_vs_healthy_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

############UIP_dif_tables###############################
table_with_gene_id <- merge(dif_table_uip, annotLookup, by.x="ID", by.y="affy_hg_u133a_2")

table_with_gene_id_1<-table_with_gene_id[,c(2:(length(colnames(table_with_gene_id))-2))] # REMOVE excess columns
aggr_table<- table_with_gene_id_1 %>% group_by(ensembl_gene_id) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
#aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$ensembl_gene_id
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(dif_table_scild)[c(1:(length(colnames(dif_table_scild))-1))]

dif_table_ipf_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_ipf_filtered, file = "DEG_results_GSE40839_uip_vs_healthy_filtered_ensembl.csv", sep = "\t", col.names =TRUE)

write.table(aggr_table, file = "DEG_results_GSE40839_uip_vs_healthy_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)

table_with_gene_id_1<-table_with_gene_id[,c(2:(length(colnames(table_with_gene_id))-3),11)] # REMOVE excess columns
aggr_table<- table_with_gene_id_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(dif_table_scild)[c(1:(length(colnames(dif_table_scild))-1))]

dif_table_ipf_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_ipf_filtered, file = "DEG_results_GSE40839_uip_vs_healthy_filtered_symbol.csv", sep = "\t", col.names =TRUE)


write.table(aggr_table, file = "DEG_results_GSE40839_uip_vs_healthy_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)


###############################Expression_matrix#################################################
############Biomart Ensg

setwd("/Lung_data/Fibroblast/Microarray/GSE40839/expression_data/expression_matrices")

raw_matrix <- read.csv("RAW_Expression_Matrix_Normalized_2023-10-27.txt", sep="\t")

metadata <- read.xlsx("/Lung_data/Fibroblast/Microarray/GSE40839/phenodata/GSE40839_curated.xlsx", sheetIndex = 1)

require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

annotLookup <- getBM(
  mart=mart,
  attributes=c(
    "affy_hg_u133a_2",
    "ensembl_gene_id",
    "gene_biotype",
    "external_gene_name"),
  filter = "affy_hg_u133a_2",
  values = rownames(raw_matrix),
  uniqueRows = TRUE)

raw_matrix$probes<-rownames(raw_matrix)

table_with_gene_id <- merge(raw_matrix, annotLookup, by.x="probes", by.y="affy_hg_u133a_2")

table_with_gene_id_1<-table_with_gene_id[,c(2:(length(colnames(table_with_gene_id))-2))] # REMOVE excess columns
aggr_table<- table_with_gene_id_1 %>% group_by(ensembl_gene_id) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
#aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$ensembl_gene_id
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(raw_matrix)[c(1:(length(colnames(raw_matrix))-1))]
write.table(aggr_table, file = "expression_matrix_ensembl_GSE40839.csv", sep = "\t", col.names =TRUE)

table_with_gene_id_1<-table_with_gene_id[,c(2:(length(colnames(table_with_gene_id))-3),25)] # REMOVE excess columns
aggr_table<- table_with_gene_id_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(raw_matrix)[c(1:(length(colnames(raw_matrix))-1))]
write.table(aggr_table, file = "expression_matrix_symbol_GSE40839.csv", sep = "\t", col.names =TRUE)
