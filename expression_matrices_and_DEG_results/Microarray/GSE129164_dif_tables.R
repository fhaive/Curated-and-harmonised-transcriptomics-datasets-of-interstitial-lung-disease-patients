###########Expression_matrix ensg#######################
setwd("/Lung_data/Fibroblast/Microarray/GSE129164/expression_data/expression_matrices")

library(xlsx)


expression_matrix<-read.table("RAW_Expression_Matrix_Normalized_2023-10-27.txt", sep="\t", header=T, row.names = NULL)



rownames_table<-c()
for (i in 1:(length(expression_matrix$row.names))) {
  rowname<-strsplit(as.character(expression_matrix$row.names[i]), "_")[[1]][1]
  rownames_table<-c(rownames_table, rowname)
}

expression_matrix$ensg<-rownames_table

expression_matrix<-expression_matrix[-grep("AFFX", expression_matrix$ensg),]

rownames(expression_matrix)<-expression_matrix$ensg

expression_matrix<-expression_matrix[,2:(length(colnames(expression_matrix))-1)]

write.table(expression_matrix, file = "expression_matrix_ensembl_GSE129164.csv", sep = "\t", col.names =TRUE)

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


table_annot_1<-table_annot[,c(2:21,23)]
aggr_table<- table_annot_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(expression_matrix)[c(1:(length(colnames(expression_matrix))-1))]

write.table(aggr_table, file = "expression_matrix_symbol_GSE129164.csv", sep = "\t", col.names =TRUE)


###########Microarray datasets: Ensembl to geneid################

setwd("/Lung_data/Fibroblast/Microarray/GSE129164/expression_data/DEG_results")


dif_table_ipf_healthy <- as.data.frame(readxl::read_excel("ALL_IPF_healthy_Differential_Expression_Tables_2023-10-27.xlsx", sheet= 2))

dif_table_non_stim_ipf_non_stim_health<-as.data.frame(readxl::read_excel("ALL_treatments_Differential_Expression_Tables_2023-10-27.xlsx", sheet = 2))

dif_table_tgf_ipf_tgf_heal<-as.data.frame(readxl::read_excel("ALL_treatments_Differential_Expression_Tables_2023-10-27.xlsx", sheet = 3))

dif_table_tgf_ipf_non_stim_ipf<-as.data.frame(readxl::read_excel("ALL_treatments_Differential_Expression_Tables_2023-10-27.xlsx", sheet = 4))

dif_table_tgf_heal_non_stim_heal<-as.data.frame(readxl::read_excel("ALL_treatments_Differential_Expression_Tables_2023-10-27.xlsx", sheet = 5))


#########################Ensembl_IDS##############################################
rownames_table<-c()
for (i in 1:(length(dif_table_ipf_healthy$ID))) {
  rowname<-strsplit(as.character(dif_table_ipf_healthy$ID[i]), "_")[[1]][1]
  rownames_table<-c(rownames_table, rowname)
}

dif_table_ipf_healthy$ensg<-rownames_table

dif_table_ipf_healthy<-dif_table_ipf_healthy[-grep("AFFX", dif_table_ipf_healthy$ensg),]

rownames(dif_table_ipf_healthy)<-dif_table_ipf_healthy$ensg

dif_table_ipf_healthy<-dif_table_ipf_healthy[,1:(length(colnames(dif_table_ipf_healthy))-2)]

dif_table_ipf_healthy_filtered<-dif_table_ipf_healthy[which(dif_table_ipf_healthy$adj.P.Val<=0.01 & abs(dif_table_ipf_healthy$logFC)>=0.58),]

write.table(dif_table_ipf_healthy_filtered, file = "DEG_results_GSE129164_all_ipf_vs_all_healthy_filtered_ensembl.csv", sep = "\t", col.names =TRUE)

write.table(dif_table_ipf_healthy, file = "DEG_results_GSE129164_all_ipf_vs_all_healthy_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)


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
  values = rownames(dif_table_ipf_healthy),
  uniqueRows = TRUE)

dif_table_ipf_healthy$ensg<-rownames(dif_table_ipf_healthy)

table_annot<- merge(dif_table_ipf_healthy, annotLookup, by.x="ensg", by.y="ensembl_gene_id")


table_annot_1<-table_annot[,c(2:8,10)]
aggr_table<- table_annot_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(dif_table_ipf_healthy)[c(1:(length(colnames(dif_table_ipf_healthy))-1))]

dif_table_ipf_healthy_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_ipf_healthy_filtered, file = "DEG_results_GSE129164_all_ipf_vs_all_healthy_filtered_symbol.csv", sep = "\t", col.names =TRUE)

write.table(aggr_table, file = "DEG_results_GSE129164_all_ipf_vs_all_healthy_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)


###############################################################################################################
#########################Ensembl_IDS-Nonstim_ipf_non_stim_heal##############################################
rownames_table<-c()
for (i in 1:(length(dif_table_non_stim_ipf_non_stim_health$ID))) {
  rowname<-strsplit(as.character(dif_table_non_stim_ipf_non_stim_health$ID[i]), "_")[[1]][1]
  rownames_table<-c(rownames_table, rowname)
}

dif_table_non_stim_ipf_non_stim_health$ensg<-rownames_table

dif_table_non_stim_ipf_non_stim_health<-dif_table_non_stim_ipf_non_stim_health[-grep("AFFX", dif_table_non_stim_ipf_non_stim_health$ensg),]

rownames(dif_table_non_stim_ipf_non_stim_health)<-dif_table_non_stim_ipf_non_stim_health$ensg

dif_table_non_stim_ipf_non_stim_health<-dif_table_non_stim_ipf_non_stim_health[,1:(length(colnames(dif_table_non_stim_ipf_non_stim_health))-2)]

dif_table_non_stim_ipf_non_stim_health_filtered<-dif_table_non_stim_ipf_non_stim_health[which(dif_table_non_stim_ipf_non_stim_health$adj.P.Val<=0.01 & abs(dif_table_non_stim_ipf_non_stim_health$logFC)>=0.58),]

write.table(dif_table_non_stim_ipf_non_stim_health_filtered, file = "DEG_results_GSE129164_non_stim_ipf_vs_non_stim_healthy_filtered_ensembl.csv", sep = "\t", col.names =TRUE)


write.table(dif_table_non_stim_ipf_non_stim_health, file = "DEG_results_GSE129164_non_stim_ipf_vs_non_stim_healthy_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)


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
  values = rownames(dif_table_non_stim_ipf_non_stim_health),
  uniqueRows = TRUE)

dif_table_non_stim_ipf_non_stim_health$ensg<-rownames(dif_table_non_stim_ipf_non_stim_health)

table_annot<- merge(dif_table_non_stim_ipf_non_stim_health, annotLookup, by.x="ensg", by.y="ensembl_gene_id")


table_annot_1<-table_annot[,c(2:8,10)]
aggr_table<- table_annot_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(dif_table_non_stim_ipf_non_stim_health)[c(1:(length(colnames(dif_table_non_stim_ipf_non_stim_health))-1))]

dif_table_non_stim_ipf_non_stim_health_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_non_stim_ipf_non_stim_health_filtered, file = "DEG_results_GSE129164_non_stim_ipf_vs_non_stim_healthy_filtered_symbol.csv", sep = "\t", col.names =TRUE)

write.table(aggr_table, file = "DEG_results_GSE129164_non_stim_ipf_vs_non_stim_healthy_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

#########################Ensembl_IDS-dif_table_tgf_ipf_tgf_heal##############################################
rownames_table<-c()
for (i in 1:(length(dif_table_tgf_ipf_tgf_heal$ID))) {
  rowname<-strsplit(as.character(dif_table_tgf_ipf_tgf_heal$ID[i]), "_")[[1]][1]
  rownames_table<-c(rownames_table, rowname)
}

dif_table_tgf_ipf_tgf_heal$ensg<-rownames_table

dif_table_tgf_ipf_tgf_heal<-dif_table_tgf_ipf_tgf_heal[-grep("AFFX", dif_table_tgf_ipf_tgf_heal$ensg),]

rownames(dif_table_tgf_ipf_tgf_heal)<-dif_table_tgf_ipf_tgf_heal$ensg

dif_table_tgf_ipf_tgf_heal<-dif_table_tgf_ipf_tgf_heal[,1:(length(colnames(dif_table_tgf_ipf_tgf_heal))-2)]

dif_table_tgf_ipf_tgf_heal_filtered<-dif_table_tgf_ipf_tgf_heal[which(dif_table_tgf_ipf_tgf_heal$adj.P.Val<=0.01 & abs(dif_table_tgf_ipf_tgf_heal$logFC)>=0.58),]

write.table(dif_table_tgf_ipf_tgf_heal_filtered, file = "DEG_results_GSE129164_tgf_ipf_vs_tgf_healthy_filtered_ensembl.csv", sep = "\t", col.names =TRUE)


write.table(dif_table_tgf_ipf_tgf_heal, file = "DEG_results_GSE129164_tgf_ipf_vs_tgf_healthy_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)


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
  values = rownames(dif_table_tgf_ipf_tgf_heal),
  uniqueRows = TRUE)

dif_table_tgf_ipf_tgf_heal$ensg<-rownames(dif_table_tgf_ipf_tgf_heal)

table_annot<- merge(dif_table_tgf_ipf_tgf_heal, annotLookup, by.x="ensg", by.y="ensembl_gene_id")


table_annot_1<-table_annot[,c(2:8,10)]
aggr_table<- table_annot_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(dif_table_tgf_ipf_tgf_heal)[c(1:(length(colnames(dif_table_tgf_ipf_tgf_heal))-1))]

dif_table_tgf_ipf_tgf_heal_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_tgf_ipf_tgf_heal_filtered, file = "DEG_results_GSE129164_tgf_ipf_vs_tgf_healthy_filtered_symbol.csv", sep = "\t", col.names =TRUE)

write.table(aggr_table, file = "DEG_results_GSE129164_tgf_ipf_vs_tgf_healthy_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

#########################Ensembl_IDS-dif_table_dif_table_tgf_ipf_non_stim_ipf##############################################
rownames_table<-c()
for (i in 1:(length(dif_table_tgf_ipf_non_stim_ipf$ID))) {
  rowname<-strsplit(as.character(dif_table_tgf_ipf_non_stim_ipf$ID[i]), "_")[[1]][1]
  rownames_table<-c(rownames_table, rowname)
}

dif_table_tgf_ipf_non_stim_ipf$ensg<-rownames_table

dif_table_tgf_ipf_non_stim_ipf<-dif_table_tgf_ipf_non_stim_ipf[-grep("AFFX", dif_table_tgf_ipf_non_stim_ipf$ensg),]

rownames(dif_table_tgf_ipf_non_stim_ipf)<-dif_table_tgf_ipf_non_stim_ipf$ensg

dif_table_tgf_ipf_non_stim_ipf<-dif_table_tgf_ipf_non_stim_ipf[,1:(length(colnames(dif_table_tgf_ipf_non_stim_ipf))-2)]

dif_table_tgf_ipf_non_stim_ipf_filtered<-dif_table_tgf_ipf_non_stim_ipf[which(dif_table_tgf_ipf_non_stim_ipf$adj.P.Val<=0.01 & abs(dif_table_tgf_ipf_non_stim_ipf$logFC)>=0.58),]

write.table(dif_table_tgf_ipf_non_stim_ipf_filtered, file = "DEG_results_GSE129164_tgf_ipf_vs_non_stim_ipf_filtered_ensembl.csv", sep = "\t", col.names =TRUE)

write.table(dif_table_tgf_ipf_non_stim_ipf, file = "DEG_results_GSE129164_tgf_ipf_vs_non_stim_ipf_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)


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
  values = rownames(dif_table_tgf_ipf_non_stim_ipf),
  uniqueRows = TRUE)

dif_table_tgf_ipf_non_stim_ipf$ensg<-rownames(dif_table_tgf_ipf_non_stim_ipf)

table_annot<- merge(dif_table_tgf_ipf_non_stim_ipf, annotLookup, by.x="ensg", by.y="ensembl_gene_id")


table_annot_1<-table_annot[,c(2:8,10)]
aggr_table<- table_annot_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(dif_table_tgf_ipf_non_stim_ipf)[c(1:(length(colnames(dif_table_tgf_ipf_non_stim_ipf))-1))]

dif_table_tgf_ipf_non_stim_ipf_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_tgf_ipf_non_stim_ipf_filtered, file = "DEG_results_GSE129164_tgf_ipf_vs_non_stim_ipf_filtered_symbol.csv", sep = "\t", col.names =TRUE)

write.table(aggr_table, file = "DEG_results_GSE129164_tgf_ipf_vs_non_stim_ipf_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

#########################Ensembl_IDS-dif_table_tgf_heal_non_stim_heal##############################################
rownames_table<-c()
for (i in 1:(length(dif_table_tgf_heal_non_stim_heal$ID))) {
  rowname<-strsplit(as.character(dif_table_tgf_heal_non_stim_heal$ID[i]), "_")[[1]][1]
  rownames_table<-c(rownames_table, rowname)
}

dif_table_tgf_heal_non_stim_heal$ensg<-rownames_table

dif_table_tgf_heal_non_stim_heal<-dif_table_tgf_heal_non_stim_heal[-grep("AFFX", dif_table_tgf_heal_non_stim_heal$ensg),]

rownames(dif_table_tgf_heal_non_stim_heal)<-dif_table_tgf_heal_non_stim_heal$ensg

dif_table_tgf_heal_non_stim_heal<-dif_table_tgf_heal_non_stim_heal[,1:(length(colnames(dif_table_tgf_heal_non_stim_heal))-2)]

dif_table_tgf_heal_non_stim_heal_filtered<-dif_table_tgf_heal_non_stim_heal[which(dif_table_tgf_heal_non_stim_heal$adj.P.Val<=0.01 & abs(dif_table_tgf_heal_non_stim_heal$logFC)>=0.58),]

write.table(dif_table_tgf_heal_non_stim_heal_filtered, file = "DEG_results_GSE129164_tgf_heal_vs_non_stim_healthy_filtered_ensembl.csv", sep = "\t", col.names =TRUE)

write.table(dif_table_tgf_heal_non_stim_heal, file = "DEG_results_GSE129164_tgf_healthy_vs_non_stim_healthy_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)

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
  values = rownames(dif_table_tgf_heal_non_stim_heal),
  uniqueRows = TRUE)

dif_table_tgf_heal_non_stim_heal$ensg<-rownames(dif_table_tgf_heal_non_stim_heal)

table_annot<- merge(dif_table_tgf_heal_non_stim_heal, annotLookup, by.x="ensg", by.y="ensembl_gene_id")


table_annot_1<-table_annot[,c(2:8,10)]
aggr_table<- table_annot_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(dif_table_tgf_heal_non_stim_heal)[c(1:(length(colnames(dif_table_tgf_heal_non_stim_heal))-1))]

dif_table_tgf_heal_non_stim_heal_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]



write.table(dif_table_tgf_heal_non_stim_heal_filtered, file = "DEG_results_GSE129164_tgf_healthy_vs_non_stim_healthy_filtered_symbol.csv", sep = "\t", col.names =TRUE)


write.table(aggr_table, file = "DEG_results_GSE129164_tgf_healthy_vs_non_stim_healthy_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)
