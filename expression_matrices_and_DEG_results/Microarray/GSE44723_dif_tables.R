###########Microarraydatasets_Ensembl to geneid################


setwd("/Lung_data/Fibroblast/Microarray/GSE44723/expression_data/expression_matrices")

library(xlsx)

expression_matrix<-read.table("RAW_Expression_Matrix_Normalized_2023-01-24.txt", sep="\t", header=T, row.names = NULL)

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

write.table(expression_matrix, file = "expression_matrix_ensembl_GSE44723.csv", sep = "\t", col.names =TRUE)

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

table_annot_ipf <- merge(expression_matrix, annotLookup, by.x="ensg", by.y="ensembl_gene_id")


table_annot_ipf_1<-table_annot_ipf[,c(2:15,17)]
aggr_table<- table_annot_ipf_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(expression_matrix)[c(1:(length(colnames(expression_matrix))-1))]

write.table(aggr_table, file = "expression_matrix_symbol_GSE44723.csv", sep = "\t", col.names =TRUE)
###############################################################################################################

setwd("/Lung_data/Fibroblast/Microarray/GSE44723/expression_data/DEG_results")

dif_table_ipf <- as.data.frame(readxl::read_excel("ALL_Differential_Expression_Tables_2023-01-24.xlsx", sheet = 2))
dif_table_rapid_healthy<-as.data.frame(readxl::read_excel("ALL_disease_state_Differential_Expression_Tables_2023-11-15.xlsx", sheet = 2))
dif_table_slow_healthy<-as.data.frame(readxl::read_excel("ALL_disease_state_Differential_Expression_Tables_2023-11-15.xlsx", sheet = 3))
dif_table_rapid_slow<-as.data.frame(readxl::read_excel("ALL_disease_state_Differential_Expression_Tables_2023-11-15.xlsx", sheet = 4))


###################DEG
#######ipf
#########################Ensembl_IDS##############################################
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

write.table(dif_table_ipf_filtered, file = "DEG_results_GSE44723_ipf_vs_healthy_filtered_ensembl.csv", sep = "\t", col.names =TRUE)


write.table(dif_table_ipf, file = "DEG_results_GSE44723_ipf_vs_healthy_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)


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

table_annot_ipf <- merge(dif_table_ipf, annotLookup, by.x="ensg", by.y="ensembl_gene_id")


table_annot_ipf_1<-table_annot_ipf[,c(2:8,10)]
aggr_table<- table_annot_ipf_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(dif_table_ipf)[c(1:(length(colnames(dif_table_ipf))-1))]

dif_table_ipf_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_ipf_filtered, file = "DEG_results_GSE44723_ipf_vs_healthy_filtered_symbol.csv", sep = "\t", col.names =TRUE)

write.table(aggr_table, file = "DEG_results_GSE44723_ipf_vs_healthy_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

###################DEG
#######dif_table_rapid_healthy
#########################Ensembl_IDS##############################################
rownames_table<-c()
for (i in 1:(length(dif_table_rapid_healthy$ID))) {
  rowname<-strsplit(as.character(dif_table_rapid_healthy$ID[i]), "_")[[1]][1]
  rownames_table<-c(rownames_table, rowname)
}

dif_table_rapid_healthy$ensg<-rownames_table

dif_table_rapid_healthy<-dif_table_rapid_healthy[-grep("AFFX", dif_table_rapid_healthy$ensg),]

rownames(dif_table_rapid_healthy)<-dif_table_rapid_healthy$ensg

dif_table_rapid_healthy<-dif_table_rapid_healthy[,1:(length(colnames(dif_table_rapid_healthy))-2)]

dif_table_rapid_healthy_filtered<-dif_table_rapid_healthy[which(dif_table_rapid_healthy$adj.P.Val<=0.01 & abs(dif_table_rapid_healthy$logFC)>=0.58),]

write.table(dif_table_rapid_healthy_filtered, file = "DEG_results_GSE44723_rapid_vs_healthy_filtered_ensembl.csv", sep = "\t", col.names =TRUE)


write.table(dif_table_rapid_healthy, file = "DEG_results_GSE44723_rapid_vs_healthy_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)


#########Gene symbols##############################


dif_table_rapid_healthy$ensg<-rownames(dif_table_rapid_healthy)

table_annot_ipf <- merge(dif_table_rapid_healthy, annotLookup, by.x="ensg", by.y="ensembl_gene_id")


table_annot_ipf_1<-table_annot_ipf[,c(2:8,10)]
aggr_table<- table_annot_ipf_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(dif_table_rapid_healthy)[c(1:(length(colnames(dif_table_rapid_healthy))-1))]

dif_table_rapid_healthy_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_rapid_healthy_filtered, file = "DEG_results_GSE44723_rapid_vs_healthy_filtered_symbol.csv", sep = "\t", col.names =TRUE)

write.table(aggr_table, file = "DEG_results_GSE44723_rapid_vs_healthy_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

###################DEG
#######dif_table_rapid_slow
#########################Ensembl_IDS##############################################
rownames_table<-c()
for (i in 1:(length(dif_table_rapid_slow$ID))) {
  rowname<-strsplit(as.character(dif_table_rapid_slow$ID[i]), "_")[[1]][1]
  rownames_table<-c(rownames_table, rowname)
}

dif_table_rapid_slow$ensg<-rownames_table

dif_table_rapid_slow<-dif_table_rapid_slow[-grep("AFFX", dif_table_rapid_slow$ensg),]

rownames(dif_table_rapid_slow)<-dif_table_rapid_slow$ensg

dif_table_rapid_slow<-dif_table_rapid_slow[,1:(length(colnames(dif_table_rapid_slow))-2)]

dif_table_rapid_slow_filtered<-dif_table_rapid_slow[which(dif_table_rapid_slow$adj.P.Val<=0.01 & abs(dif_table_rapid_slow$logFC)>=0.58),]

write.table(dif_table_rapid_slow_filtered, file = "DEG_results_GSE44723_rapid_vs_slow_filtered_ensembl.csv", sep = "\t", col.names =TRUE)


write.table(dif_table_rapid_slow, file = "DEG_results_GSE44723_rapid_vs_slow_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)


#########Gene symbols##############################


dif_table_rapid_slow$ensg<-rownames(dif_table_rapid_slow)

table_annot_ipf <- merge(dif_table_rapid_slow, annotLookup, by.x="ensg", by.y="ensembl_gene_id")


table_annot_ipf_1<-table_annot_ipf[,c(2:8,10)]
aggr_table<- table_annot_ipf_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(dif_table_rapid_slow)[c(1:(length(colnames(dif_table_rapid_slow))-1))]

dif_table_rapid_slow_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_rapid_slow_filtered, file = "DEG_results_GSE44723_rapid_vs_slow_filtered_symbol.csv", sep = "\t", col.names =TRUE)

write.table(aggr_table, file = "DEG_results_GSE44723_rapid_vs_slow_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

###################DEG
#######dif_table_slow_healthy
#########################Ensembl_IDS##############################################
rownames_table<-c()
for (i in 1:(length(dif_table_slow_healthy$ID))) {
  rowname<-strsplit(as.character(dif_table_slow_healthy$ID[i]), "_")[[1]][1]
  rownames_table<-c(rownames_table, rowname)
}

dif_table_slow_healthy$ensg<-rownames_table

dif_table_slow_healthy<-dif_table_slow_healthy[-grep("AFFX", dif_table_slow_healthy$ensg),]

rownames(dif_table_slow_healthy)<-dif_table_slow_healthy$ensg

dif_table_slow_healthy<-dif_table_slow_healthy[,1:(length(colnames(dif_table_slow_healthy))-2)]

dif_table_slow_healthy_filtered<-dif_table_slow_healthy[which(dif_table_slow_healthy$adj.P.Val<=0.01 & abs(dif_table_slow_healthy$logFC)>=0.58),]

write.table(dif_table_slow_healthy_filtered, file = "DEG_results_GSE44723_slow_vs_healthy_filtered_ensembl.csv", sep = "\t", col.names =TRUE)


write.table(dif_table_slow_healthy, file = "DEG_results_GSE44723_slow_vs_healthy_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)


#########Gene symbols##############################


dif_table_slow_healthy$ensg<-rownames(dif_table_slow_healthy)

table_annot_ipf <- merge(dif_table_slow_healthy, annotLookup, by.x="ensg", by.y="ensembl_gene_id")


table_annot_ipf_1<-table_annot_ipf[,c(2:8,10)]
aggr_table<- table_annot_ipf_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(dif_table_slow_healthy)[c(1:(length(colnames(dif_table_slow_healthy))-1))]

dif_table_slow_healthy_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_slow_healthy_filtered, file = "DEG_results_GSE44723_slow_vs_healthy_filtered_symbol.csv", sep = "\t", col.names =TRUE)

write.table(aggr_table, file = "DEG_results_GSE44723_slow_vs_healthy_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

