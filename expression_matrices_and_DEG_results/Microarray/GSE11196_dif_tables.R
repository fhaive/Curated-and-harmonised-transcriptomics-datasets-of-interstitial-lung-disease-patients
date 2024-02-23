###########Microarraydatasets_Ensembl to geneid################

setwd("/Lung_data/Fibroblast/Microarray/GSE11196/expression_data/DEG_results")

library(xlsx)


dif_table_ipf <- read.xlsx("ALL_Differential_Expression_Tables_2023-10-30.xlsx", sheetIndex = 2)

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

write.table(dif_table_ipf_filtered, file = "DEG_results_GSE11196_ipf_vs_healthy_filtered_ensembl.csv", sep = "\t", col.names =TRUE)


write.table(dif_table_ipf, file = "DEG_results_GSE11196_ipf_vs_healthy_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)


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

write.table(dif_table_ipf_filtered, file = "DEG_results_GSE11196_ipf_vs_healthy_filtered_symbol.csv", sep = "\t", col.names =TRUE)


write.table(aggr_table, file = "DEG_results_GSE11196_ipf_vs_healthy_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

###########Expression_matrix ensg#######################

setwd("/Lung_data/Fibroblast/Microarray/GSE11196/expression_data/expression_matrices")

raw_matrix<-read.table("RAW_Expression_Matrix_Normalized_2023-10-30.txt", sep="\t")

rownames_table<-c()
for (i in 1:(length(rownames(raw_matrix)))) {
  rowname<-strsplit(as.character(rownames(raw_matrix)[i]), "_")[[1]][1]
  rownames_table<-c(rownames_table, rowname)
}

raw_matrix$ensg<-rownames_table

raw_matrix<-raw_matrix[-grep("AFFX", raw_matrix$ensg),]

rownames(raw_matrix)<-raw_matrix$ensg

raw_matrix<-raw_matrix[,1:(length(colnames(raw_matrix))-1)]

write.table(raw_matrix, file = "expression_matrix_ensembl_GSE11196.csv", sep = "\t", col.names =TRUE)

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
  values = rownames(raw_matrix),
  uniqueRows = TRUE)

raw_matrix$ensg<-rownames(raw_matrix)

table_annot_ipf <- merge(raw_matrix, annotLookup, by.x="ensg", by.y="ensembl_gene_id")


table_annot_ipf_1<-table_annot_ipf[,c(2:25,27)]
aggr_table<- table_annot_ipf_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(raw_matrix)[c(1:(length(colnames(raw_matrix))-1))]

write.table(aggr_table, file = "expression_matrix_symbol_GSE11196.csv", sep = "\t", col.names =TRUE)
###############################################################################################################
