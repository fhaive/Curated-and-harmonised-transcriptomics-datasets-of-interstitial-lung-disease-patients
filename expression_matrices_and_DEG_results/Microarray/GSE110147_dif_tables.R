###########Microarraydatasets_Ensembl to geneid################

setwd("/Lung_data/Biopsy/Microarray/GSE110147/expression_data/DEG_results")

library(xlsx)


dif_table_ipf<- as.data.frame(readxl::read_excel("RAW_Differential_Expression_Tables_2023-11-12.xlsx", sheet = 2))

dif_table_mixed_ipf_nsip <-as.data.frame(readxl::read_excel("RAW_Differential_Expression_Tables_2023-11-12.xlsx", sheet = 3))

dif_table_nsip<-as.data.frame(readxl::read_excel("RAW_Differential_Expression_Tables_2023-11-12.xlsx", sheet = 4))




setwd("/Lung_data/Biopsy/Microarray/GSE110147/expression_data/expression_matrices")

expr_mat<- read.table("RAW_Expression_Matrix_Normalized_2023-09-14.txt", sep="\t", header=T, row.names = NULL)

####################################


#rownames_table<-c()
#for (i in 1:(length(GSE110147_Mixed$ID))) {
  #rowname<-strsplit(as.character(GSE110147_Mixed$ID[i]), "_")[[1]][1]
  #rownames_table<-c(rownames_table, rowname)
#}

#rownames_table


#rownames(GSE110147_Mixed)<-rownames_table

#GSE110147_ipf<-GSE110147_ipf[,1:length(colnames(GSE110147_ipf))-1]

#GSE11196_ipf<-GSE110147_NsIp[- grep("AFFX", rownames(GSE110147_NsIp)),]

require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

annotLookup <- getBM(
  mart=mart,
  attributes=c(
    "affy_hugene_1_0_st_v1",
    "ensembl_gene_id",
    "gene_biotype",
    "external_gene_name"),
  filter = "affy_hugene_1_0_st_v1",
  values = expr_mat$row.names,
  uniqueRows = TRUE)

ensembl_ids = annotLookup$ensembl_gene_id
ensembl_ids = as.vector(ensembl_ids)

probes<-expr_mat$row.names

expr_mat$probes<-probes


table_with_gene_id <- merge(expr_mat, annotLookup, by.x="probes", by.y="affy_hugene_1_0_st_v1")

table_with_gene_id_1<-table_with_gene_id[,c(3:50, 53)] # REMOVE THE ENTREZ COLUMN
aggr_table<- table_with_gene_id_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(expr_mat)[c(2:49)]


write.table(aggr_table, file = "expression_matrix_symbol_GSE110147.csv", sep = "\t", col.names =TRUE)

table_with_gene_id_1<-table_with_gene_id[,c(3:51)] # REMOVE THE ENTREZ COLUMN
aggr_table<- table_with_gene_id_1 %>% group_by(ensembl_gene_id) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
#aggr_table<-aggr_table[-which(aggr_table$ensembl_gene_id == ""), ]
rownames(aggr_table)<-aggr_table$ensembl_gene_id
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(expr_mat)[c(2:49)]

write.table(aggr_table, file = "expression_matrix_ensembl_GSE110147.csv", sep = "\t", col.names =TRUE)


###########DIF_tables###########################################################################################

####GENE_ID

require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

annotLookup <- getBM(
  mart=mart,
  attributes=c(
    "affy_hugene_1_0_st_v1",
    "ensembl_gene_id",
    "gene_biotype",
    "external_gene_name"),
  filter = "affy_hugene_1_0_st_v1",
  values = dif_table_ipf$ID,
  uniqueRows = TRUE)

table_with_gene_id <- merge(dif_table_ipf, annotLookup, by.x="ID", by.y="affy_hugene_1_0_st_v1")

table_with_gene_id_1<-table_with_gene_id[,c(2:8, 11)] # REMOVE THE EXCESS COLUMN
aggr_table<- table_with_gene_id_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(dif_table_ipf[c(1:7)])

dif_table_ipf_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]


setwd("/Lung_data/Biopsy/Microarray/GSE110147/expression_data/DEG_results")
write.table(aggr_table, file = "DEG_results_GSE110147_ipf_vs_healthy_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)
write.table(dif_table_ipf_filtered, file = "DEG_results_GSE110147_ipf_vs_healthy_filtered_symbol.csv", sep = "\t", col.names =TRUE)

####ENSEMBL



table_with_gene_id_1<-table_with_gene_id[,c(2:9)] # REMOVE THE ENTREZ COLUMN
aggr_table<- table_with_gene_id_1 %>% group_by(ensembl_gene_id) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
#aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$ensembl_gene_id
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(dif_table_ipf[c(1:7)])

head(aggr_table)

dif_table_ipf_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]
write.table(aggr_table, file = "DEG_results_GSE110147_ipf_vs_healthy_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)
write.table(dif_table_ipf_filtered, file = "DEG_results_GSE110147_ipf_vs_healthy_filtered_ensembl.csv", sep = "\t", col.names =TRUE)

####################MIXED-IPF_NSIP########################################################

table_with_gene_id <- merge(dif_table_mixed_ipf_nsip, annotLookup, by.x="ID", by.y="affy_hugene_1_0_st_v1")

table_with_gene_id_1<-table_with_gene_id[,c(2:8, 11)] # REMOVE THE EXCESS COLUMN
aggr_table<- table_with_gene_id_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(dif_table_ipf[c(1:7)])

dif_table_mixed_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(aggr_table, file = "DEG_results_GSE110147_mixed_vs_healthy_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)
write.table(dif_table_mixed_filtered, file = "DEG_results_GSE110147_mixed_vs_healthy_filtered_symbol.csv", sep = "\t", col.names =TRUE)

####ENSEMBL



table_with_gene_id_1<-table_with_gene_id[,c(2:9)] # REMOVE THE ENTREZ COLUMN
aggr_table<- table_with_gene_id_1 %>% group_by(ensembl_gene_id) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
#aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$ensembl_gene_id
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(dif_table_ipf[c(1:7)])

head(aggr_table)

dif_table_mixed_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]
write.table(aggr_table, file = "DEG_results_GSE110147_mixed_vs_healthy_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)
write.table(dif_table_mixed_filtered, file = "DEG_results_GSE110147_mixed_vs_healthy_filtered_ensembl.csv", sep = "\t", col.names =TRUE)

####################NSIP########################################################

table_with_gene_id <- merge(dif_table_nsip, annotLookup, by.x="ID", by.y="affy_hugene_1_0_st_v1")

table_with_gene_id_1<-table_with_gene_id[,c(2:8, 11)] # REMOVE THE EXCESS COLUMN
aggr_table<- table_with_gene_id_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(dif_table_ipf[c(1:7)])

dif_table_nsip_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(aggr_table, file = "DEG_results_GSE110147_nsip_vs_healthy_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)
write.table(dif_table_nsip_filtered, file = "DEG_results_GSE110147_nsip_vs_healthy_filtered_symbol.csv", sep = "\t", col.names =TRUE)

####ENSEMBL



table_with_gene_id_1<-table_with_gene_id[,c(2:9)] # REMOVE THE ENTREZ COLUMN
aggr_table<- table_with_gene_id_1 %>% group_by(ensembl_gene_id) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
#aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$ensembl_gene_id
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(dif_table_ipf[c(1:7)])

head(aggr_table)

dif_table_NSIP_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]
write.table(aggr_table, file = "DEG_results_GSE110147_nsip_vs_healthy_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)
write.table(dif_table_NSIP_filtered, file = "DEG_results_GSE110147_nsip_vs_healthy_filtered_ensembl.csv", sep = "\t", col.names =TRUE)



