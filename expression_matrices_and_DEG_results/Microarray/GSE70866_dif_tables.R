###########Microarraydatasets_Ensembl to geneid################

setwd("/Lung_data/BAL/Microarray/GSE70866/expression_data/DEG_results")

library(xlsx)


dif_table_ipf <- read.xlsx("RAW_Differential_Expression_Tables_2023-11-06.xlsx", sheetIndex = 2)


#########Dif_table_ensemble##################################3
require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

annotLookup <- getBM(
  mart=mart,
  attributes=c(
    "agilent_sureprint_g3_ge_8x60k",
    "ensembl_gene_id",
    "gene_biotype",
    "external_gene_name"),
  filter = "agilent_sureprint_g3_ge_8x60k",
  values = dif_table_ipf$eListRaw.genes.ProbeName.ncIdx.,
  uniqueRows = TRUE)


table_with_gene_id <- merge(dif_table_ipf, annotLookup, by.x="eListRaw.genes.ProbeName.ncIdx.", by.y="agilent_sureprint_g3_ge_8x60k")

table_with_gene_id_1<-table_with_gene_id[,c(2,3,4,5,6,7,9,10)] # REMOVE THE ENTREZ COLUMN
aggr_table<- table_with_gene_id_1 %>% group_by(ensembl_gene_id) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
#aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$ensembl_gene_id
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]

colnames(aggr_table)<-colnames(dif_table_ipf)[c(1:6,9)]

dif_table_ipf_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_ipf_filtered, file = "DEG_results_GSE70866_ipf_vs_healthy_filtered_ensembl.csv", sep = "\t", col.names =TRUE)


write.table(aggr_table, file = "DEG_results_GSE70866_ipf_vs_healthy_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)

########Symbol

table_with_gene_id_1<-table_with_gene_id[,c(2,3,4,5,6,7,9,12)] # REMOVE THE ENTREZ COLUMN
aggr_table<- table_with_gene_id_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(dif_table_ipf)[c(1:6,9)]

dif_table_ipf_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_ipf_filtered, file = "DEG_results_GSE70866_ipf_vs_healthy_filtered_symbol.csv", sep = "\t", col.names =TRUE)


write.table(aggr_table, file = "DEG_results_GSE70866_ipf_vs_healthy_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

###########expression_matrix###########################################################################################

setwd("/Lung_data/BAL/Microarray/GSE70866/expression_data/expression_matrices")

expr_mat <- read.table("RAW_Expression_Matrix_Aggregated_2023-11-06.txt", sep="\t")

test_frame<-data.frame(rownames(expr_mat), dif_table_ipf$eListRaw.genes.SystematicName.ncIdx., dif_table_ipf$eListRaw.genes.ProbeName.ncIdx.)
identical(test_frame$rownames.expr_mat., test_frame$dif_table_ipf.eListRaw.genes.SystematicName.ncIdx.)
#true

expr_mat$probes<-dif_table_ipf$eListRaw.genes.ProbeName.ncIdx.


table_with_gene_id <- merge(expr_mat, annotLookup, by.x=c("probes"), by.y=c("agilent_sureprint_g3_ge_8x60k"))

table_with_gene_id_1<-table_with_gene_id[,c(2:84)] # REMOVE THE ENTREZ COLUMN
aggr_table<- table_with_gene_id_1 %>% group_by(ensembl_gene_id) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
#aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$ensembl_gene_id
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(expr_mat)[c(1:82)]
write.table(aggr_table, file = "expression_matrix_ensembl_GSE70866.csv", sep = "\t", col.names =TRUE)

##########symbol
table_with_gene_id_1<-table_with_gene_id[,c(2:83, 86)] # REMOVE THE ENTREZ COLUMN
aggr_table<- table_with_gene_id_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(expr_mat)[c(1:82)]
write.table(aggr_table, file = "expression_matrix_symbol_GSE70866.csv", sep = "\t", col.names =TRUE)
