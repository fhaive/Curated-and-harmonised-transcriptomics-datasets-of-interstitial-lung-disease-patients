###########Microarraydatasets_Ensembl to geneid################

setwd("/Lung_data/Biopsy/Microarray/GSE72073/expression_data/DEG_results")


GSE72073_ipf<- as.data.frame(readxl::read_excel("RAW_Differential_Expression_Tables_2023-02-03.xlsx", sheet = 2))

setwd("/Lung_data/Biopsy/Microarray/GSE72073/expression_data/expression_matrices")

expression_matrix<- read.table("RAW_Expression_Matrix_Normalized_2023-02-03.txt", sep="\t", header=T, row.names = NULL)


###############DEG_results#####################

setwd("/Lung_data/Biopsy/Microarray/GSE72073/expression_data/DEG_results")

rownames_table<-c()
for (i in 1:(length(GSE72073_ipf$ID))) {
  rowname<-strsplit(as.character(GSE72073_ipf$ID[i]), "_")[[1]][1]
  rownames_table<-c(rownames_table, rowname)
}

#rownames_table

GSE72073_ipf$ensg<-rownames_table

GSE72073_ipf<-GSE72073_ipf[-grep("AFFX", GSE72073_ipf$ensg),]

rownames(GSE72073_ipf)<-GSE72073_ipf$ensg

GSE72073_ipf<-GSE72073_ipf[,1:(length(colnames(GSE72073_ipf))-2)]

dif_table_ipf_filtered<-GSE72073_ipf[which(GSE72073_ipf$adj.P.Val<=0.01 & abs(GSE72073_ipf$logFC)>=0.58),]

write.table(dif_table_ipf_filtered, file = "DEG_results_GSE72073_ipf_vs_healthy_filtered_ensembl.csv", sep = "\t", col.names =TRUE)

write.table(GSE72073_ipf, file = "DEG_results_GSE72073_ipf_vs_healthy_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)



###########GENE-IDS###########################################################################################
###############BIOMART_GENE_SYMBOL


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
  values = rownames(GSE72073_ipf),
  uniqueRows = TRUE)

ensembl_ids = annotLookup$ensembl_gene_id
ensembl_ids = as.vector(ensembl_ids)

GSE72073_ipf$ensg<-rownames(GSE72073_ipf)

table_with_gene_id <- merge(GSE72073_ipf, annotLookup, by.x="ensg", by.y="ensembl_gene_id")

table_with_gene_id_1<-table_with_gene_id[,c(2,3,4,5,6,7,8,10)] # REMOVE THE ENTREZ COLUMN
aggr_table<- table_with_gene_id_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(GSE72073_ipf)[c(1:7)]

dif_table_ipf_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]
write.table(dif_table_ipf_filtered, file = "DEG_results_GSE72073_ipf_vs_healthy_filtered_symbol.csv", sep = "\t", col.names =TRUE)

write.table(aggr_table, file = "DEG_results_GSE72073_ipf_vs_healty_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

################Expression_Matrix##################

setwd("/Lung_data/Biopsy/Microarray/GSE72073/expression_data/expression_matrices")

rownames_table<-c()
for (i in 1:(length(expression_matrix$row.names))) {
  rowname<-strsplit(as.character(expression_matrix$row.names[i]), "_")[[1]][1]
  rownames_table<-c(rownames_table, rowname)
}

expression_matrix$ensg<-rownames_table

expression_matrix<-expression_matrix[-grep("AFFX", expression_matrix$ensg),]


#aggr_table<- expression_matrix %>% group_by(ensg) %>% dplyr::summarise_all(.funs = c(median = "median"))
#aggr_table <- as.data.frame(aggr_table)
#aggr_table <- na.omit(aggr_table)

rownames(expression_matrix)<-expression_matrix$ensg
#aggr_table<-aggr_table[,3:length(colnames(aggr_table))]

expr_mat_ensemble<-expression_matrix[,2:(length(colnames(expression_matrix))-1)]

write.table(expr_mat_ensemble, file = "expression_matrix_ensembl_GSE72073.csv", sep = "\t", col.names =TRUE)

#############Symbol_With_biomart##############3

require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

annotLookup <- getBM(mart=mart,
                     attributes=c(
                       "ensembl_gene_id",
                       "gene_biotype",
                       "external_gene_name"),
                     filter = "ensembl_gene_id",
                     values = rownames(expr_mat_ensemble),
                     uniqueRows = TRUE)

symbol = annotLookup$external_gene_name
symbol = as.vector(symbol)

expr_mat_ensemble$probes<-rownames(expr_mat_ensemble)


table_with_gene_id <- merge(expr_mat_ensemble, annotLookup, by.x=c("probes"), by.y=c("ensembl_gene_id"))

table_with_gene_id_1<-table_with_gene_id[,c(2:(length(colnames(table_with_gene_id))-2),length(colnames(table_with_gene_id)))] # REMOVE THE ENTREZ COLUMN
aggr_table<- table_with_gene_id_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(expression_matrix)[c(2:(length(colnames(expression_matrix))-1))]
write.table(aggr_table, file = "expression_matrix_symbol_GSE72073.csv", sep = "\t", col.names =TRUE)



