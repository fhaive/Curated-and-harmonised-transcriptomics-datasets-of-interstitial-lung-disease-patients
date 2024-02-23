###########Microarraydatasets_Ensembl to geneid################

setwd("/Lung_data/Biopsy/Microarray/GSE53845/expression_data/DEG_results")

library(xlsx)


dif_table_ipf <- read.xlsx("RAW_Differential_Expression_Tables_2023-09-19.xlsx", sheetIndex = 2)

require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

annotLookup <- getBM(
  mart=mart,
  attributes=c(
    "agilent_wholegenome_4x44k_v1",
    "ensembl_gene_id",
    "gene_biotype",
    "external_gene_name"),
  filter = "agilent_wholegenome_4x44k_v1",
  values = as.character(dif_table_ipf$ProbeName),
  uniqueRows = TRUE)



ensembl_ids = rownames(annotLookup$ensembl_gene_id)
ensembl_ids = as.vector(ensembl_ids)

probes<-dif_table_ipf$ProbeName

table_with_gene_id <- merge(dif_table_ipf, annotLookup, by.x="ProbeName", by.y="agilent_wholegenome_4x44k_v1")

table_with_gene_id_1<-table_with_gene_id[,c(2,3,4,5,6,7,9,10)] # REMOVE excess columns
aggr_table<- table_with_gene_id_1 %>% group_by(ensembl_gene_id) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
#aggr_table<-aggr_table[-which(aggr_table$ensembl_gene_id == ""), ]
rownames(aggr_table)<-aggr_table$ensembl_gene_id
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(dif_table_ipf)[c(1:6,9)]

write.table(aggr_table, file = "DEG_results_GSE53845_ipf_vs_healthy_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)

dif_table_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_filtered, file = "DEG_results_GSE53845_ipf_vs_healthy_filtered_ensembl.csv", sep = "\t", col.names =TRUE)


###########GENE-IDS_FROM_Biomart###########################################################################################


table_with_gene_id_1<-table_with_gene_id[,c(2,3,4,5,6,7,9,12)]

aggr_table<- table_with_gene_id_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]

colnames(aggr_table)<-colnames(dif_table_ipf)[c(1:6,9)]


write.table(aggr_table, file = "DEG_results_GSE53845_ipf_vs_healthy_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

dif_table_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_filtered, file = "DEG_results_GSE53845_ipf_vs_healthy_filtered_symbol.csv", sep = "\t", col.names =TRUE)


##################expression_matrix##############################

setwd("/Lung_data/Biopsy/Microarray/GSE53845/expression_data/expression_matrices")

expression_mat_ipf <- read.table("RAW_Expression_Matrix_Aggregated_2023-09-20.txt", sep="\t", header=T)

####Note to self 
##Here is an dataframe where there are from differential expression table the "Probenames" and
#"Systemaic names" coming from eUTOPIA. The systematic names do not fully correspond (first 2744
#do correspond) to the probe names and that is why the agilent probenames doesnt map into 
#the rownames of the expression matrix which are the "systematic names".

check_frame<-data.frame(as.character(dif_table_ipf$ProbeName), as.character(dif_table_ipf$SystematicName),
                        as.character(rownames(expression_mat_ipf)))

colnames(check_frame)<-c("probename", "sys_name", "expr_mat_rowname")

sum(as.character(check_frame$probename)==as.character(check_frame$sys_name)) #2744

sum(as.character(check_frame$probename)==as.character(check_frame$expr_mat_rowname)) #2744

sum(as.character(check_frame$sys_name)==as.character(check_frame$expr_mat_rowname)) #30846

length(check_frame$expr_mat_rowname) #30846

##This means that the sys name and the expression matrix rownames are the same and both correspond 
#similar way to the probe names so I can change the rownames of the expression matrix to the 
#probenames from the differential expression table and take the gene annotations from there

rownames(expression_mat_ipf)<-dif_table_ipf$ProbeName



require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

annotLookup <- getBM(
  mart=mart,
  attributes=c(
    "agilent_wholegenome_4x44k_v1",
    "ensembl_gene_id",
    "gene_biotype",
    "external_gene_name"),
  filter = "agilent_wholegenome_4x44k_v1",
  values = rownames(expression_mat_ipf),
  uniqueRows = T)




expression_mat_ipf$probes<-rownames(expression_mat_ipf)

table_with_gene_id <- merge(expression_mat_ipf, annotLookup, by.x=c("probes"), by.y=c("agilent_wholegenome_4x44k_v1"))

table_with_gene_id_1<-table_with_gene_id[,c(2:(length(colnames(table_with_gene_id))-2))] # Remove not relevant columns
aggr_table<- table_with_gene_id_1 %>% group_by(ensembl_gene_id) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
#aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$ensembl_gene_id
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(expression_mat_ipf)[c(2:length(colnames(expression_mat_ipf))-1)]
write.table(aggr_table, file = "expression_matrix_ensembl_GSE53845.csv", sep = "\t", col.names =TRUE)


############SYMBOL_EXPRMAT_BIOMART


table_with_gene_id <- merge(expression_mat_ipf, annotLookup, by.x=c("probes"), by.y=c("agilent_wholegenome_4x44k_v1"))

table_with_gene_id_1<-table_with_gene_id[,c(2:(length(colnames(table_with_gene_id))-3), length(colnames(table_with_gene_id)))] # Remove not relevant columns
aggr_table<- table_with_gene_id_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(expression_mat_ipf)[c(2:length(colnames(expression_mat_ipf))-1)]
write.table(aggr_table, file = "expression_matrix_symbol_GSE53845.csv", sep = "\t", col.names =TRUE)

