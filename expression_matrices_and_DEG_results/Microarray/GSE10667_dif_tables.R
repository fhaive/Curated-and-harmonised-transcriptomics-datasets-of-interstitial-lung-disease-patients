###########Microarraydatasets_Ensembl to geneid################

setwd("/Lung_data/Biopsy/Microarray/GSE10667/expression_data/DEG_results")

library(xlsx)


dif_table_ipf <- read.xlsx("RAW_Differential_Expression_Tables_2023-10-03.xlsx", sheetIndex = 2)

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
  values = as.character(dif_table_ipf$eListRaw.genes.ProbeName.ncIdx.),
  uniqueRows = TRUE)



table_with_gene_id <- merge(dif_table_ipf, annotLookup, by.x="eListRaw.genes.ProbeName.ncIdx.", by.y="agilent_wholegenome_4x44k_v1")

table_with_gene_id_1<-table_with_gene_id[,c(2,3,4,5,6,7,9,10)] # REMOVE excess columns
aggr_table<- table_with_gene_id_1 %>% group_by(ensembl_gene_id) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
#aggr_table<-aggr_table[-which(aggr_table$ensembl_gene_id == ""), ]
rownames(aggr_table)<-aggr_table$ensembl_gene_id
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(dif_table_ipf)[c(1:6,9)]
dif_table_ipf_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]
write.table(dif_table_ipf_filtered, file = "DEG_results_GSE10667_ipf_vs_healthy_filtered_ensembl.csv", sep = "\t", col.names =TRUE)
write.table(aggr_table, file = "DEG_results_GSE10667_ipf_vs_healthy_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)

###########GENE-IDS_FROM_Biomart###########################################################################################


table_with_gene_id_1<-table_with_gene_id[,c(2,3,4,5,6,7,9,12)]

aggr_table<- table_with_gene_id_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]

colnames(aggr_table)<-colnames(dif_table_ipf)[c(1:6,9)]
dif_table_ipf_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]
write.table(dif_table_ipf_filtered, file = "DEG_results_GSE10667_ipf_vs_healthy_filtered_symbol.csv", sep = "\t", col.names =TRUE)


write.table(aggr_table, file = "DEG_results_GSE10667_ipf_vs_healthy_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

####DIF_tables_Disease_state######################
#########Acute_healthy

setwd("/Lung_data/Biopsy/Microarray/GSE10667/expression_data/DEG_results")

acute_healthy <- as.data.frame(readxl::read_excel("RAW_Differential_Expression_Tables_2023-12-19_disease_state.xlsx", sheet = 2))

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
  values = as.character(acute_healthy$eListRaw.genes.ProbeName.ncIdx.),
  uniqueRows = TRUE)



table_with_gene_id <- merge(acute_healthy, annotLookup, by.x="eListRaw.genes.ProbeName.ncIdx.", by.y="agilent_wholegenome_4x44k_v1")

table_with_gene_id_1<-table_with_gene_id[,c(2,3,4,5,6,7,9,10)] # REMOVE excess columns
aggr_table<- table_with_gene_id_1 %>% group_by(ensembl_gene_id) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
#aggr_table<-aggr_table[-which(aggr_table$ensembl_gene_id == ""), ]
rownames(aggr_table)<-aggr_table$ensembl_gene_id
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(acute_healthy)[c(1:6,9)]
dif_table_acute_healthy_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]
write.table(dif_table_acute_healthy_filtered, file = "DEG_results_GSE10667_acute_vs_healthy_filtered_ensembl.csv", sep = "\t", col.names =TRUE)
write.table(aggr_table, file = "DEG_results_GSE10667_acute_vs_healthy_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)

###########GENE-IDS_FROM_Biomart###########################################################################################


table_with_gene_id_1<-table_with_gene_id[,c(2,3,4,5,6,7,9,12)]

aggr_table<- table_with_gene_id_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]

colnames(aggr_table)<-colnames(acute_healthy)[c(1:6,9)]
dif_table_acute_healthy_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]
write.table(dif_table_acute_healthy_filtered, file = "DEG_results_GSE10667_acute_vs_healthy_filtered_symbol.csv", sep = "\t", col.names =TRUE)


write.table(aggr_table, file = "DEG_results_GSE10667_acute_vs_healthy_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

####################################################################################
##########Usual_healthy

usual_healthy <- as.data.frame(readxl::read_excel("RAW_Differential_Expression_Tables_2023-12-19_disease_state.xlsx", sheet = 3))


##### probe names match with previous "annotLookup"
#require("biomaRt")
#mart <- useMart("ENSEMBL_MART_ENSEMBL")
#mart <- useDataset("hsapiens_gene_ensembl", mart)

#annotLookup <- getBM(
  #mart=mart,
  #attributes=c(
    #"agilent_wholegenome_4x44k_v1",
    #"ensembl_gene_id",
    #"gene_biotype",
    #"external_gene_name"),
  #filter = "agilent_wholegenome_4x44k_v1",
  #values = as.character(usual_healthy$eListRaw.genes.ProbeName.ncIdx.),
  #uniqueRows = TRUE)



table_with_gene_id <- merge(usual_healthy, annotLookup, by.x="eListRaw.genes.ProbeName.ncIdx.", by.y="agilent_wholegenome_4x44k_v1")

table_with_gene_id_1<-table_with_gene_id[,c(2,3,4,5,6,7,9,10)] # REMOVE excess columns
aggr_table<- table_with_gene_id_1 %>% group_by(ensembl_gene_id) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
#aggr_table<-aggr_table[-which(aggr_table$ensembl_gene_id == ""), ]
rownames(aggr_table)<-aggr_table$ensembl_gene_id
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(usual_healthy)[c(1:6,9)]
dif_table_usual_healthy_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]
write.table(dif_table_usual_healthy_filtered, file = "DEG_results_GSE10667_usual_vs_healthy_filtered_ensembl.csv", sep = "\t", col.names =TRUE)
write.table(aggr_table, file = "DEG_results_GSE10667_usual_vs_healthy_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)

###########GENE-IDS_FROM_Biomart###########################################################################################


table_with_gene_id_1<-table_with_gene_id[,c(2,3,4,5,6,7,9,12)]

aggr_table<- table_with_gene_id_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]

colnames(aggr_table)<-colnames(usual_healthy)[c(1:6,9)]
dif_table_usual_healthy_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]
write.table(dif_table_usual_healthy_filtered, file = "DEG_results_GSE10667_usual_vs_healthy_filtered_symbol.csv", sep = "\t", col.names =TRUE)


write.table(aggr_table, file = "DEG_results_GSE10667_usual_vs_healthy_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

####################################################################################
##########acute_usual

acute_usual <- as.data.frame(readxl::read_excel("RAW_Differential_Expression_Tables_2023-12-19_disease_state.xlsx", sheet = 4))

#require("biomaRt")
#mart <- useMart("ENSEMBL_MART_ENSEMBL")
#mart <- useDataset("hsapiens_gene_ensembl", mart)

#annotLookup <- getBM(
#mart=mart,
#attributes=c(
#"agilent_wholegenome_4x44k_v1",
#"ensembl_gene_id",
#"gene_biotype",
#"external_gene_name"),
#filter = "agilent_wholegenome_4x44k_v1",
#values = as.character(acute_usual$eListRaw.genes.ProbeName.ncIdx.),
#uniqueRows = TRUE)



table_with_gene_id <- merge(acute_usual, annotLookup, by.x="eListRaw.genes.ProbeName.ncIdx.", by.y="agilent_wholegenome_4x44k_v1")

table_with_gene_id_1<-table_with_gene_id[,c(2,3,4,5,6,7,9,10)] # REMOVE excess columns
aggr_table<- table_with_gene_id_1 %>% group_by(ensembl_gene_id) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
#aggr_table<-aggr_table[-which(aggr_table$ensembl_gene_id == ""), ]
rownames(aggr_table)<-aggr_table$ensembl_gene_id
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(acute_usual)[c(1:6,9)]
dif_table_acute_usual_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]
write.table(dif_table_acute_usual_filtered, file = "DEG_results_GSE10667_acute_vs_usual_filtered_ensembl.csv", sep = "\t", col.names =TRUE)
write.table(aggr_table, file = "DEG_results_GSE10667_acute_vs_usual_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)

###########GENE-IDS_FROM_Biomart###########################################################################################


table_with_gene_id_1<-table_with_gene_id[,c(2,3,4,5,6,7,9,12)]

aggr_table<- table_with_gene_id_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]

colnames(aggr_table)<-colnames(acute_usual)[c(1:6,9)]
dif_table_acute_usual_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]
write.table(dif_table_acute_usual_filtered, file = "DEG_results_GSE10667_acute_vs_usual_filtered_symbol.csv", sep = "\t", col.names =TRUE)


write.table(aggr_table, file = "DEG_results_GSE10667_acute_vs_usual_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)


##################expression_matrix##############################

setwd("/Lung_data/Biopsy/Microarray/GSE10667/expression_data/expression_matrices")

expression_mat_ipf <- read.table("RAW_Expression_Matrix_Aggregated_2023-10-03.txt", sep="\t", header=T)

####Note to self 
##Here is an dataframe where there are from differential expression table the "Probenames" and
#"Systemaic names" coming from eUTOPIA. The systematic names do not fully correspond (first 2744
#do correspond) to the probe names and that is why the agilent probenames doesnt map into 
#the rownames of the expression matrix which are the "systematic names".

check_frame<-data.frame(as.character(dif_table_ipf$eListRaw.genes.ProbeName.ncIdx.), as.character(dif_table_ipf$eListRaw.genes.SystematicName.ncIdx.),
                        as.character(rownames(expression_mat_ipf)))

colnames(check_frame)<-c("probename", "sys_name", "expr_mat_rowname")

sum(as.character(check_frame$probename)==as.character(check_frame$sys_name)) #1946

sum(as.character(check_frame$probename)==as.character(check_frame$expr_mat_rowname)) #1946

sum(as.character(check_frame$sys_name)==as.character(check_frame$expr_mat_rowname)) #28739

length(check_frame$expr_mat_rowname) #28739

##This means that the sys name and the expression matrix rownames are the same and both correspond 
#similar way to the probe names so I can change the rownames of the expression matrix to the 
#probenames from the differential expression table and take the gene annotations from there
#same results with identical()

rownames(expression_mat_ipf)<-dif_table_ipf$eListRaw.genes.ProbeName.ncIdx.

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
write.table(aggr_table, file = "expression_matrix_ensembl_GSE106675.csv", sep = "\t", col.names =TRUE)


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
write.table(aggr_table, file = "expression_matrix_symbol_GSE10667.csv", sep = "\t", col.names =TRUE)

