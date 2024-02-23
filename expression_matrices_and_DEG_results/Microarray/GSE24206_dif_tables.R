###########Microarraydatasets_Ensembl to geneid################

setwd("/Lung_data/Biopsy/Microarray/GSE24206/expression_data/DEG_results")


GSE24206_ipf<- as.data.frame(readxl::read_excel("RAW_Differential_Expression_Tables_2023-02-02.xlsx", sheet = 2))

GSE24206_upper_vs_lower_ipf<-as.data.frame(readxl::read_excel("RAW_upper_lobe_vs_lower_lobe_Differential_Expression_Tables_2023-11-14.xlsx", sheet = 2))

GSE24206_upper_vs_healthy<-as.data.frame(readxl::read_excel("RAW_upper_lower_vs_healthy_Differential_Expression_Tables_2023-11-14.xlsx", sheet = 3))

GSE24206_lower_vs_healthy<-as.data.frame(readxl::read_excel("RAW_upper_lower_vs_healthy_Differential_Expression_Tables_2023-11-14.xlsx", sheet = 2))

GSE24206_advanced_vs_early_ipf<-as.data.frame(readxl::read_excel("RAW_disease_state_Differential_Expression_Tables_2023-11-14.xlsx", sheet = 2))

GSE24206_advanced_vs_healthy<-as.data.frame(readxl::read_excel("RAW_disease_state_Differential_Expression_Tables_2023-11-14.xlsx", sheet = 3))

GSE24206_early_vs_healthy<-as.data.frame(readxl::read_excel("RAW_disease_state_Differential_Expression_Tables_2023-11-14.xlsx", sheet = 4))

setwd("/Lung_data/Biopsy/Microarray/GSE24206/expression_data/expression_matrices")

expression_matrix <- read.table("RAW_Expression_Matrix_Normalized_2023-02-02.txt", sep="\t", header=T, row.names = NULL)

################Expression_Matrix##################3

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


write.table(expr_mat_ensemble, file = "expression_matrix_ensembl_GSE24206.csv", sep = "\t", col.names =TRUE)

####SYMBOL_EXPR_MAT_BIOMART

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
write.table(aggr_table, file = "expression_matrix_symbol_GSE24206.csv", sep = "\t", col.names =TRUE)

####################################
####################################DIF_exp
#IPF_all


rm(list=ls())
setwd("/Lung_data/Biopsy/Microarray/GSE24206/expression_data/DEG_results")

#IPF_vs_healthy

GSE24206_ipf<- as.data.frame(readxl::read_excel("RAW_Differential_Expression_Tables_2023-02-02.xlsx", sheet = 2))

rownames_table<-c()
for (i in 1:(length(GSE24206_ipf$ID))) {
  rowname<-strsplit(as.character(GSE24206_ipf$ID[i]), "_")[[1]][1]
  rownames_table<-c(rownames_table, rowname)
}

#rownames_table


GSE24206_ipf$ensg<-rownames_table

GSE24206_ipf<-GSE24206_ipf[-grep("AFFX", GSE24206_ipf$ensg),]

rownames(GSE24206_ipf)<-GSE24206_ipf$ensg

GSE24206_ipf<-GSE24206_ipf[,1:(length(colnames(GSE24206_ipf))-2)]

dif_table_ipf_filtered<-GSE24206_ipf[which(GSE24206_ipf$adj.P.Val<=0.01 & abs(GSE24206_ipf$logFC)>=0.58),]

write.table(dif_table_ipf_filtered, file = "DEG_results_GSE24206_ipf_all_vs_healthy_filtered_ensembl.csv", sep = "\t", col.names =TRUE)

write.table(GSE24206_ipf, file = "DEG_results_GSE24206_ipf_all_vs_healthy_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)

###############################Symbols_With_Biomart

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
  values = rownames(GSE24206_ipf),
  uniqueRows = TRUE)


GSE24206_ipf$probes<-rownames(GSE24206_ipf)

table_with_gene_id <- merge(GSE24206_ipf, annotLookup, by.x="probes", by.y="ensembl_gene_id")

table_with_gene_id_1<-table_with_gene_id[,c(2,3,4,5,6,7,8,10)] # REMOVE THE ENTREZ COLUMN
aggr_table<- table_with_gene_id_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(GSE24206_ipf)[c(1:7)]

dif_table_ipf_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_ipf_filtered, file = "DEG_results_GSE24206_ipf_all_vs_healthy_filtered_symbol.csv", sep = "\t", col.names =TRUE)

write.table(aggr_table, file = "DEG_results_GSE24206_ipf_all_vs_healthy_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

################################################################################################
#################################################################################################
####################################
####################################DIF_exp
#upper_vs_lower

rm(list=ls())

GSE24206_upper_vs_lower_ipf<-as.data.frame(readxl::read_excel("RAW_upper_lobe_vs_lower_lobe_Differential_Expression_Tables_2023-11-14.xlsx", sheet = 2))

rownames_table<-c()
for (i in 1:(length(GSE24206_upper_vs_lower_ipf$ID))) {
  rowname<-strsplit(as.character(GSE24206_upper_vs_lower_ipf$ID[i]), "_")[[1]][1]
  rownames_table<-c(rownames_table, rowname)
}

#rownames_table


GSE24206_upper_vs_lower_ipf$ensg<-rownames_table

GSE24206_upper_vs_lower_ipf<-GSE24206_upper_vs_lower_ipf[-grep("AFFX", GSE24206_upper_vs_lower_ipf$ensg),]

rownames(GSE24206_upper_vs_lower_ipf)<-GSE24206_upper_vs_lower_ipf$ensg

GSE24206_upper_vs_lower_ipf<-GSE24206_upper_vs_lower_ipf[,1:(length(colnames(GSE24206_upper_vs_lower_ipf))-2)]

dif_table_GSE24206_upper_vs_lower_ipf_filtered<-GSE24206_upper_vs_lower_ipf[which(GSE24206_upper_vs_lower_ipf$adj.P.Val<=0.01 & abs(GSE24206_upper_vs_lower_ipf$logFC)>=0.58),]

write.table(dif_table_GSE24206_upper_vs_lower_ipf_filtered, file = "DEG_results_GSE24206_upper_ipf_vs_lower_ipf_filtered_ensembl.csv", sep = "\t", col.names =TRUE)

write.table(GSE24206_upper_vs_lower_ipf, file = "DEG_results_GSE24206_upper_ipf_vs_lower_ipf_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)

###############################Symbols_With_Biomart

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
  values = rownames(GSE24206_upper_vs_lower_ipf),
  uniqueRows = TRUE)


GSE24206_upper_vs_lower_ipf$probes<-rownames(GSE24206_upper_vs_lower_ipf)

table_with_gene_id <- merge(GSE24206_upper_vs_lower_ipf, annotLookup, by.x="probes", by.y="ensembl_gene_id")

table_with_gene_id_1<-table_with_gene_id[,c(2,3,4,5,6,7,8,10)] 
aggr_table<- table_with_gene_id_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(GSE24206_upper_vs_lower_ipf)[c(1:7)]

dif_table_GSE24206_upper_vs_lower_ipf_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_GSE24206_upper_vs_lower_ipf_filtered, file = "DEG_results_GSE24206_upper_ipf_vs_lower_ipf_filtered_symbol.csv", sep = "\t", col.names =TRUE)

write.table(aggr_table, file = "DEG_results_GSE24206_upper_ipf_vs_lower_ipf_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)


################################################################################################
#################################################################################################
####################################
####################################DIF_exp
#upper_vs_healthy

rm(list=ls())

GSE24206_upper_vs_healthy<-as.data.frame(readxl::read_excel("RAW_upper_lower_vs_healthy_Differential_Expression_Tables_2023-11-14.xlsx", sheet = 3))

rownames_table<-c()
for (i in 1:(length(GSE24206_upper_vs_healthy$ID))) {
  rowname<-strsplit(as.character(GSE24206_upper_vs_healthy$ID[i]), "_")[[1]][1]
  rownames_table<-c(rownames_table, rowname)
}

#rownames_table


GSE24206_upper_vs_healthy$ensg<-rownames_table

GSE24206_upper_vs_healthy<-GSE24206_upper_vs_healthy[-grep("AFFX", GSE24206_upper_vs_healthy$ensg),]

rownames(GSE24206_upper_vs_healthy)<-GSE24206_upper_vs_healthy$ensg

GSE24206_upper_vs_healthy<-GSE24206_upper_vs_healthy[,1:(length(colnames(GSE24206_upper_vs_healthy))-2)]

dif_table_GSE24206_upper_vs_healthy_filtered<-GSE24206_upper_vs_healthy[which(GSE24206_upper_vs_healthy$adj.P.Val<=0.01 & abs(GSE24206_upper_vs_healthy$logFC)>=0.58),]

write.table(dif_table_GSE24206_upper_vs_healthy_filtered, file = "DEG_results_GSE24206_upper_ipf_vs_healthy_filtered_ensembl.csv", sep = "\t", col.names =TRUE)

write.table(GSE24206_upper_vs_healthy, file = "DEG_results_GSE24206_upper_ipf_vs_healthy_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)

###############################Symbols_With_Biomart


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
  values = rownames(GSE24206_upper_vs_healthy),
  uniqueRows = TRUE)

GSE24206_upper_vs_healthy$probes<-rownames(GSE24206_upper_vs_healthy)

table_with_gene_id <- merge(GSE24206_upper_vs_healthy, annotLookup, by.x="probes", by.y="ensembl_gene_id")

table_with_gene_id_1<-table_with_gene_id[,c(2,3,4,5,6,7,8,10)] # REMOVE THE ENTREZ COLUMN
aggr_table<- table_with_gene_id_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(GSE24206_upper_vs_healthy)[c(1:7)]

dif_table_GSE24206_upper_vs_healthy_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_GSE24206_upper_vs_healthy_filtered, file = "DEG_results_GSE24206_upper_ipf_vs_healthy_filtered_symbol.csv", sep = "\t", col.names =TRUE)

write.table(aggr_table, file = "DEG_results_GSE24206_upper_ipf_vs_healthy_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

#################################################################################################
####################################
####################################DIF_exp
#lower_vs_healthy


rm(list=ls())

GSE24206_lower_vs_healthy<-as.data.frame(readxl::read_excel("RAW_upper_lower_vs_healthy_Differential_Expression_Tables_2023-11-14.xlsx", sheet = 2))


rownames_table<-c()
for (i in 1:(length(GSE24206_lower_vs_healthy$ID))) {
  rowname<-strsplit(as.character(GSE24206_lower_vs_healthy$ID[i]), "_")[[1]][1]
  rownames_table<-c(rownames_table, rowname)
}

#rownames_table


GSE24206_lower_vs_healthy$ensg<-rownames_table

GSE24206_lower_vs_healthy<-GSE24206_lower_vs_healthy[-grep("AFFX", GSE24206_lower_vs_healthy$ensg),]

rownames(GSE24206_lower_vs_healthy)<-GSE24206_lower_vs_healthy$ensg

GSE24206_lower_vs_healthy<-GSE24206_lower_vs_healthy[,1:(length(colnames(GSE24206_lower_vs_healthy))-2)]

dif_table_GSE24206_lower_vs_healthy_filtered<-GSE24206_lower_vs_healthy[which(GSE24206_lower_vs_healthy$adj.P.Val<=0.01 & abs(GSE24206_lower_vs_healthy$logFC)>=0.58),]

write.table(dif_table_GSE24206_lower_vs_healthy_filtered, file = "DEG_results_GSE24206_lower_ipf_vs_healthy_filtered_ensembl.csv", sep = "\t", col.names =TRUE)

write.table(GSE24206_lower_vs_healthy, file = "DEG_results_GSE24206_lower_ipf_vs_healthy_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)

###############################Symbols_With_Biomart


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
  values = rownames(GSE24206_lower_vs_healthy),
  uniqueRows = TRUE)

GSE24206_lower_vs_healthy$probes<-rownames(GSE24206_lower_vs_healthy)

table_with_gene_id <- merge(GSE24206_lower_vs_healthy, annotLookup, by.x="probes", by.y="ensembl_gene_id")

table_with_gene_id_1<-table_with_gene_id[,c(2,3,4,5,6,7,8,10)] # REMOVE THE ENTREZ COLUMN
aggr_table<- table_with_gene_id_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(GSE24206_lower_vs_healthy)[c(1:7)]

dif_table_GSE24206_lower_vs_healthy_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_GSE24206_lower_vs_healthy_filtered, file = "DEG_results_GSE24206_lower_ipf_vs_healthy_filtered_symbol.csv", sep = "\t", col.names =TRUE)

write.table(aggr_table, file = "DEG_results_GSE24206_lower_ipf_vs_healthy_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

#################################################################################################
####################################
####################################DIF_exp
#GSE24206_advanced_vs_early_ipf
rm(list=ls())

GSE24206_advanced_vs_early_ipf<-as.data.frame(readxl::read_excel("RAW_disease_state_Differential_Expression_Tables_2023-11-14.xlsx", sheet = 2))

rownames_table<-c()
for (i in 1:(length(GSE24206_advanced_vs_early_ipf$ID))) {
  rowname<-strsplit(as.character(GSE24206_advanced_vs_early_ipf$ID[i]), "_")[[1]][1]
  rownames_table<-c(rownames_table, rowname)
}

#rownames_table


GSE24206_advanced_vs_early_ipf$ensg<-rownames_table

GSE24206_advanced_vs_early_ipf<-GSE24206_advanced_vs_early_ipf[-grep("AFFX", GSE24206_advanced_vs_early_ipf$ensg),]

rownames(GSE24206_advanced_vs_early_ipf)<-GSE24206_advanced_vs_early_ipf$ensg

GSE24206_advanced_vs_early_ipf<-GSE24206_advanced_vs_early_ipf[,1:(length(colnames(GSE24206_advanced_vs_early_ipf))-2)]

dif_table_GSE24206_advanced_vs_early_ipf_filtered<-GSE24206_advanced_vs_early_ipf[which(GSE24206_advanced_vs_early_ipf$adj.P.Val<=0.01 & abs(GSE24206_advanced_vs_early_ipf$logFC)>=0.58),]

write.table(dif_table_GSE24206_advanced_vs_early_ipf_filtered, file = "DEG_results_GSE24206_advanced_ipf_vs_early_ipf_filtered_ensembl.csv", sep = "\t", col.names =TRUE)

write.table(GSE24206_advanced_vs_early_ipf, file = "DEG_results_GSE24206_advanced_ipf_vs_early_ipf_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)

###############################Symbols_With_Biomart


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
  values = rownames(GSE24206_advanced_vs_early_ipf),
  uniqueRows = TRUE)

GSE24206_advanced_vs_early_ipf$probes<-rownames(GSE24206_advanced_vs_early_ipf)

table_with_gene_id <- merge(GSE24206_advanced_vs_early_ipf, annotLookup, by.x="probes", by.y="ensembl_gene_id")

table_with_gene_id_1<-table_with_gene_id[,c(2,3,4,5,6,7,8,10)] # REMOVE THE ENTREZ COLUMN
aggr_table<- table_with_gene_id_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(GSE24206_advanced_vs_early_ipf)[c(1:7)]

dif_table_GSE24206_advanced_vs_early_ipf_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_GSE24206_advanced_vs_early_ipf_filtered, file = "DEG_results_GSE24206_advanced_ipf_vs_early_ipf_filtered_symbol.csv", sep = "\t", col.names =TRUE)

write.table(aggr_table, file = "DEG_results_GSE24206_advanced_ipf_vs_early_ipf_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

#################################################################################################
####################################
####################################DIF_exp
#GSE24206_advanced_vs_healthy
rm(list=ls())

GSE24206_advanced_vs_healthy<-as.data.frame(readxl::read_excel("RAW_disease_state_Differential_Expression_Tables_2023-11-14.xlsx", sheet = 3))


rownames_table<-c()
for (i in 1:(length(GSE24206_advanced_vs_healthy$ID))) {
  rowname<-strsplit(as.character(GSE24206_advanced_vs_healthy$ID[i]), "_")[[1]][1]
  rownames_table<-c(rownames_table, rowname)
}

#rownames_table


GSE24206_advanced_vs_healthy$ensg<-rownames_table

GSE24206_advanced_vs_healthy<-GSE24206_advanced_vs_healthy[-grep("AFFX", GSE24206_advanced_vs_healthy$ensg),]

rownames(GSE24206_advanced_vs_healthy)<-GSE24206_advanced_vs_healthy$ensg

GSE24206_advanced_vs_healthy<-GSE24206_advanced_vs_healthy[,1:(length(colnames(GSE24206_advanced_vs_healthy))-2)]

dif_table_GSE24206_advanced_vs_healthy_filtered<-GSE24206_advanced_vs_healthy[which(GSE24206_advanced_vs_healthy$adj.P.Val<=0.01 & abs(GSE24206_advanced_vs_healthy$logFC)>=0.58),]

write.table(dif_table_GSE24206_advanced_vs_healthy_filtered, file = "DEG_results_GSE24206_advanced_ipf_vs_healthy_filtered_ensembl.csv", sep = "\t", col.names =TRUE)

write.table(GSE24206_advanced_vs_healthy, file = "DEG_results_GSE24206_advanced_ipf_vs_healthy_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)

###############################Symbols_With_Biomart

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
  values = rownames(GSE24206_advanced_vs_healthy),
  uniqueRows = TRUE)


GSE24206_advanced_vs_healthy$probes<-rownames(GSE24206_advanced_vs_healthy)

table_with_gene_id <- merge(GSE24206_advanced_vs_healthy, annotLookup, by.x="probes", by.y="ensembl_gene_id")

table_with_gene_id_1<-table_with_gene_id[,c(2,3,4,5,6,7,8,10)] # REMOVE THE ENTREZ COLUMN
aggr_table<- table_with_gene_id_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(GSE24206_advanced_vs_healthy)[c(1:7)]

dif_table_GSE24206_advanced_vs_healthy_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_GSE24206_advanced_vs_healthy_filtered, file = "DEG_results_GSE24206_advanced_ipf_vs_healthy_filtered_symbol.csv", sep = "\t", col.names =TRUE)

write.table(aggr_table, file = "DEG_results_GSE24206_advanced_ipf_vs_healthy_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

#################################################################################################
####################################
####################################DIF_exp
#GSE24206_early_vs_healthy
rm(list=ls())

GSE24206_early_vs_healthy<-as.data.frame(readxl::read_excel("RAW_disease_state_Differential_Expression_Tables_2023-11-14.xlsx", sheet = 4))


rownames_table<-c()
for (i in 1:(length(GSE24206_early_vs_healthy$ID))) {
  rowname<-strsplit(as.character(GSE24206_early_vs_healthy$ID[i]), "_")[[1]][1]
  rownames_table<-c(rownames_table, rowname)
}

#rownames_table


GSE24206_early_vs_healthy$ensg<-rownames_table

GSE24206_early_vs_healthy<-GSE24206_early_vs_healthy[-grep("AFFX", GSE24206_early_vs_healthy$ensg),]

rownames(GSE24206_early_vs_healthy)<-GSE24206_early_vs_healthy$ensg

GSE24206_early_vs_healthy<-GSE24206_early_vs_healthy[,1:(length(colnames(GSE24206_early_vs_healthy))-2)]

dif_table_GSE24206_early_vs_healthy_filtered<-GSE24206_early_vs_healthy[which(GSE24206_early_vs_healthy$adj.P.Val<=0.01 & abs(GSE24206_early_vs_healthy$logFC)>=0.58),]

write.table(dif_table_GSE24206_early_vs_healthy_filtered, file = "DEG_results_GSE24206_early_ipf_vs_healthy_filtered_ensembl.csv", sep = "\t", col.names =TRUE)

write.table(GSE24206_early_vs_healthy, file = "DEG_results_GSE24206_early_ipf_vs_healthy_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)

###############################Symbols_With_Biomart

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
  values = rownames(GSE24206_early_vs_healthy),
  uniqueRows = TRUE)


GSE24206_early_vs_healthy$probes<-rownames(GSE24206_early_vs_healthy)

table_with_gene_id <- merge(GSE24206_early_vs_healthy, annotLookup, by.x="probes", by.y="ensembl_gene_id")

table_with_gene_id_1<-table_with_gene_id[,c(2,3,4,5,6,7,8,10)] # REMOVE THE ENTREZ COLUMN
aggr_table<- table_with_gene_id_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(GSE24206_early_vs_healthy)[c(1:7)]

dif_table_GSE24206_early_vs_healthy_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_GSE24206_early_vs_healthy_filtered, file = "DEG_results_GSE24206_early_ipf_vs_healthy_filtered_symbol.csv", sep = "\t", col.names =TRUE)

write.table(aggr_table, file = "DEG_results_GSE24206_early_ipf_vs_healthy_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)



