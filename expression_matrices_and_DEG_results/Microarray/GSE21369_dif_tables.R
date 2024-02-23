###########Microarraydatasets_Ensembl to geneid################

setwd("/Lung_data/Biopsy/Microarray/GSE21369/expression_data/DEG_results")

GSE21369_ipf<- as.data.frame(readxl::read_excel("RAW_Differential_Expression_Tables_2023-09-13.xlsx", sheet = 2))

GSE21369_cryp_org_pneu<- as.data.frame(readxl::read_excel("RAW_Differential_Expression_Tables_2023-09-13.xlsx", sheet = 3))

GSE21369_hypersens<- as.data.frame(readxl::read_excel("RAW_Differential_Expression_Tables_2023-09-13.xlsx", sheet = 4))

GSE21369_nsip<- as.data.frame(readxl::read_excel("RAW_Differential_Expression_Tables_2023-09-13.xlsx", sheet = 5))

GSE21369_rbild<- as.data.frame(readxl::read_excel("RAW_Differential_Expression_Tables_2023-09-13.xlsx", sheet = 6))

GSE21369_unchar_fib<- as.data.frame(readxl::read_excel("RAW_Differential_Expression_Tables_2023-09-13.xlsx", sheet = 7))

setwd("/Lung_data/Biopsy/Microarray/GSE21369/expression_data/expression_matrices")

expression_matrix<- read.table("RAW_Expression_Matrix_Normalized_2023-09-15.txt", sep="\t", header=T, row.names = NULL)



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
expression_matrix<-expression_matrix[,2:(length(colnames(expression_matrix))-1)]

expr_mat_ensemble<-expression_matrix

write.table(expr_mat_ensemble, file = "expression_matrix_ensembl_GSE21369.csv", sep = "\t", col.names =TRUE)

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
colnames(aggr_table)<-colnames(expression_matrix)
write.table(aggr_table, file = "expression_matrix_symbol_GSE21369.csv", sep = "\t", col.names =TRUE)

###############################DIF_TABLES################################
#################################IPF#####################################

###########
rm(list=ls())
setwd("/Lung_data/Biopsy/Microarray/GSE21369/expression_data/DEG_results")
GSE21369_ipf<- as.data.frame(readxl::read_excel("RAW_Differential_Expression_Tables_2023-09-13.xlsx", sheet = 2))

rownames_table<-c()
for (i in 1:(length(GSE21369_ipf$ID))) {
  rowname<-strsplit(as.character(GSE21369_ipf$ID[i]), "_")[[1]][1]
  rownames_table<-c(rownames_table, rowname)
}

GSE21369_ipf$ensg<-rownames_table

GSE21369_ipf<-GSE21369_ipf[-grep("AFFX", GSE21369_ipf$ensg),]


rownames(GSE21369_ipf)<-GSE21369_ipf$ensg

GSE21369_ipf<-GSE21369_ipf[,1:(length(colnames(GSE21369_ipf))-2)]

GSE21369_ipf<-GSE21369_ipf[order(rownames(GSE21369_ipf)),]

dif_table_ipf_filtered<-GSE21369_ipf[which(GSE21369_ipf$adj.P.Val<=0.01 & abs(GSE21369_ipf$logFC)>=0.58),]

write.table(dif_table_ipf_filtered, file = "DEG_results_GSE21369_ipf_vs_healthy_filtered_ensembl.csv", sep = "\t", col.names =TRUE)

write.table(GSE21369_ipf, file = "DEG_results_GSE21369_ipf_vs_healthy_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)


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
  values = rownames(GSE21369_ipf),
  uniqueRows = TRUE)

GSE21369_ipf$probes<-rownames(GSE21369_ipf)

table_with_gene_id <- merge(GSE21369_ipf, annotLookup, by.x="probes", by.y="ensembl_gene_id")

table_with_gene_id_1<-table_with_gene_id[,c(2,3,4,5,6,7,8,10)] # REMOVE THE ENTREZ COLUMN
aggr_table<- table_with_gene_id_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(GSE21369_ipf)[c(1:7)]

dif_table_ipf_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_ipf_filtered, file = "DEG_results_GSE21369_ipf_vs_healthy_filtered_symbol.csv", sep = "\t", col.names =TRUE)

write.table(aggr_table, file = "DEG_results_GSE21369_ipf_vs_healthy_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

###############################DIF_TABLES################################
#################################CRYP#####################################

###########
rm(list=ls())
GSE21369_cryp_org_pneu<- as.data.frame(readxl::read_excel("RAW_Differential_Expression_Tables_2023-09-13.xlsx", sheet = 3))

rownames_table<-c()
for (i in 1:(length(GSE21369_cryp_org_pneu$ID))) {
  rowname<-strsplit(as.character(GSE21369_cryp_org_pneu$ID[i]), "_")[[1]][1]
  rownames_table<-c(rownames_table, rowname)
}

GSE21369_cryp_org_pneu$ensg<-rownames_table

GSE21369_cryp_org_pneu<-GSE21369_cryp_org_pneu[-grep("AFFX", GSE21369_cryp_org_pneu$ensg),]


rownames(GSE21369_cryp_org_pneu)<-GSE21369_cryp_org_pneu$ensg

GSE21369_cryp_org_pneu<-GSE21369_cryp_org_pneu[,1:(length(colnames(GSE21369_cryp_org_pneu))-2)]

GSE21369_cryp_org_pneu<-GSE21369_cryp_org_pneu[order(rownames(GSE21369_cryp_org_pneu)),]

dif_table_cryp_filtered<-GSE21369_cryp_org_pneu[which(GSE21369_cryp_org_pneu$adj.P.Val<=0.01 & abs(GSE21369_cryp_org_pneu$logFC)>=0.58),]

write.table(dif_table_cryp_filtered, file = "DEG_results_GSE21369_cryp_org_pneu_vs_healthy_filtered_ensembl.csv", sep = "\t", col.names =TRUE)

write.table(GSE21369_cryp_org_pneu, file = "DEG_results_GSE21369_cryp_org_pneu_vs_healthy_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)


###########Symbol

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
  values = rownames(GSE21369_cryp_org_pneu),
  uniqueRows = TRUE)

GSE21369_cryp_org_pneu$probes<-rownames(GSE21369_cryp_org_pneu)

table_with_gene_id <- merge(GSE21369_cryp_org_pneu, annotLookup, by.x="probes", by.y="ensembl_gene_id")

table_with_gene_id_1<-table_with_gene_id[,c(2,3,4,5,6,7,8,10)] # REMOVE THE ENTREZ COLUMN
aggr_table<- table_with_gene_id_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(GSE21369_cryp_org_pneu)[c(1:7)]

dif_table_org_pneu_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_org_pneu_filtered, file = "DEG_results_GSE21369_cryp_org_pneu_vs_healthy_filtered_symbol.csv", sep = "\t", col.names =TRUE)

write.table(aggr_table, file = "DEG_results_GSE21369_cryp_org_pneu_vs_healthy_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

###############################DIF_TABLES################################
#################################hypersens#####################################

###########

rm(list=ls())
GSE21369_hypersens<- as.data.frame(readxl::read_excel("RAW_Differential_Expression_Tables_2023-09-13.xlsx", sheet = 4))

rownames_table<-c()
for (i in 1:(length(GSE21369_hypersens$ID))) {
  rowname<-strsplit(as.character(GSE21369_hypersens$ID[i]), "_")[[1]][1]
  rownames_table<-c(rownames_table, rowname)
}

GSE21369_hypersens$ensg<-rownames_table

GSE21369_hypersens<-GSE21369_hypersens[-grep("AFFX", GSE21369_hypersens$ensg),]


rownames(GSE21369_hypersens)<-GSE21369_hypersens$ensg

GSE21369_hypersens<-GSE21369_hypersens[,1:(length(colnames(GSE21369_hypersens))-2)]

GSE21369_hypersens<-GSE21369_hypersens[order(rownames(GSE21369_hypersens)),]

dif_table_hypersens_filtered<-GSE21369_hypersens[which(GSE21369_hypersens$adj.P.Val<=0.01 & abs(GSE21369_hypersens$logFC)>=0.58),]

write.table(dif_table_hypersens_filtered, file = "DEG_results_GSE21369_hypersens_vs_healthy_filtered_ensembl.csv", sep = "\t", col.names =TRUE)

write.table(GSE21369_hypersens, file = "DEG_results_GSE21369_hypersens_vs_healthy_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)

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
  values = rownames(GSE21369_hypersens),
  uniqueRows = TRUE)

GSE21369_hypersens$probes<-rownames(GSE21369_hypersens)

table_with_gene_id <- merge(GSE21369_hypersens, annotLookup, by.x="probes", by.y="ensembl_gene_id")

table_with_gene_id_1<-table_with_gene_id[,c(2,3,4,5,6,7,8,10)] # REMOVE THE ENTREZ COLUMN
aggr_table<- table_with_gene_id_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(GSE21369_hypersens)[c(1:7)]

dif_table_hypersens_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_hypersens_filtered, file = "DEG_results_GSE21369_hypersens_vs_healthy_filtered_symbol.csv", sep = "\t", col.names =TRUE)

write.table(aggr_table, file = "DEG_results_GSE21369_hypersens_vs_healthy_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

###############################DIF_TABLES################################
#################################nsip#####################################

###########
rm(list=ls())
GSE21369_nsip<- as.data.frame(readxl::read_excel("RAW_Differential_Expression_Tables_2023-09-13.xlsx", sheet = 5))

rownames_table<-c()
for (i in 1:(length(GSE21369_nsip$ID))) {
  rowname<-strsplit(as.character(GSE21369_nsip$ID[i]), "_")[[1]][1]
  rownames_table<-c(rownames_table, rowname)
}

GSE21369_nsip$ensg<-rownames_table

GSE21369_nsip<-GSE21369_nsip[-grep("AFFX", GSE21369_nsip$ensg),]


rownames(GSE21369_nsip)<-GSE21369_nsip$ensg

GSE21369_nsip<-GSE21369_nsip[,1:(length(colnames(GSE21369_nsip))-2)]

GSE21369_nsip<-GSE21369_nsip[order(rownames(GSE21369_nsip)),]

dif_table_nsip_filtered<-GSE21369_nsip[which(GSE21369_nsip$adj.P.Val<=0.01 & abs(GSE21369_nsip$logFC)>=0.58),]

write.table(dif_table_nsip_filtered, file = "DEG_results_GSE21369_nsip_vs_healthy_filtered_ensembl.csv", sep = "\t", col.names =TRUE)

write.table(GSE21369_nsip, file = "DEG_results_GSE21369_nsip_vs_healthy_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)


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
  values = rownames(GSE21369_nsip),
  uniqueRows = TRUE)

GSE21369_nsip$probes<-rownames(GSE21369_nsip)

table_with_gene_id <- merge(GSE21369_nsip, annotLookup, by.x="probes", by.y="ensembl_gene_id")

table_with_gene_id_1<-table_with_gene_id[,c(2,3,4,5,6,7,8,10)] # REMOVE THE ENTREZ COLUMN
aggr_table<- table_with_gene_id_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(GSE21369_nsip)[c(1:7)]

dif_table_nsip_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_nsip_filtered, file = "DEG_results_GSE21369_nsip_vs_healthy_filtered_symbol.csv", sep = "\t", col.names =TRUE)

write.table(aggr_table, file = "DEG_results_GSE21369_nsip_vs_healthy_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

###############################DIF_TABLES################################
#################################rbild#####################################

###########
rm(list=ls())
GSE21369_rbild<- as.data.frame(readxl::read_excel("RAW_Differential_Expression_Tables_2023-09-13.xlsx", sheet = 6))


rownames_table<-c()
for (i in 1:(length(GSE21369_rbild$ID))) {
  rowname<-strsplit(as.character(GSE21369_rbild$ID[i]), "_")[[1]][1]
  rownames_table<-c(rownames_table, rowname)
}

GSE21369_rbild$ensg<-rownames_table

GSE21369_rbild<-GSE21369_rbild[-grep("AFFX", GSE21369_rbild$ensg),]


rownames(GSE21369_rbild)<-GSE21369_rbild$ensg

GSE21369_rbild<-GSE21369_rbild[,1:(length(colnames(GSE21369_rbild))-2)]

GSE21369_rbild<-GSE21369_rbild[order(rownames(GSE21369_rbild)),]

dif_table_rbild_filtered<-GSE21369_rbild[which(GSE21369_rbild$adj.P.Val<=0.01 & abs(GSE21369_rbild$logFC)>=0.58),]

write.table(dif_table_rbild_filtered, file = "DEG_results_GSE21369_rbild_vs_healthy_filtered_ensembl.csv", sep = "\t", col.names =TRUE)

write.table(GSE21369_rbild, file = "DEG_results_GSE21369_rbild_vs_healthy_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)

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
  values = rownames(GSE21369_rbild),
  uniqueRows = TRUE)

GSE21369_rbild$probes<-rownames(GSE21369_rbild)

table_with_gene_id <- merge(GSE21369_rbild, annotLookup, by.x="probes", by.y="ensembl_gene_id")

table_with_gene_id_1<-table_with_gene_id[,c(2,3,4,5,6,7,8,10)] # REMOVE THE ENTREZ COLUMN
aggr_table<- table_with_gene_id_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(GSE21369_rbild)[c(1:7)]

dif_table_rbild_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_rbild_filtered, file = "DEG_results_GSE21369_rbild_vs_healthy_filtered_symbol.csv", sep = "\t", col.names =TRUE)

write.table(aggr_table, file = "DEG_results_GSE21369_rbild_vs_healthy_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

###############################DIF_TABLES################################
#################################unchar_fib#####################################

###########

rm(list=ls())
GSE21369_unchar_fib<- as.data.frame(readxl::read_excel("RAW_Differential_Expression_Tables_2023-09-13.xlsx", sheet = 7))



rownames_table<-c()
for (i in 1:(length(GSE21369_unchar_fib$ID))) {
  rowname<-strsplit(as.character(GSE21369_unchar_fib$ID[i]), "_")[[1]][1]
  rownames_table<-c(rownames_table, rowname)
}

GSE21369_unchar_fib$ensg<-rownames_table

GSE21369_unchar_fib<-GSE21369_unchar_fib[-grep("AFFX", GSE21369_unchar_fib$ensg),]


rownames(GSE21369_unchar_fib)<-GSE21369_unchar_fib$ensg

GSE21369_unchar_fib<-GSE21369_unchar_fib[,1:(length(colnames(GSE21369_unchar_fib))-2)]

GSE21369_unchar_fib<-GSE21369_unchar_fib[order(rownames(GSE21369_unchar_fib)),]

dif_table_unchar_fib_filtered<-GSE21369_unchar_fib[which(GSE21369_unchar_fib$adj.P.Val<=0.01 & abs(GSE21369_unchar_fib$logFC)>=0.58),]

write.table(dif_table_unchar_fib_filtered, file = "DEG_results_GSE21369_unchar_fib_vs_healthy_filtered_ensembl.csv", sep = "\t", col.names =TRUE)

write.table(GSE21369_unchar_fib, file = "DEG_results_GSE21369_unchar_fib_vs_healthy_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)

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
  values = rownames(GSE21369_unchar_fib),
  uniqueRows = TRUE)

GSE21369_unchar_fib$probes<-rownames(GSE21369_unchar_fib)

table_with_gene_id <- merge(GSE21369_unchar_fib, annotLookup, by.x="probes", by.y="ensembl_gene_id")

table_with_gene_id_1<-table_with_gene_id[,c(2,3,4,5,6,7,8,10)] # REMOVE THE ENTREZ COLUMN
aggr_table<- table_with_gene_id_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(GSE21369_unchar_fib)[c(1:7)]

dif_table_unchar_fib_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_unchar_fib_filtered, file = "DEG_results_GSE21369_unchar_fib_vs_healthy_filtered_symbol.csv", sep = "\t", col.names =TRUE)

write.table(aggr_table, file = "DEG_results_GSE21369_unchar_fib_vs_healthy_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)
