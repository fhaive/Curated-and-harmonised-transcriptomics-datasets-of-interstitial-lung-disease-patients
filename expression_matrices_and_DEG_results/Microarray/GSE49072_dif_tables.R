###########Microarraydatasets_Ensembl to geneid################

setwd("/Lung_data/Macrophage/Microarray/GSE49072/expression_data/DEG_results")


dif_table_ipf_healthy_vol <- as.data.frame(readxl::read_excel("ALL_Differential_Expression_Tables_2023-11-01.xlsx", sheet = 2))

dif_table_fpf_health_rela <- as.data.frame(readxl::read_excel("ALL_Differential_Expression_Tables_2023-11-01.xlsx", sheet = 3))

dif_table_ipf_vs_fpf <- as.data.frame(readxl::read_excel("ALL_Differential_Expression_Tables_2023-11-01.xlsx", sheet = 4))

dif_table_health_rela_health_vol <- as.data.frame(readxl::read_excel("ALL_Differential_Expression_Tables_2023-11-01.xlsx", sheet = 5))

dif_table_fpf_health_volu <- as.data.frame(readxl::read_excel("ALL_Differential_Expression_Tables_2023-11-01.xlsx", sheet = 6))

dif_table_ipf_health_rela <- as.data.frame(readxl::read_excel("ALL_Differential_Expression_Tables_2023-11-01.xlsx", sheet = 7))



#############dif_table_ipf_healthy_vol#################
#########################Ensembl_IDS##############################################
#########ENSG##############################

require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

annotLookup <- getBM(
  mart=mart,
  attributes=c(
    "affy_hg_u133a_2",
    "ensembl_gene_id",
    "gene_biotype",
    "external_gene_name"),
  filter = "affy_hg_u133a_2",
  values = dif_table_ipf_healthy_vol$ID,
  uniqueRows = TRUE)


table_annot<- merge(dif_table_ipf_healthy_vol, annotLookup, by.x="ID", by.y="affy_hg_u133a_2")

#ensembl
table_annot_1<-table_annot[,c(2:9)]
aggr_table<- table_annot_1 %>% group_by(ensembl_gene_id) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
#aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$ensembl_gene_id
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(dif_table_ipf_healthy_vol)[c(1:(length(colnames(dif_table_ipf_healthy_vol))-1))]

dif_table_ipf_healthy_vol_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_ipf_healthy_vol_filtered, file = "DEG_results_GSE49072_ipf_vs_healthy_volunteer_filtered_ensembl.csv", sep = "\t", col.names =TRUE)

write.table(aggr_table, file = "DEG_results_GSE49072_ipf_vs_healthy_volunteer_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)

table_annot_1<-table_annot[,c(2:8, 11)]
aggr_table<- table_annot_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(dif_table_ipf_healthy_vol)[c(1:(length(colnames(dif_table_ipf_healthy_vol))-1))]

dif_table_ipf_healthy_vol_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_ipf_healthy_vol_filtered, file = "DEG_results_GSE49072_ipf_vs_healthy_volunteer_filtered_symbol.csv", sep = "\t", col.names =TRUE)

write.table(aggr_table, file = "DEG_results_GSE49072_ipf_vs_healthy_volunteer_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)


#############dif_table_fpf_health_rela#################
#########################Ensembl_IDS##############################################
#########ENSG##############################



table_annot<- merge(dif_table_fpf_health_rela, annotLookup, by.x="ID", by.y="affy_hg_u133a_2")

#ensembl
table_annot_1<-table_annot[,c(2:9)]
aggr_table<- table_annot_1 %>% group_by(ensembl_gene_id) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
#aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$ensembl_gene_id
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(dif_table_fpf_health_rela)[c(1:(length(colnames(dif_table_fpf_health_rela))-1))]

dif_table_fpf_health_rela_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

DEG_results_GSE49072_ipf_vs_healthy_vol_unfiltered_symbol
write.table(dif_table_fpf_health_rela_filtered, file = "DEG_results_GSE49072_fpf_vs_healthy_relative_filtered_ensembl.csv", sep = "\t", col.names =TRUE)

write.table(aggr_table, file = "DEG_results_GSE49072_fpf_vs_healthy_relative_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)

table_annot_1<-table_annot[,c(2:8, 11)]
aggr_table<- table_annot_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(dif_table_fpf_health_rela)[c(1:(length(colnames(dif_table_fpf_health_rela))-1))]

dif_table_fpf_health_rela_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_fpf_health_rela_filtered, file = "DEG_results_GSE49072_fpf_vs_healthy_relative_filtered_symbol.csv", sep = "\t", col.names =TRUE)


write.table(aggr_table, file = "DEG_results_GSE49072_fpf_vs_healthy_relative_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

#############dif_table_ipf_vs_fpf#################
#########################Ensembl_IDS##############################################
#########ENSG##############################

identical(dif_table_fpf_health_rela$ID[order(dif_table_fpf_health_rela$ID)], dif_table_ipf_vs_fpf$ID[order(dif_table_ipf_vs_fpf$ID)])
###TRUE -> I can use the same annotLookup

#require("biomaRt")
#mart <- useMart("ENSEMBL_MART_ENSEMBL")
#mart <- useDataset("hsapiens_gene_ensembl", mart)

#annotLookup <- getBM(
  #mart=mart,
  #attributes=c(
    #"affy_hg_u133a_2",
    #"ensembl_gene_id",
    #"gene_biotype",
    #"external_gene_name"),
  #filter = "affy_hg_u133a_2",
  #values = dif_table_fpf_health_rela$ID,
  #uniqueRows = TRUE)


table_annot<- merge(dif_table_ipf_vs_fpf, annotLookup, by.x="ID", by.y="affy_hg_u133a_2")

#ensembl
table_annot_1<-table_annot[,c(2:9)]
aggr_table<- table_annot_1 %>% group_by(ensembl_gene_id) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
#aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$ensembl_gene_id
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(dif_table_ipf_vs_fpf)[c(1:(length(colnames(dif_table_ipf_vs_fpf))-1))]

dif_table_ipf_fpf_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_ipf_fpf_filtered, file = "DEG_results_GSE49072_ipf_vs_fpf_filtered_ensembl.csv", sep = "\t", col.names =TRUE)


write.table(aggr_table, file = "DEG_results_GSE49072_ipf_vs_fpf_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)

table_annot_1<-table_annot[,c(2:8, 11)]
aggr_table<- table_annot_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(dif_table_ipf_vs_fpf)[c(1:(length(colnames(dif_table_ipf_vs_fpf))-1))]

dif_table_ipf_fpf_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_ipf_fpf_filtered, file = "DEG_results_GSE49072_ipf_vs_fpf_filtered_symbol.csv", sep = "\t", col.names =TRUE)

write.table(aggr_table, file = "DEG_results_GSE49072_ipf_vs_fpf_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

#############dif_table_health_rela_health_vol#################
#########################Ensembl_IDS##############################################
#########ENSG##############################

identical(dif_table_fpf_health_rela$ID[order(dif_table_fpf_health_rela$ID)], dif_table_health_rela_health_vol$ID[order(dif_table_health_rela_health_vol$ID)])
###TRUE -> I can use the same annotLookup

#require("biomaRt")
#mart <- useMart("ENSEMBL_MART_ENSEMBL")
#mart <- useDataset("hsapiens_gene_ensembl", mart)

#annotLookup <- getBM(
#mart=mart,
#attributes=c(
#"affy_hg_u133a_2",
#"ensembl_gene_id",
#"gene_biotype",
#"external_gene_name"),
#filter = "affy_hg_u133a_2",
#values = dif_table_fpf_health_rela$ID,
#uniqueRows = TRUE)


table_annot<- merge(dif_table_health_rela_health_vol, annotLookup, by.x="ID", by.y="affy_hg_u133a_2")

#ensembl
table_annot_1<-table_annot[,c(2:9)]
aggr_table<- table_annot_1 %>% group_by(ensembl_gene_id) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
#aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$ensembl_gene_id
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(dif_table_health_rela_health_vol)[c(1:(length(colnames(dif_table_health_rela_health_vol))-1))]

dif_table_health_rela_health_vol_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_health_rela_health_vol_filtered, file = "DEG_results_GSE49072_healthy_relative_vs_healthy_volunteer_filtered_ensembl.csv", sep = "\t", col.names =TRUE)

write.table(aggr_table, file = "DEG_results_GSE49072_healthy_relative_vs_healthy_volunteer_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)

table_annot_1<-table_annot[,c(2:8, 11)]
aggr_table<- table_annot_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(dif_table_health_rela_health_vol)[c(1:(length(colnames(dif_table_health_rela_health_vol))-1))]

dif_table_health_rela_health_vol_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_health_rela_health_vol_filtered, file = "DEG_results_GSE49072_healthy_relative_vs_healthy_volunteer_filtered_symbol.csv", sep = "\t", col.names =TRUE)


write.table(aggr_table, file = "DEG_results_GSE49072_healthy_relative_vs_healthy_volunteer_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

#############dif_table_fpf_health_volu#################
#########################Ensembl_IDS##############################################
#########ENSG##############################

identical(dif_table_fpf_health_rela$ID[order(dif_table_fpf_health_rela$ID)], dif_table_fpf_health_volu$ID[order(dif_table_fpf_health_volu$ID)])
###TRUE -> I can use the same annotLookup

#require("biomaRt")
#mart <- useMart("ENSEMBL_MART_ENSEMBL")
#mart <- useDataset("hsapiens_gene_ensembl", mart)

#annotLookup <- getBM(
#mart=mart,
#attributes=c(
#"affy_hg_u133a_2",
#"ensembl_gene_id",
#"gene_biotype",
#"external_gene_name"),
#filter = "affy_hg_u133a_2",
#values = dif_table_fpf_health_rela$ID,
#uniqueRows = TRUE)


table_annot<- merge(dif_table_fpf_health_volu, annotLookup, by.x="ID", by.y="affy_hg_u133a_2")

#ensembl
table_annot_1<-table_annot[,c(2:9)]
aggr_table<- table_annot_1 %>% group_by(ensembl_gene_id) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
#aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$ensembl_gene_id
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(dif_table_fpf_health_volu)[c(1:(length(colnames(dif_table_fpf_health_volu))-1))]

dif_table_fpf_health_volu_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_fpf_health_volu_filtered, file = "DEG_results_GSE49072_fpf_vs_health_volunteer_filtered_ensembl.csv", sep = "\t", col.names =TRUE)


write.table(aggr_table, file = "DEG_results_GSE49072_fpf_vs_healthy_volunteer_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)

table_annot_1<-table_annot[,c(2:8, 11)]
aggr_table<- table_annot_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(dif_table_fpf_health_volu)[c(1:(length(colnames(dif_table_fpf_health_volu))-1))]

dif_table_fpf_health_volu_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_fpf_health_volu_filtered, file = "DEG_results_GSE49072_fpf_vs_healthy_volunteer_filtered_symbol.csv", sep = "\t", col.names =TRUE)


write.table(aggr_table, file = "DEG_results_GSE49072_fpf_vs_healthy_volunteer_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

#############dif_table_ipf_health_rela#################
#########################Ensembl_IDS##############################################
#########ENSG##############################

identical(dif_table_fpf_health_rela$ID[order(dif_table_fpf_health_rela$ID)], dif_table_ipf_health_rela$ID[order(dif_table_ipf_health_rela$ID)])
###TRUE -> I can use the same annotLookup

#require("biomaRt")
#mart <- useMart("ENSEMBL_MART_ENSEMBL")
#mart <- useDataset("hsapiens_gene_ensembl", mart)

#annotLookup <- getBM(
#mart=mart,
#attributes=c(
#"affy_hg_u133a_2",
#"ensembl_gene_id",
#"gene_biotype",
#"external_gene_name"),
#filter = "affy_hg_u133a_2",
#values = dif_table_fpf_health_rela$ID,
#uniqueRows = TRUE)


table_annot<- merge(dif_table_ipf_health_rela, annotLookup, by.x="ID", by.y="affy_hg_u133a_2")

#ensembl
table_annot_1<-table_annot[,c(2:9)]
aggr_table<- table_annot_1 %>% group_by(ensembl_gene_id) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
#aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$ensembl_gene_id
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(dif_table_ipf_health_rela)[c(1:(length(colnames(dif_table_ipf_health_rela))-1))]

dif_table_ipf_health_rela_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_ipf_health_rela_filtered, file = "DEG_results_GSE49072_ipf_vs_healthy_relative_filtered_ensembl.csv", sep = "\t", col.names =TRUE)

write.table(aggr_table, file = "DEG_results_GSE49072_ipf_vs_healthy_relative_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)

table_annot_1<-table_annot[,c(2:8, 11)]
aggr_table<- table_annot_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(dif_table_ipf_health_rela)[c(1:(length(colnames(dif_table_ipf_health_rela))-1))]


dif_table_ipf_health_rela_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_ipf_health_rela_filtered, file = "DEG_results_GSE49072_ipf_vs_healthy_relative_filtered_symbol.csv", sep = "\t", col.names =TRUE)

write.table(aggr_table, file = "DEG_results_GSE49072_ipf_vs_healthy_relative_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

#############expression_matrix#################
#########################Ensembl_IDS##############################################
#########ENSG##############################

setwd("/Lung_data/Macrophage/Microarray/GSE49072/expression_data/expression_matrices")

expr_mat<-read.table("RAW_Expression_Matrix_Normalized_2023-11-01.txt", sep="\t")

identical(dif_table_fpf_health_rela$ID[order(dif_table_fpf_health_rela$ID)], rownames(expr_mat)[order(rownames(expr_mat))])
###FALSE -> Lets make a new annotLookup just in case

expr_mat$probe<-rownames(expr_mat)
require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

annotLookup <- getBM(
mart=mart,
attributes=c(
"affy_hg_u133a_2",
"ensembl_gene_id",
"gene_biotype",
"external_gene_name"),
filter = "affy_hg_u133a_2",
values = expr_mat$probe,
uniqueRows = TRUE)


table_annot<- merge(expr_mat, annotLookup, by.x="probe", by.y="affy_hg_u133a_2")

#ensembl
table_annot_1<-table_annot[,c(2:86)]
aggr_table<- table_annot_1 %>% group_by(ensembl_gene_id) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
#aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$ensembl_gene_id
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(expr_mat)[c(1:(length(colnames(expr_mat))-1))]

write.table(aggr_table, file = "expression_matrix_ensembl_GSE49072.csv", sep = "\t", col.names =TRUE)

table_annot_1<-table_annot[,c(2:85, 88)]
aggr_table<- table_annot_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(expr_mat)[c(1:(length(colnames(expr_mat))-1))]

write.table(aggr_table, file = "expression_matrix_symbol_GSE49072.csv", sep = "\t", col.names =TRUE)

