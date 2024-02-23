###########Microarraydatasets_Ensembl to geneid################

setwd("/Lung_data/Macrophage/Microarray/GSE90010/expression_data/DEG_results")

library(xlsx)


am_ipf_MDM_control <- as.data.frame(readxl::read_excel("ALL_Differential_Expression_Tables_2023-11-01.xlsx", sheet  = 2))

am_rb_ild_MDM_control <-  as.data.frame(readxl::read_excel("ALL_Differential_Expression_Tables_2023-11-01.xlsx", sheet  = 3))

MDM_ap_neut_vs_control <- as.data.frame(readxl::read_excel("ALL_Differential_Expression_Tables_2023-11-01.xlsx", sheet = 4))

MDM_ap_neut_with_lps_vs_control<- as.data.frame(readxl::read_excel("ALL_Differential_Expression_Tables_2023-11-01.xlsx", sheet = 5))

MDM_with_lps_vs_control <-as.data.frame(readxl::read_excel("ALL_Differential_Expression_Tables_2023-11-01.xlsx", sheet = 6))

MDM_ap_neut_with_lps_vs_MDM_ap_neutr<- as.data.frame(readxl::read_excel("ALL_Differential_Expression_Tables_2023-11-01.xlsx", sheet= 7))


############################################################################################################
############am_ipf_mdm_control####################################################
rownames_table<-c()
for (i in 1:(length(am_ipf_MDM_control$ID))) {
  rowname<-strsplit(as.character(am_ipf_MDM_control$ID[i]), "_")[[1]][1]
  rownames_table<-c(rownames_table, rowname)
}

#rownames_table


am_ipf_MDM_control$ensg<-rownames_table

#am_ipf_MDM_control<-am_ipf_MDM_control[-grep("AFFX", am_ipf_MDM_control$ensg),]

rownames(am_ipf_MDM_control)<-am_ipf_MDM_control$ensg

am_ipf_MDM_control<-am_ipf_MDM_control[,1:(length(colnames(am_ipf_MDM_control))-2)]

am_ipf_MDM_control<-am_ipf_MDM_control[order(rownames(am_ipf_MDM_control), decreasing = F),]

dif_table_am_ipf_MDM_control_filtered<-am_ipf_MDM_control[which(am_ipf_MDM_control$adj.P.Val<=0.01 & abs(am_ipf_MDM_control$logFC)>=0.58),]

write.table(dif_table_am_ipf_MDM_control_filtered, file = "DEG_results_GSE90010_am_ipf_vs_mdm_control_filtered_ensembl.csv", sep = "\t", col.names =TRUE)

write.table(am_ipf_MDM_control, file = "DEG_results_GSE90010_am_ipf_vs_mdm_control_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)

########symbol

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
  values = rownames(am_ipf_MDM_control),
  uniqueRows = TRUE)

am_ipf_MDM_control$ensg<-rownames(am_ipf_MDM_control)

table_annot<- merge(am_ipf_MDM_control, annotLookup, by.x="ensg", by.y="ensembl_gene_id")

table_annot_1<-table_annot[,c(2:8, 10)]
aggr_table<- table_annot_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(am_ipf_MDM_control)[c(1:(length(colnames(am_ipf_MDM_control))-1))]

dif_table_am_ipf_MDM_control_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_am_ipf_MDM_control_filtered, file = "DEG_results_GSE90010_am_ipf_vs_mdm_control_filtered_symbol.csv", sep = "\t", col.names =TRUE)


write.table(aggr_table, file = "DEG_results_GSE90010_am_ipf_vs_mdm_control_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

############################################################################################################
############am_rb_ild_MDM_control####################################################
rownames_table<-c()
for (i in 1:(length(am_rb_ild_MDM_control$ID))) {
  rowname<-strsplit(as.character(am_rb_ild_MDM_control$ID[i]), "_")[[1]][1]
  rownames_table<-c(rownames_table, rowname)
}

#rownames_table


am_rb_ild_MDM_control$ensg<-rownames_table

#am_rb_ild_MDM_control<-am_rb_ild_MDM_control[-grep("AFFX", am_rb_ild_MDM_control$ensg),]

rownames(am_rb_ild_MDM_control)<-am_rb_ild_MDM_control$ensg

am_rb_ild_MDM_control<-am_rb_ild_MDM_control[,1:(length(colnames(am_rb_ild_MDM_control))-2)]

am_rb_ild_MDM_control<-am_rb_ild_MDM_control[order(rownames(am_rb_ild_MDM_control), decreasing = F),]

dif_table_am_rb_ild_MDM_control_filtered<-am_rb_ild_MDM_control[which(am_rb_ild_MDM_control$adj.P.Val<=0.01 & abs(am_rb_ild_MDM_control$logFC)>=0.58),]

write.table(dif_table_am_rb_ild_MDM_control_filtered, file = "DEG_results_GSE90010_am_rb_ild_vs_mdm_control_filtered_ensembl.csv", sep = "\t", col.names =TRUE)

write.table(am_rb_ild_MDM_control, file = "DEG_results_GSE90010_am_rb_ild_vs_mdm_control_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)

########symbol


am_rb_ild_MDM_control$ensg<-rownames(am_rb_ild_MDM_control)

table_annot<- merge(am_rb_ild_MDM_control, annotLookup, by.x="ensg", by.y="ensembl_gene_id")

table_annot_1<-table_annot[,c(2:8, 10)]
aggr_table<- table_annot_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(am_rb_ild_MDM_control)[c(1:(length(colnames(am_rb_ild_MDM_control))-1))]

dif_table_am_rb_ild_MDM_control_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_am_rb_ild_MDM_control_filtered, file = "DEG_results_GSE90010_am_rb_ild_vs_mdm_control_filtered_symbol.csv", sep = "\t", col.names =TRUE)


write.table(aggr_table, file = "DEG_results_GSE90010_am_rb_ild_vs_mdm_control_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

############################################################################################################
############MDM_ap_neut_vs_control####################################################
rownames_table<-c()
for (i in 1:(length(MDM_ap_neut_vs_control$ID))) {
  rowname<-strsplit(as.character(MDM_ap_neut_vs_control$ID[i]), "_")[[1]][1]
  rownames_table<-c(rownames_table, rowname)
}

#rownames_table


MDM_ap_neut_vs_control$ensg<-rownames_table

#MDM_ap_neut_vs_control<-MDM_ap_neut_vs_control[-grep("AFFX", MDM_ap_neut_vs_control$ensg),]

rownames(MDM_ap_neut_vs_control)<-MDM_ap_neut_vs_control$ensg

MDM_ap_neut_vs_control<-MDM_ap_neut_vs_control[,1:(length(colnames(MDM_ap_neut_vs_control))-2)]

MDM_ap_neut_vs_control<-MDM_ap_neut_vs_control[order(rownames(MDM_ap_neut_vs_control), decreasing = F),]

dif_table_MDM_ap_neut_vs_control_filtered<-MDM_ap_neut_vs_control[which(MDM_ap_neut_vs_control$adj.P.Val<=0.01 & abs(MDM_ap_neut_vs_control$logFC)>=0.58),]

write.table(dif_table_MDM_ap_neut_vs_control_filtered, file = "DEG_results_GSE90010_mdm_ap_neut_vs_mdm_control_filtered_ensembl.csv", sep = "\t", col.names =TRUE)

write.table(MDM_ap_neut_vs_control, file = "DEG_results_GSE90010_mdm_ap_neut_vs_mdm_control_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)

########symbol

MDM_ap_neut_vs_control$ensg<-rownames(MDM_ap_neut_vs_control)

table_annot<- merge(MDM_ap_neut_vs_control, annotLookup, by.x="ensg", by.y="ensembl_gene_id")

table_annot_1<-table_annot[,c(2:8, 10)]
aggr_table<- table_annot_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(MDM_ap_neut_vs_control)[c(1:(length(colnames(MDM_ap_neut_vs_control))-1))]

dif_table_MDM_ap_neut_vs_control_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_MDM_ap_neut_vs_control_filtered, file = "DEG_results_GSE90010_mdm_ap_neut_vs_mdm_control_filtered_symbol.csv", sep = "\t", col.names =TRUE)


write.table(aggr_table, file = "DEG_results_GSE90010_mdm_ap_neut_vs_mdm_control_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

############################################################################################################
############MDM_ap_neut_with_lps_vs_control####################################################
rownames_table<-c()
for (i in 1:(length(MDM_ap_neut_with_lps_vs_control$ID))) {
  rowname<-strsplit(as.character(MDM_ap_neut_with_lps_vs_control$ID[i]), "_")[[1]][1]
  rownames_table<-c(rownames_table, rowname)
}

#rownames_table


MDM_ap_neut_with_lps_vs_control$ensg<-rownames_table

#MDM_ap_neut_with_lps_vs_control<-MDM_ap_neut_with_lps_vs_control[-grep("AFFX", MDM_ap_neut_with_lps_vs_control$ensg),]

rownames(MDM_ap_neut_with_lps_vs_control)<-MDM_ap_neut_with_lps_vs_control$ensg

MDM_ap_neut_with_lps_vs_control<-MDM_ap_neut_with_lps_vs_control[,1:(length(colnames(MDM_ap_neut_with_lps_vs_control))-2)]

MDM_ap_neut_with_lps_vs_control<-MDM_ap_neut_with_lps_vs_control[order(rownames(MDM_ap_neut_with_lps_vs_control), decreasing = F),]

dif_table_MDM_ap_neut_with_lps_vs_control_filtered<-MDM_ap_neut_with_lps_vs_control[which(MDM_ap_neut_with_lps_vs_control$adj.P.Val<=0.01 & abs(MDM_ap_neut_with_lps_vs_control$logFC)>=0.58),]

write.table(dif_table_MDM_ap_neut_with_lps_vs_control_filtered, file = "DEG_results_GSE90010_mdm_ap_neut_with_lps_vs_mdm_control_filtered_ensembl.csv", sep = "\t", col.names =TRUE)


write.table(MDM_ap_neut_with_lps_vs_control, file = "DEG_results_GSE90010_mdm_ap_neut_with_lps_vs_mdm_control_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)

########symbol

MDM_ap_neut_with_lps_vs_control$ensg<-rownames(MDM_ap_neut_with_lps_vs_control)

table_annot<- merge(MDM_ap_neut_with_lps_vs_control, annotLookup, by.x="ensg", by.y="ensembl_gene_id")

table_annot_1<-table_annot[,c(2:8, 10)]
aggr_table<- table_annot_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(MDM_ap_neut_with_lps_vs_control)[c(1:(length(colnames(MDM_ap_neut_with_lps_vs_control))-1))]

dif_table_MDM_ap_neut_with_lps_vs_control_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_MDM_ap_neut_with_lps_vs_control_filtered, file = "DEG_results_GSE90010_mdm_ap_neut_with_lps_vs_mdm_control_filtered_symbol.csv", sep = "\t", col.names =TRUE)


write.table(aggr_table, file = "DEG_results_GSE90010_mdm_ap_neut_with_lps_vs_mdm_control_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

############################################################################################################
############MDM_with_lps_vs_control####################################################
rownames_table<-c()
for (i in 1:(length(MDM_with_lps_vs_control$ID))) {
  rowname<-strsplit(as.character(MDM_with_lps_vs_control$ID[i]), "_")[[1]][1]
  rownames_table<-c(rownames_table, rowname)
}

#rownames_table


MDM_with_lps_vs_control$ensg<-rownames_table

#MDM_with_lps_vs_control<-MDM_with_lps_vs_control[-grep("AFFX", MDM_with_lps_vs_control$ensg),]

rownames(MDM_with_lps_vs_control)<-MDM_with_lps_vs_control$ensg

MDM_with_lps_vs_control<-MDM_with_lps_vs_control[,1:(length(colnames(MDM_with_lps_vs_control))-2)]

MDM_with_lps_vs_control<-MDM_with_lps_vs_control[order(rownames(MDM_with_lps_vs_control), decreasing = F),]

dif_table_MDM_with_lps_vs_control_filtered<-MDM_with_lps_vs_control[which(MDM_with_lps_vs_control$adj.P.Val<=0.01 & abs(MDM_with_lps_vs_control$logFC)>=0.58),]

write.table(dif_table_MDM_with_lps_vs_control_filtered, file = "DEG_results_GSE90010_mdm_with_lps_vs_mdm_control_filtered_ensembl.csv", sep = "\t", col.names =TRUE)

write.table(MDM_with_lps_vs_control, file = "DEG_results_GSE90010_mdm_with_lps_vs_mdm_control_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)

########symbol


MDM_with_lps_vs_control$ensg<-rownames(MDM_with_lps_vs_control)

table_annot<- merge(MDM_with_lps_vs_control, annotLookup, by.x="ensg", by.y="ensembl_gene_id")

table_annot_1<-table_annot[,c(2:8, 10)]
aggr_table<- table_annot_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(MDM_with_lps_vs_control)[c(1:(length(colnames(MDM_with_lps_vs_control))-1))]

dif_table_MDM_with_lps_vs_control_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_MDM_with_lps_vs_control_filtered, file = "DEG_results_GSE90010_mdm_with_lps_vs_mdm_control_filtered_symbol.csv", sep = "\t", col.names =TRUE)


write.table(aggr_table, file = "DEG_results_GSE90010_mdm_with_lps_vs_mdm_control_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

############################################################################################################
############MDM_ap_neut_with_lps_vs_MDM_ap_neutr####################################################
rownames_table<-c()
for (i in 1:(length(MDM_ap_neut_with_lps_vs_MDM_ap_neutr$ID))) {
  rowname<-strsplit(as.character(MDM_ap_neut_with_lps_vs_MDM_ap_neutr$ID[i]), "_")[[1]][1]
  rownames_table<-c(rownames_table, rowname)
}

#rownames_table


MDM_ap_neut_with_lps_vs_MDM_ap_neutr$ensg<-rownames_table

#MDM_ap_neut_with_lps_vs_MDM_ap_neutr<-MDM_ap_neut_with_lps_vs_MDM_ap_neutr[-grep("AFFX", MDM_ap_neut_with_lps_vs_MDM_ap_neutr$ensg),]

rownames(MDM_ap_neut_with_lps_vs_MDM_ap_neutr)<-MDM_ap_neut_with_lps_vs_MDM_ap_neutr$ensg

MDM_ap_neut_with_lps_vs_MDM_ap_neutr<-MDM_ap_neut_with_lps_vs_MDM_ap_neutr[,1:(length(colnames(MDM_ap_neut_with_lps_vs_MDM_ap_neutr))-2)]

MDM_ap_neut_with_lps_vs_MDM_ap_neutr<-MDM_ap_neut_with_lps_vs_MDM_ap_neutr[order(rownames(MDM_ap_neut_with_lps_vs_MDM_ap_neutr), decreasing = F),]

dif_table_MDM_ap_neut_with_lps_vs_MDM_ap_neutr_filtered<-MDM_ap_neut_with_lps_vs_MDM_ap_neutr[which(MDM_ap_neut_with_lps_vs_MDM_ap_neutr$adj.P.Val<=0.01 & abs(MDM_ap_neut_with_lps_vs_MDM_ap_neutr$logFC)>=0.58),]

write.table(dif_table_MDM_ap_neut_with_lps_vs_MDM_ap_neutr_filtered, file = "DEG_results_GSE90010_mdm_ap_neut_with_lps_vs_mdm_ap_neutr_filtered_ensembl.csv", sep = "\t", col.names =TRUE)


write.table(MDM_ap_neut_with_lps_vs_MDM_ap_neutr, file = "DEG_results_GSE90010_mdm_ap_neut_with_lps_vs_mdm_ap_neutr_unfiltered_ensembl.csv", sep = "\t", col.names =TRUE)

########symbol

MDM_ap_neut_with_lps_vs_MDM_ap_neutr$ensg<-rownames(MDM_ap_neut_with_lps_vs_MDM_ap_neutr)

table_annot<- merge(MDM_ap_neut_with_lps_vs_MDM_ap_neutr, annotLookup, by.x="ensg", by.y="ensembl_gene_id")

table_annot_1<-table_annot[,c(2:8, 10)]
aggr_table<- table_annot_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(MDM_ap_neut_with_lps_vs_MDM_ap_neutr)[c(1:(length(colnames(MDM_ap_neut_with_lps_vs_MDM_ap_neutr))-1))]

dif_table_MDM_ap_neut_with_lps_vs_MDM_ap_neutr_filtered<-aggr_table[which(aggr_table$adj.P.Val<=0.01 & abs(aggr_table$logFC)>=0.58),]

write.table(dif_table_MDM_ap_neut_with_lps_vs_MDM_ap_neutr_filtered, file = "DEG_results_GSE90010_mdm_ap_neut_with_lps_vs_mdm_ap_neutr_filtered_symbol.csv", sep = "\t", col.names =TRUE)

write.table(aggr_table, file = "DEG_results_GSE90010_mdm_ap_neut_with_lps_vs_mdm_ap_neutr_unfiltered_symbol.csv", sep = "\t", col.names =TRUE)

############################################################################################################
############expr_mat####################################################

setwd("/Lung_data/Macrophage/Microarray/GSE90010/expression_data/expression_matrices")

expr_mat<-read.table("RAW_Expression_Matrix_Normalized_2023-11-01.txt", sep="\t")

ensembles<-rownames(expr_mat)

rownames_table <- sub("_.*", "", ensembles)

#rownames_table


expr_mat$ensg<-rownames_table

#expr_mat<-expr_mat[-grep("AFFX", expr_mat$ensg),]

rownames(expr_mat)<-expr_mat$ensg

expr_mat<-expr_mat[,1:(length(colnames(expr_mat))-1)]

expr_mat<-expr_mat[order(rownames(expr_mat), decreasing = F),]

write.table(expr_mat, file = "expression_matrix_GSE90010_ensembl.csv", sep = "\t", col.names =TRUE)

########symbol

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
  values = rownames(expr_mat),
  uniqueRows = TRUE)

expr_mat$ensg<-rownames(expr_mat)

table_annot<- merge(expr_mat, annotLookup, by.x="ensg", by.y="ensembl_gene_id")

table_annot_1<-table_annot[,c(2:25, 27)]
aggr_table<- table_annot_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table<- as.data.frame(aggr_table)
aggr_table  <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,c(2:length(colnames(aggr_table)))]
colnames(aggr_table)<-colnames(expr_mat)[c(1:(length(colnames(expr_mat))-1))]

write.table(aggr_table, file = "expression_matrix_GSE90010_symbol.csv", sep = "\t", col.names =TRUE)
