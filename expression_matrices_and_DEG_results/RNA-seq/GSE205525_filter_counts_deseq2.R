#NO RAW DATA AVAILABLE FOR THIS DATASET. RAW DATA ANNOTATED TO GENE SYMBOLS


setwd("/Lung_data/Epithelial/RNA-seq/GSE205525/expression_data/expression_matrices")


raw_matrix <- read.csv("GSE205525_coculture_counts_data.csv", sep=",")

raw_matrix<-raw_matrix[,2:length(colnames(raw_matrix))]

# Sample vector
target_id <- raw_matrix$target_id

ensembl <- c()

# Loop through the vector and modify the strings
for (i in 1:length(target_id)) {
  ensembl[i] <- sub("\\..*", "", target_id[i])
}

raw_matrix$ensembl<-ensembl

raw_matrix<-raw_matrix[,2:length(colnames(raw_matrix))]


aggr_table<- raw_matrix %>% group_by(ensembl) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
#aggr_table<-aggr_table[-which(aggr_table$ensembl_gene_id == ""), ]
rownames(aggr_table)<-aggr_table$ensembl
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(raw_matrix)[c(1:(length(colnames(raw_matrix))-1))]




write.table(aggr_table, file = "raw_matrix_GSE205525_ensembl.csv", sep = "\t", col.names =TRUE)


###################################################

GSE205525_raw_data<-read.table("raw_matrix_GSE205525_ensembl.csv", sep = "\t")
library(xlsx)

metadata_GSE205525 <- read.xlsx("/Lung_data/Epithelial/RNA-seq/GSE205525/phenodata/GSE205525_curated.xlsx", sheetIndex = 1)


###metadata and data samples in different order, lets reorder

GSE205525_raw_data<-GSE205525_raw_data[,as.character(metadata_GSE205525$title)]

conditions<-metadata_GSE205525$disease


#########################Filter the low counts#############################################

filter_low_counts <- function(counts.matrix, conditions, method = "cpm", normalized=FALSE, depth=NULL, cpm=1, p.adj = "fdr"){
  
  if(is.null(counts.matrix)){stop("Error: please provide a numeric count matrix!")}
  if(is.null(conditions)){stop("Error: please provide a factor or a vector indicating the conditions!")}
  if(!method %in% c("cpm", "wilcoxon", "proportion")) {stop("Error: Please type in one of valid methods!")}
  
  
  if (method=="cpm"){
    filtered.counts = NOISeq::filtered.data(counts.matrix, factor = conditions, norm = normalized, method = 1, cv.cutoff = 100, cpm = cpm, p.adj = p.adj)
    
  }else if(method=="wilcoxon"){
    filtered.counts = NOISeq::filtered.data(counts.matrix, factor = conditions, norm = normalized, method = 2, cv.cutoff = 100, p.adj = p.adj)
    
  }else if(method=="proportion"){
    
    if(is.null(depth)){stop("Error: indicate a numeric vector indicating per sample library depths")}
    if(!class(depth)=="numeric"){stop("Error: please provide the depth argument with a numeric vector!")}
    ### Compute librarary depth
    
    
    filtered.counts = NOISeq::filtered.data(counts.matrix, factor = conditions, norm = normalized, depth = depth, method = 3, cv.cutoff = 100, cpm = cpm, p.adj = p.adj)
  }
  return(filtered.counts)
}
filtered_data <- filter_low_counts(counts.matrix =GSE205525_raw_data, conditions = conditions, method = "proportion", normalized = FALSE, p.adj = "fdr", depth = as.numeric(apply(GSE205525_raw_data, 2, sum)))

write.table(filtered_data, file = "Filtered_expression_matrix_ensembl_GSE205525.csv", sep = "\t", col.names =TRUE)

############################################

############Biomart Symbol

require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

annotLookup <- getBM(
  mart=mart,
  attributes=c(
    "ensembl_gene_id",
    "external_gene_name"),
  filter = "ensembl_gene_id",
  values = rownames(filtered_data),
  uniqueRows = TRUE)



filtered_data$ensembl<-rownames(filtered_data)


table_with_gene_id <- merge(filtered_data, annotLookup, by.x="ensembl", by.y="ensembl_gene_id")

table_with_gene_id_1<-table_with_gene_id[,c(2:length(colnames(table_with_gene_id)))] # REMOVE excess columns
aggr_table<- table_with_gene_id_1 %>% group_by(external_gene_name) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
aggr_table<-aggr_table[-which(aggr_table$external_gene_name == ""), ]
rownames(aggr_table)<-aggr_table$external_gene_name
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(filtered_data)[c(1:(length(colnames(filtered_data))-1))]
write.table(aggr_table, file = "Filtered_expression_matrix_symbol_GSE205525.csv", sep = "\t", col.names =TRUE)



###PCA before normalisation###

setwd("/Lung_data/supplements/RNA-seq_pca_plots")

table_transpose<-as.data.frame(t(aggr_table))
df_pca <- prcomp(table_transpose, center = TRUE, scale. = TRUE)
p<-ggplot(as.data.frame(df_pca$x), aes(x=PC1, y=PC2, color=conditions, label = colnames(aggr_table)))
p <- p+geom_point()+ggtitle("Before normalization")
p

ggsave(file="GSE205525_pca_before_normalization_condition.svg", width=10, height=8)

timepoint<-metadata_GSE205525$timepoint_days

table_transpose<-as.data.frame(t(aggr_table))
df_pca <- prcomp(table_transpose, center = TRUE, scale. = TRUE)
p<-ggplot(as.data.frame(df_pca$x), aes(x=PC1, y=PC2, color=timepoint, label = colnames(aggr_table)))
p <- p+geom_point()+ggtitle("Before normalization")
p
ggsave(file="GSE205525_pca_before_normalization_timepoint.svg", width=10, height=8)


########################################################################################

### Differential expression analysis IPF ALL Symbol ###

setwd("/Lung_data/Epithelial/RNA-seq/GSE205525/expression_data/expression_matrices")
expr_mat<-read.table("Filtered_expression_matrix_symbol_GSE205525.csv", sep = "\t")
metadata_GSE205525 <- read.xlsx("/Lung_data/Epithelial/RNA-seq/GSE205525/phenodata/GSE205525_curated.xlsx", sheetIndex = 1)

#######variables##########################

conditions<-metadata_GSE205525$disease

sex<-metadata_GSE205525$sex

fibro_status<-metadata_GSE205525$fibroblast_disease_status

#sample_source<-metadata_GSE205525$sample_source #fibro_status and sample_source linear combinations. Model doesn't approve

time<-metadata_GSE205525$timepoint_days

###############################################################

table_1<-sapply(expr_mat,as.integer)

rownames(table_1)<-rownames(expr_mat)


colData <- data.frame(condition=as.vector(conditions),
                       sex=as.vector(sex), fibro_status=as.vector(fibro_status),time=as.vector(time))

rownames(colData) <- colnames(table_1)




ddsMat <- DESeq2::DESeqDataSetFromMatrix(countData = table_1,
                                         colData = colData,
                                         design = ~sex+fibro_status+time+condition)

dds <- DESeq2::DESeq(ddsMat)
total_norm_counts <- DESeq2::counts(dds, normalized=TRUE)
write.table(total_norm_counts, file = "Normalized_expression_matrix_symbol_GSE205525.txt", quote = FALSE, sep = "\t", row.names = TRUE, col=NA)

####PCA after normalization#########################3

setwd("/Lung_data/supplements/RNA-seq_pca_plots")

table_transpose<-as.data.frame(t(total_norm_counts))
df_pca <- prcomp(table_transpose, center = TRUE, scale. = TRUE)
p<-ggplot(as.data.frame(df_pca$x), aes(x=PC1, y=PC2, color=conditions, label = colnames(total_norm_counts)))
p <- p+geom_point()+ggtitle("After normalization")
p

ggsave(file="GSE205525_pca_after_normalization_condition.svg", width=10, height=8)

timepoint<-metadata_GSE205525$timepoint_days

table_transpose<-as.data.frame(t(total_norm_counts))
df_pca <- prcomp(table_transpose, center = TRUE, scale. = TRUE)
p<-ggplot(as.data.frame(df_pca$x), aes(x=PC1, y=PC2, color=timepoint, label = colnames(total_norm_counts)))
p <- p+geom_point()+ggtitle("After normalization")
p
ggsave(file="GSE205525_pca_after_normalization_timepoint.svg", width=10, height=8)


##

setwd("/Lung_data/Epithelial/RNA-seq/GSE205525/expression_data/DEG_results")

res1 <- DESeq2::results(dds, contrast = c("condition", "IPF", "healthy"), pAdjustMethod = "fdr", independentFiltering = FALSE)
write.table(res1, file = "DEG_results_DESeq2_IPF_ALL_vs_healthy_ALL_unfiltered_symbol_GSE205525.txt", row.names = TRUE, sep = "\t", quote = FALSE)
res1_adj <- res1[which(res1$padj<=0.01 & abs(res1$log2FoldChange)>=0.58),]
print(dim(res1_adj))
write.table(res1_adj, file = "DEG_results_DESeq2_IPF_ALL_vs_healthy_ALL_filtered_symbol_GSE205525.txt", row.names = TRUE, sep = "\t", quote = FALSE)

### Differential expression analysis IPF ALL Ensembl ###

setwd("/Lung_data/Epithelial/RNA-seq/GSE205525/expression_data/expression_matrices")
expr_mat<-read.table("Filtered_expression_matrix_ensembl_GSE205525.csv", sep = "\t")
metadata_GSE205525 <- read.xlsx("/Lung_data/Epithelial/RNA-seq/GSE205525/phenodata/GSE205525_curated.xlsx", sheetIndex = 1)

#######variables##########################

conditions<-metadata_GSE205525$disease

sex<-metadata_GSE205525$sex

fibro_status<-metadata_GSE205525$fibroblast_disease_status

#sample_source<-metadata_GSE205525$sample_source fibro_status and sample_source linear combinations. Model doesnt approve

time<-metadata_GSE205525$timepoint_days

###############################################################

table_1<-sapply(expr_mat,as.integer)

rownames(table_1)<-rownames(expr_mat)


colData <- data.frame(condition=as.vector(conditions),
                      sex=as.vector(sex), fibro_status=as.vector(fibro_status),time=as.vector(time))

rownames(colData) <- colnames(table_1)




ddsMat <- DESeq2::DESeqDataSetFromMatrix(countData = table_1,
                                         colData = colData,
                                         design = ~sex+fibro_status+time+condition)

dds <- DESeq2::DESeq(ddsMat)
total_norm_counts <- DESeq2::counts(dds, normalized=TRUE)
write.table(total_norm_counts, file = "Normalized_expression_matrix_ensembl_GSE205525.txt", quote = FALSE, sep = "\t", row.names = TRUE, col=NA)

setwd("/Lung_data/Epithelial/RNA-seq/GSE205525/expression_data/DEG_results")
res1 <- DESeq2::results(dds, contrast = c("condition", "IPF", "healthy"), pAdjustMethod = "fdr", independentFiltering = FALSE)
write.table(res1, file = "DEG_results_DESeq2_IPF_ALL_vs_healthy_ALL_unfiltered_ensembl_GSE205525.txt", row.names = TRUE, sep = "\t", quote = FALSE)
res1_adj <- res1[which(res1$padj<=0.01 & abs(res1$log2FoldChange)>=0.58),]
print(dim(res1_adj))
write.table(res1_adj, file = "DEG_results_DESeq2_IPF_ALL_vs_healthy_ALL_filtered_ensembl_GSE205525.txt", row.names = TRUE, sep = "\t", quote = FALSE)

#################################################################################################
#################################################################################################
### Differential expression analysis_fibroblast_status ###
##################################################################################################
################################Symbol##################################################################

setwd("/Lung_data/Epithelial/RNA-seq/GSE205525/expression_data/expression_matrices")
expr_mat<-read.table("Filtered_expression_matrix_symbol_GSE205525.csv", sep = "\t")
metadata_GSE205525 <- read.xlsx("/Lung_data/Epithelial/RNA-seq/GSE205525/phenodata/GSE205525_curated.xlsx", sheetIndex = 1)

#######variables##########################

conditions<-metadata_GSE205525$sample_description

sex<-metadata_GSE205525$sex

#fibro_status<-metadata_GSE205525$fibroblast_disease_status not relevant here

#sample_source<-metadata_GSE205525$sample_source #sample_source and time and sample_source and condition linear combinations. Model doesnt approve

time<-metadata_GSE205525$timepoint_days

###############################################################

table_1<-sapply(expr_mat,as.integer)

rownames(table_1)<-rownames(expr_mat)


colData <- data.frame(condition=as.vector(conditions),
                      sex=as.vector(sex),time=as.vector(time))

rownames(colData) <- colnames(table_1)




ddsMat <- DESeq2::DESeqDataSetFromMatrix(countData = table_1,
                                         colData = colData,
                                         design = ~sex+time+condition)

dds <- DESeq2::DESeq(ddsMat)

setwd("/Lung_data/Epithelial/RNA-seq/GSE205525/expression_data/DEG_results")


res1 <- DESeq2::results(dds, contrast = c("condition", "IPF_no_fibroblast", "healthy_no_fibroblast"), pAdjustMethod = "fdr", independentFiltering = FALSE)
write.table(res1, file = "DEG_results_DESeq2_IPF_NO_FIBRO_vs_healthy_NO_FIBRO_unfiltered_symbol_GSE205525.txt", row.names = TRUE, sep = "\t", quote = FALSE)
res1_adj <- res1[which(res1$padj<=0.01 & abs(res1$log2FoldChange)>=0.58),]
print(dim(res1_adj))
write.table(res1_adj, file = "DEG_results_DESeq2_IPF_NO_FIBRO_vs_healthy_NO_FIBRO_filtered_symbol_GSE205525.txt", row.names = TRUE, sep = "\t", quote = FALSE)

res1 <- DESeq2::results(dds, contrast = c("condition", "healthy_IPF_fibroblast", "healthy_no_fibroblast"), pAdjustMethod = "fdr", independentFiltering = FALSE)
write.table(res1, file = "DEG_results_DESeq2_healthy_IPF_FIBRO_vs_healthy_NO_FIBRO_unfiltered_symbol_GSE205525.txt", row.names = TRUE, sep = "\t", quote = FALSE)
res1_adj <- res1[which(res1$padj<=0.01 & abs(res1$log2FoldChange)>=0.58),]
print(dim(res1_adj))
write.table(res1_adj, file = "DEG_results_DESeq2_healthy_IPF_FIBRO_vs_healthy_NO_FIBRO_filtered_symbol_GSE205525.txt", row.names = TRUE, sep = "\t", quote = FALSE)

res1 <- DESeq2::results(dds, contrast = c("condition", "healthy_IPF_fibroblast", "healthy_healthy_fibroblast"), pAdjustMethod = "fdr", independentFiltering = FALSE)
write.table(res1, file = "DEG_results_DESeq2_healthy_IPF_FIBRO_vs_healthy_healthy_FIBRO_unfiltered_symbol_GSE205525.txt", row.names = TRUE, sep = "\t", quote = FALSE)
res1_adj <- res1[which(res1$padj<=0.01 & abs(res1$log2FoldChange)>=0.58),]
print(dim(res1_adj))
write.table(res1_adj, file = "DEG_results_DESeq2_healthy_IPF_FIBRO_vs_healthy_healthy_FIBRO_filtered_symbol_GSE205525.txt", row.names = TRUE, sep = "\t", quote = FALSE)

res1 <- DESeq2::results(dds, contrast = c("condition", "IPF_no_fibroblast", "healthy_healthy_fibroblast"), pAdjustMethod = "fdr", independentFiltering = FALSE)
write.table(res1, file = "DEG_results_DESeq2_IPF_NO_FIBRO_vs_healthy_healthy_FIBRO_unfiltered_symbol_GSE205525.txt", row.names = TRUE, sep = "\t", quote = FALSE)
res1_adj <- res1[which(res1$padj<=0.01 & abs(res1$log2FoldChange)>=0.58),]
print(dim(res1_adj))
write.table(res1_adj, file = "DEG_results_DESeq2_IPF_NO_FIBRO_vs_healthy_healthy_FIBRO_filtered_symbol_GSE205525.txt", row.names = TRUE, sep = "\t", quote = FALSE)

res1 <- DESeq2::results(dds, contrast = c("condition", "IPF_no_fibroblast", "healthy_IPF_fibroblast"), pAdjustMethod = "fdr", independentFiltering = FALSE)
write.table(res1, file = "DEG_results_DESeq2_IPF_NO_FIBRO_vs_healthy_IPF_FIBRO_unfiltered_symbol_GSE205525.txt", row.names = TRUE, sep = "\t", quote = FALSE)
res1_adj <- res1[which(res1$padj<=0.01 & abs(res1$log2FoldChange)>=0.58),]
print(dim(res1_adj))
write.table(res1_adj, file = "DEG_results_DESeq2_IPF_NO_FIBRO_vs_healthy_IPF_FIBRO_filtered_symbol_GSE205525.txt", row.names = TRUE, sep = "\t", quote = FALSE)


### Differential expression analysis fibro status Ensembl ###

setwd("/Lung_data/Epithelial/RNA-seq/GSE205525/expression_data/expression_matrices")
expr_mat<-read.table("Filtered_expression_matrix_ensembl_GSE205525.csv", sep = "\t")
metadata_GSE205525 <- read.xlsx("/Lung_data/Epithelial/RNA-seq/GSE205525/phenodata/GSE205525_curated.xlsx", sheetIndex = 1)

#######variables##########################

conditions<-metadata_GSE205525$sample_description

sex<-metadata_GSE205525$sex

#fibro_status<-metadata_GSE205525$fibroblast_disease_status not relevant here

#sample_source<-metadata_GSE205525$sample_source #sample_source and time and sample_source and condition linear combinations. Model doesnt approve

time<-metadata_GSE205525$timepoint_days


###############################################################

table_1<-sapply(expr_mat,as.integer)

rownames(table_1)<-rownames(expr_mat)


colData <- data.frame(condition=as.vector(conditions),
                      sex=as.vector(sex),time=as.vector(time))

rownames(colData) <- colnames(table_1)




ddsMat <- DESeq2::DESeqDataSetFromMatrix(countData = table_1,
                                         colData = colData,
                                         design = ~sex+time+condition)

dds <- DESeq2::DESeq(ddsMat)

setwd("/Lung_data/Epithelial/RNA-seq/GSE205525/expression_data/DEG_results")


res1 <- DESeq2::results(dds, contrast = c("condition", "IPF_no_fibroblast", "healthy_no_fibroblast"), pAdjustMethod = "fdr", independentFiltering = FALSE)
write.table(res1, file = "DEG_results_DESeq2_IPF_NO_FIBRO_vs_healthy_NO_FIBRO_unfiltered_ensembl_GSE205525.txt", row.names = TRUE, sep = "\t", quote = FALSE)
res1_adj <- res1[which(res1$padj<=0.01 & abs(res1$log2FoldChange)>=0.58),]
print(dim(res1_adj))
write.table(res1_adj, file = "DEG_results_DESeq2_IPF_NO_FIBRO_vs_healthy_NO_FIBRO_filtered_ensembl_GSE205525.txt", row.names = TRUE, sep = "\t", quote = FALSE)

res1 <- DESeq2::results(dds, contrast = c("condition", "healthy_IPF_fibroblast", "healthy_no_fibroblast"), pAdjustMethod = "fdr", independentFiltering = FALSE)
write.table(res1, file = "DEG_results_DESeq2_healthy_IPF_FIBRO_vs_healthy_NO_FIBRO_unfiltered_ensembl_GSE205525.txt", row.names = TRUE, sep = "\t", quote = FALSE)
res1_adj <- res1[which(res1$padj<=0.01 & abs(res1$log2FoldChange)>=0.58),]
print(dim(res1_adj))
write.table(res1_adj, file = "DEG_results_DESeq2_healthy_IPF_FIBRO_vs_healthy_NO_FIBRO_filtered_ensembl_GSE205525.txt", row.names = TRUE, sep = "\t", quote = FALSE)

res1 <- DESeq2::results(dds, contrast = c("condition", "healthy_IPF_fibroblast", "healthy_healthy_fibroblast"), pAdjustMethod = "fdr", independentFiltering = FALSE)
write.table(res1, file = "DEG_results_DESeq2_healthy_IPF_FIBRO_vs_healthy_healthy_FIBRO_unfiltered_ensembl_GSE205525.txt", row.names = TRUE, sep = "\t", quote = FALSE)
res1_adj <- res1[which(res1$padj<=0.01 & abs(res1$log2FoldChange)>=0.58),]
print(dim(res1_adj))
write.table(res1_adj, file = "DEG_results_DESeq2_healthy_IPF_FIBRO_vs_healthy_healthy_FIBRO_filtered_ensembl_GSE205525.txt", row.names = TRUE, sep = "\t", quote = FALSE)

res1 <- DESeq2::results(dds, contrast = c("condition", "IPF_no_fibroblast", "healthy_healthy_fibroblast"), pAdjustMethod = "fdr", independentFiltering = FALSE)
write.table(res1, file = "DEG_results_DESeq2_IPF_NO_FIBRO_vs_healthy_healthy_FIBRO_unfiltered_ensembl_GSE205525.txt", row.names = TRUE, sep = "\t", quote = FALSE)
res1_adj <- res1[which(res1$padj<=0.01 & abs(res1$log2FoldChange)>=0.58),]
print(dim(res1_adj))
write.table(res1_adj, file = "DEG_results_DESeq2_IPF_NO_FIBRO_vs_healthy_healthy_FIBRO_filtered_ensembl_GSE205525.txt", row.names = TRUE, sep = "\t", quote = FALSE)

res1 <- DESeq2::results(dds, contrast = c("condition", "IPF_no_fibroblast", "healthy_IPF_fibroblast"), pAdjustMethod = "fdr", independentFiltering = FALSE)
write.table(res1, file = "DEG_results_DESeq2_IPF_NO_FIBRO_vs_healthy_IPF_FIBRO_unfiltered_ensembl_GSE205525.txt", row.names = TRUE, sep = "\t", quote = FALSE)
res1_adj <- res1[which(res1$padj<=0.01 & abs(res1$log2FoldChange)>=0.58),]
print(dim(res1_adj))
write.table(res1_adj, file = "DEG_results_DESeq2_IPF_NO_FIBRO_vs_healthy_IPF_FIBRO_filtered_ensembl_GSE205525.txt", row.names = TRUE, sep = "\t", quote = FALSE)


#################################################################################################
#################################################################################################
