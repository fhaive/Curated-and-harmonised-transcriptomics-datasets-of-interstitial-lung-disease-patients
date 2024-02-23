#NO RAW DATA AVAILABLE FOR THIS DATASET. AVAILABLE DATA ANNOTATED TO GENE SYMBOLS


setwd("/Lung_data/Biopsy/RNA-seq/GSE150910/expression_data/expression_matrices")

raw_matrix <- read.csv("GSE150910_gene-level_count_file.csv", sep=",")

rownames(raw_matrix) <- raw_matrix$symbol

raw_matrix<-raw_matrix[,2:length(colnames(raw_matrix))]

#change the colnames to SRRxxxxxxxx

#column_names<-c()
#for (i in 1:length(colnames(raw_matrix))) {
#new_name<-substring(colnames(raw_matrix)[i], 9, 18)
#column_names<-c(column_names, new_name)
#}
#colnames(raw_matrix)<-column_names

###Read in the metadata to check the conditions

library(xlsx)


metadata <- read.xlsx("/Lung_data/Biopsy/RNA-seq/GSE150910/phenodata/GSE150910_curated.xlsx", sheetIndex = 1)

conditions<-(metadata$disease)



###metadata and data samples in different order, lets reorder

raw_matrix<-raw_matrix[,as.character(metadata$title)]


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
filtered_data <- filter_low_counts(counts.matrix =raw_matrix, conditions = conditions, method = "proportion", normalized = FALSE, p.adj = "fdr", depth = as.numeric(apply(raw_matrix, 2, sum)))

write.table(filtered_data, file = "Filtered_expression_matrix_symbol_GSE150910.csv", sep = "\t", col.names =TRUE)

############################################

###PCA before normalisation###


setwd("/Lung_data/supplements/RNA-seq_pca_plots")

table_transpose<-as.data.frame(t(filtered_data))
df_pca <- prcomp(table_transpose, center = TRUE, scale. = TRUE)
p<-ggplot(as.data.frame(df_pca$x), aes(x=PC1, y=PC2, color=conditions, label = colnames(filtered_data)))
p <- p+geom_point()+ggtitle("Before normalization")
p

ggsave(file="GSE150910_pca_before_normalization.svg", width=10, height=8)




############Biomart Ensg

require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)

annotLookup <- getBM(
  mart=mart,
  attributes=c(
    "ensembl_gene_id",
    "external_gene_name"),
  filter = "external_gene_name",
  values = rownames(filtered_data),
  uniqueRows = TRUE)



ensembl_ids = annotLookup$ensembl_gene_id
ensembl_ids = as.vector(ensembl_ids)

gene_symbols<-rownames(filtered_data)

filtered_data$gen_symb<-rownames(filtered_data)

table_with_gene_id <- merge(filtered_data, annotLookup, by.x="gen_symb", by.y="external_gene_name")

table_with_gene_id_1<-table_with_gene_id[,c(2:length(colnames(table_with_gene_id)))] # REMOVE excess columns
aggr_table<- table_with_gene_id_1 %>% group_by(ensembl_gene_id) %>% dplyr::summarise_all(.funs = c(median = "median"))
aggr_table <- as.data.frame(aggr_table)
aggr_table <- na.omit(aggr_table)
#aggr_table<-aggr_table[-which(aggr_table$ensembl_gene_id == ""), ]
rownames(aggr_table)<-aggr_table$ensembl_gene_id
aggr_table<-aggr_table[,2:length(colnames(aggr_table))]
colnames(aggr_table)<-colnames(filtered_data)[c(1:(length(colnames(filtered_data))-1))]
write.table(aggr_table, file = "Filtered_expression_matrix_ensembl_GSE150910.csv", sep = "\t", col.names =TRUE)



### Differential expression analysis symbol ###
#####SYMBOL

setwd("/Lung_data/Biopsy/RNA-seq/GSE150910/expression_data/expression_matrices")
expr_mat<-read.table("Filtered_expression_matrix_symbol_GSE150910.csv", sep = "\t")
metadata <- read.xlsx("/Lung_data/Biopsy/RNA-seq/GSE150910/phenodata/GSE150910_curated.xlsx", sheetIndex = 1)

#######Variables##########################

conditions<-metadata$disease

sex<-metadata$sex #NA

muc5b_genotype<-metadata$muc5b_genotype

#smoking_status<-metadata$smoking_status ### contains NA values

plate<-metadata$plate

antigen<-metadata$antigen_identified

batch<-metadata$batch

batch<-gsub("-", "_", batch)

#race<-metadata$race ## Contains NA values

immunosuppressant<-metadata$immunosupressant ## Contains NA

lane<-metadata$lane

#age<-metadata$age ##Contains NA values

tissue_source<-metadata$tissue_source

###############################################################

table_1<-sapply(expr_mat,as.integer)

rownames(table_1)<-rownames(expr_mat)


colData <- data.frame(condition=as.vector(conditions),  muc5b_genotype=as.vector(muc5b_genotype),
                       plate=as.vector(plate),batch=as.vector(batch), lane=as.vector(lane),  
                        tissue_source=as.vector(tissue_source), antigen=as.vector(antigen))     

rownames(colData) <- colnames(table_1)




ddsMat <- DESeq2::DESeqDataSetFromMatrix(countData = table_1,
                                         colData = colData,
                                         design = ~tissue_source+lane+muc5b_genotype+
                                           plate+condition) #This model works


dds <- DESeq2::DESeq(ddsMat)
total_norm_counts <- DESeq2::counts(dds, normalized=TRUE)
write.table(total_norm_counts, file = "Normalized_expression_matrix_symbol_GSE150910.txt", quote = FALSE, sep = "\t", row.names = TRUE, col=NA)

########PCA after normalization
setwd("/Lung_data/supplements/RNA-seq_pca_plots")

table_transpose<-as.data.frame(t(total_norm_counts))
df_pca <- prcomp(table_transpose, center = TRUE, scale. = TRUE)
p<-ggplot(as.data.frame(df_pca$x), aes(x=PC1, y=PC2, color=conditions))
p <- p+geom_point()+ggtitle("After normalization")
p

ggsave(file="GSE150910_pca_after_normalization.svg", width=10, height=8)
######################################################################

setwd("/Lung_data/Biopsy/RNA-seq/GSE150910/expression_data/DEG_results")
res1 <- DESeq2::results(dds, contrast = c("condition", "IPF", "healthy"), pAdjustMethod = "fdr", independentFiltering = FALSE)
write.table(res1, file = "DEG_results_DESeq2_IPF_vs_healthy_unfiltered_symbol_GSE150910.txt", row.names = TRUE, sep = "\t", quote = FALSE)
res1_adj <- res1[which(res1$padj<=0.01 & abs(res1$log2FoldChange)>=0.58),]
print(dim(res1_adj))
write.table(res1_adj, file = "DEG_results_DESeq2_IPF_vs_healthy_filtered_symbol_GSE150910.txt", row.names = TRUE, sep = "\t", quote = FALSE)

res1 <- DESeq2::results(dds, contrast = c("condition", "hypersensitivity_pneumonitis", "healthy"), pAdjustMethod = "fdr", independentFiltering = FALSE)
write.table(res1, file = "DEG_results_DESeq2_hypersensitivity_pneumonitis_vs_healthy_unfiltered_symbol_GSE150910.txt", row.names = TRUE, sep = "\t", quote = FALSE)
res1_adj <- res1[which(res1$padj<=0.01 & abs(res1$log2FoldChange)>=0.58),]
print(dim(res1_adj))
write.table(res1_adj, file = "DEG_results_DESeq2_hypersensitivity_pneumonitis_vs_healthy_filtered_symbol_GSE150910.txt", row.names = TRUE, sep = "\t", quote = FALSE)

##########################################################################################
#########################################################################################
### Differential expression analysis ensembl IPF ###

setwd("/Lung_data/Biopsy/RNA-seq/GSE150910/expression_data/expression_matrices")
expr_mat<-read.table("Filtered_expression_matrix_ensembl_GSE150910.csv", sep = "\t")
metadata <- read.xlsx("/Lung_data/Biopsy/RNA-seq/GSE150910/phenodata/GSE150910_curated.xlsx", sheetIndex = 1)


#Variables


conditions<-metadata$disease

sex<-metadata$sex #NA

muc5b_genotype<-metadata$muc5b_genotype

#smoking_status<-metadata$smoking_status ### contains NA values

plate<-metadata$plate

antigen<-metadata$antigen_identified

batch<-metadata$batch

batch<-gsub("-", "_", batch)

#race<-metadata$race ## Contains NA values

#immunosuppressant<-metadata$immunosupressant ## ALL NA

lane<-metadata$lane

#age<-metadata$age ##Contains NA values

tissue_source<-metadata$tissue_source

###############################################################

table_1<-sapply(expr_mat,as.integer)

rownames(table_1)<-rownames(expr_mat)


colData <- data.frame(condition=as.vector(conditions),  muc5b_genotype=as.vector(muc5b_genotype),
                      plate=as.vector(plate),batch=as.vector(batch), lane=as.vector(lane),  
                      tissue_source=as.vector(tissue_source), antigen=as.vector(antigen))

rownames(colData) <- colnames(table_1)




ddsMat <- DESeq2::DESeqDataSetFromMatrix(countData = table_1,
                                         colData = colData,
                                         design = ~tissue_source+lane+muc5b_genotype+
                                           plate+condition) #this model works



dds <- DESeq2::DESeq(ddsMat)
total_norm_counts <- DESeq2::counts(dds, normalized=TRUE)
write.table(total_norm_counts, file = "Normalized_expression_matrix_ensembl_GSE150910.txt", quote = FALSE, sep = "\t", row.names = TRUE, col=NA)


setwd("/Lung_data/Biopsy/RNA-seq/GSE150910/expression_data/DEG_results")
res1 <- DESeq2::results(dds, contrast = c("condition", "IPF", "healthy"), pAdjustMethod = "fdr", independentFiltering = FALSE)
write.table(res1, file = "DEG_results_DESeq2_IPF_vs_healthy_unfiltered_ensembl_GSE150910.txt", row.names = TRUE, sep = "\t", quote = FALSE)
res1_adj <- res1[which(res1$padj<=0.01 & abs(res1$log2FoldChange)>=0.58),]
print(dim(res1_adj))
write.table(res1_adj, file = "DEG_results_DESeq2_IPF_vs_healthy_filtered_ensembl_GSE150910.txt", row.names = TRUE, sep = "\t", quote = FALSE)

res1 <- DESeq2::results(dds, contrast = c("condition", "hypersensitivity_pneumonitis", "healthy"), pAdjustMethod = "fdr", independentFiltering = FALSE)
write.table(res1, file = "DEG_results_DESeq2_hypersensitivity_pneumonitis_vs_healthy_unfiltered_ensembl_GSE150910.txt", row.names = TRUE, sep = "\t", quote = FALSE)
res1_adj <- res1[which(res1$padj<=0.01 & abs(res1$log2FoldChange)>=0.58),]
print(dim(res1_adj))
write.table(res1_adj, file = "DEG_results_DESeq2_hypersensitivity_pneumonitis_vs_healthy_filtered_ensembl_GSE150910.txt", row.names = TRUE, sep = "\t", quote = FALSE)


