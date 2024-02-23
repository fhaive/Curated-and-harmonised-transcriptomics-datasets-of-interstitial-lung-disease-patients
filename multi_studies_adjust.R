rm(list=ls())

#### Required libraries ####
suppressMessages(library(esc))
suppressMessages(library(metafor))
suppressMessages(library(metap))
suppressMessages(library(TopKLists))
suppressMessages(library(RankProd))
suppressMessages(library(matrixStats))
suppressMessages(library(igraph))
suppressMessages(library(TopKLists))
suppressMessages(library(minet))
suppressMessages(library(foreach))
suppressMessages(library(parallel))
suppressMessages(library(doParallel))
suppressMessages(library(ggplot2))
suppressMessages(library(plyr))

######## Module 1 - Meta-analysis section ########
##################################################


#' Mean-adjustes transcriptomics data by batch (wrapper of pamr.batchadjust from the pamr CRAN package)
#' 
#' @importFrom pamr pamr.batchadjust
#' 
#' @param expr_mat A dataframe with genes on the rows and samples in the columns.
#' @param samples_label A factor of samples labels in the same order as the samples in expr_mat columns
#' @param batch_labels A factor of labels indicating the batches (the studies where the samples are coming from) 
#' @return A batch-adjusted expression matrix of the same dimension of expr_mat
#' @examples
#' \dontrun {
#' calc_effect_size_rank(meta_dataframe)
#' }
#' @export
multi_studies_adjust <- function(expr_mat, samples_label, batch_labels){
  mylist <- list(x=as.matrix(expr_mat), y=as.factor(samples_label), batchlabels=as.factor(batch_labels))
  adjusted_mat <- pamr::pamr.batchadjust(data = mylist)
  
  table_transpose<-as.data.frame(t(expr_mat))
  df_pca <- prcomp(table_transpose, center = TRUE, scale. = TRUE)
  p<-ggplot(as.data.frame(df_pca$x), aes(x=PC1, y=PC2, color=batch_labels))
  p <- p+geom_point(size=3)+guides(color = guide_legend(title = "Batch labels"))+
    ggtitle("Before batch correction")+ 
    theme_bw()+
    theme(legend.key.size = unit(5, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=16, face="bold"), #change legend title font size
          legend.text = element_text(size=14),
          plot.title = element_text(size = 20, face = "bold"),
          axis.text.x = element_text(color = "grey20", size=14),
          axis.text.y = element_text(color = "grey20", size = 14),
          axis.title=element_text(size=14,face="bold"))
  
  
  table_transpose2<-as.data.frame(t(adjusted_mat$x))
  df_pca2 <- prcomp(table_transpose2, center = TRUE, scale. = TRUE)
  p2<-ggplot(as.data.frame(df_pca2$x), aes(x=PC1, y=PC2, color=batch_labels))
  p2 <- p2+geom_point(size=3)+ guides(color = guide_legend(title = "Batch labels"))+
    ggtitle("After batch correction") +
    theme_bw()+
    theme(legend.key.size = unit(5, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text( size=16, face="bold"), #change legend title font size
          legend.text = element_text(size=14),
          plot.title = element_text(size = 20, face = "bold"),
          axis.text.x = element_text(color = "grey20", size=14),
          axis.text.y = element_text(color = "grey20", size = 14),
          axis.title=element_text(size=14,face="bold"))
  
  require(gridExtra)
  grid.arrange(p, p2, ncol=2)
  
  return(adjusted_mat)
}
