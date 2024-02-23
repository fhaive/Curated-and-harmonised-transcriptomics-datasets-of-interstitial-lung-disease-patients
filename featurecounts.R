setwd("/nasdata/sinkala/expressiondata/rnaseq/GSE124685_ion")
files <- list.files(pattern = "\\.bam")

#isPairedEnd is FALSE or TRUE depending on the dataset

counts<-Rsubread::featureCounts(files, isGTFAnnotationFile = TRUE, annot.ext="/nasdata/RNA_Seq/references/Ensembl_v108_hsapiens/gtf/Homo_sapiens.GRCh38.108.gtf", GTF.attrType="gene_id", isPairedEnd=FALSE)

