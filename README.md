# Curated-and-harmonised-transcriptomics-datasets-of-interstitial-lung-disease-patients

The code uploaded within this repository contains the processing code of the datatests of interstitial lung disease patients retrieved from the NCBI Gene Expression Omnibus and European Nucleotide Archive (ENA) 
repositories. Overall, we retrieved 30 transcriptomics datasets, produced through both DNA microarrays and RNA-Sequencing technologies, along with relative meta-data tables. After data collection, all the datasets 
underwent meta-data curation and harmonisation, data quality check, and pre-processing with standardised procedures. Beside, data model was created in order to homogenise the phenotypic data enabling the comparability 
across the datasets. In addition, for every study included in this data release we provide tables of differentially expressed genes and co-expression network models of IPF samples and healthy counterparts deriving from 
both RNA-seq and DNA microarrays.

###Explanations for the codes###

1. fastqc.sh

	-makes the fastqc quality control reports for the fastq files

2. cutadapt.sh and cutadapt_single.sh

	-Adapter trimming and read filtering. cutadapt is for paired end data and cutadapt_single for single end and ion torrent data

3. hisat2_launch.sh and launch_hisat2_single.sh

	-Hisat2 alignment. hisat2_launch for paired end data and launch_hisat2_single for single end data

4. uniq_sort.sh
  
   -Creates a BAM file with the unique reads and sorts the unique reads. 

7. featurecounts.R

	-Transcript quantification

8. Contents in expression_matrices_and_DEG_result -folder
	
	-Microarray
		-Contains the probe annotation codes for the microarray datasets. Annotations are for ensembl gene id:s and gene Symbols. 
		The R programming required substantial manual effort, largely due to the inherent characteristics of the public data. As a precautionary 
		measure to ensure the accuracy of each step in the analytical process, the analysis included significant redundancy.

	-RNA-seq
		-Contains the low counts filtering for the counts matrix and annotations for gene symbols. Contains also the gene expression matrix normalization and
		differential gene expression analysis with deseq2. The R programming required substantial manual effort, largely due to the inherent characteristics 
		of the public data. As a precautionary  measure to ensure the accuracy of each step in the analytical process, the analysis included significant redundancy.

9. Contents of combined_matrices -folder
	-Contains the code for combining the expression matrices for ipf and healthy samples for each cell type and tissue and RNA-seq and Microarray.

10. multi_studies_adjust.R
   
	-Pamr batch adjustment code for the combined matrices

11. build_networks.R
    
    -Network inference code for the batch adjusted combined matrices
