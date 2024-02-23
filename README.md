# Curated-and-harmonised-transcriptomics-datasets-of-interstitial-lung-disease-patients

The code uploaded within this repository contains the processing code of the datatests of interstitial lung disease patients retrieved from the NCBI Gene Expression Omnibus and European Nucleotide Archive (ENA) 
repositories. Overall, we retrieved 30 transcriptomics datasets, produced through both DNA microarrays and RNA-Sequencing technologies, along with relative meta-data tables. After data collection, all the datasets 
underwent meta-data curation and harmonisation, data quality check, and pre-processing with standardised procedures. Beside, data model was created in order to homogenise the phenotypic data enabling the comparability 
across the datasets. In addition, for every study included in this data release we provide tables of differentially expressed genes and co-expression network models of IPF samples and healthy counterparts deriving from 
both RNA-seq and DNA microarrays.


###Explanations for the codes###

The raw data files of RNA-seq (.fastq) files have been downloaded from European Nucleotide Archive (ENA) and the raw data files DNA microarray (.CEL and .txt files) were retrieved from GEO using GEOquery. 
The meta-data curation was performed with ESPERANTO software: https://github.com/fhaive/esperanto. 
The normalization and DEG-analysis of the microarray datasets was performed with eUTOPIA softwre https://github.com/Greco-Lab/eUTOPIA.

The pipeline goes in the following order below.

1. fastqc.sh

	-makes the fastqc quality control reports for the fastq files

2. cutadapt.sh and cutadapt_single.sh

	-Adapter trimming and read filtering. cutadapt is for paired end data and cutadapt_single for single end and ion torrent data

3. hisat2_launch.sh and launch_hisat2_single.sh

	-Hisat2 alignment. hisat2_launch for paired end data and launch_hisat2_single for single end data

4. uniq_sort.sh
  
   -Creates a BAM file with the unique reads and sorts the unique reads. 

5. featurecounts.R

	-Transcript quantification

6. Contents in expression_matrices_and_DEG_result -folder
	
	-Microarray
		-Contains the probe annotation codes for the microarray datasets. Annotations are for ensembl gene id:s and gene Symbols. 
		The R programming required substantial manual effort, largely due to the inherent characteristics of the public data. As a precautionary 
		measure to ensure the accuracy of each step in the analytical process, the analysis included significant redundancy.

	-RNA-seq
		-Contains the low counts filtering for the counts matrix and annotations for gene symbols. Contains also the gene expression matrix normalization and
		differential gene expression analysis with deseq2. The R programming required substantial manual effort, largely due to the inherent characteristics 
		of the public data. As a precautionary  measure to ensure the accuracy of each step in the analytical process, the analysis included significant redundancy.

7. Contents of combined_matrices -folder
   
    -Contains the code for combining the expression matrices for ipf and healthy samples for each cell type and tissue and RNA-seq and Microarray.

8. multi_studies_adjust.R
   
	-Pamr batch adjustment code for the combined matrices

9. build_networks.R
    
    -Network inference code for the batch adjusted combined matrices

All the preprocessed transcriptomics data, along with harmonised meta-data and networks, were submitted to Zenodo: https://doi.org/10.5281/zenodo.10692129.
