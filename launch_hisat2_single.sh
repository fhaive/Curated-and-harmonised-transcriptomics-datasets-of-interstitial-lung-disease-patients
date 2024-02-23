HISAT2_INDEXES=/path/indexes/grch38
FASTQ_FILES_R1=$(ls *.fastq.gz |cut -d "." -f1)

for i in $FASTQ_FILES_R1; do /path/hisat2-2.2.1/hisat2 -q -p 10 -x $HISAT2_INDEXES/genome -U ${i}.fastq.gz | samtools view -Sbh > ${i}.bam; done
