HISAT2_INDEXES=/path/indexes/grch38
FASTQ_FILES_R1=$(ls *1.fastq.gz | cut -d "." -f1 | sed 's/_1$//')
FASTQ_FILES_R2=$(ls *2.fastq.gz | cut -d "." -f1 | sed 's/_2$//')
for i in $FASTQ_FILES_R1; do
/path/hisat2-2.2.1/hisat2 -q -p 10 -x $HISAT2_INDEXES/genome -1 ${i}_1.fastq.gz -2 ${i}_2.fastq.gz | samtools view -Sbh > ${i}.bam;
done
