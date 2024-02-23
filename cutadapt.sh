FASTQ_FILES_R1=$(ls *_1.fastq.gz |cut -d "_" -f1)
FASTQ_FILES_R2=$(ls *_2.fastq.gz |cut -d "_" -f1)

for i in $FASTQ_FILES_R1; do ~/.local/bin/cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG  -q 20 -m 60 -j 2 -o TRIMMED_${i}_1.fastq.gz -p TRIMMED_${i}_2.fastq.gz ${i}_1.fastq.gz ${i}_2.fastq.gz; done
