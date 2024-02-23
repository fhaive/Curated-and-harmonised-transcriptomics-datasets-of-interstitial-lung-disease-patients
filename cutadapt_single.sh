FASTQ_FILES=$(ls *.fastq.gz |cut -d "." -f1)

for i in $FASTQ_FILES; do ~/.local/bin/cutadapt -a AGATCGGAAGAG -q 20 -m 60 -j 5 -o TRIMMED_${i}.fastq.gz ${i}.fastq.gz; done
