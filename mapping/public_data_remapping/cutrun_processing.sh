#!/bin/bash 

# process a single file
srr_id="$1"
sample="$2"
index="${3:-/mnt/md0/DataRepository/genomes/hg38_analysisSet/hg38.analysisSet.fa.gz}"
threads="${4:-8}"

fastq1="${srr_id}_1.fastq"
fastq2="${srr_id}_2.fastq"

# align and remove multimappers
bowtie2 --dovetail --threads $threads -x $index -1 $fastq1 -2 $fastq2 | grep -v XS: - | samtools view -bh -F 1024 -f 2 -b > "bam/${sample}.bam"

# rmdup 
#sambamba sort -t "$threads" -m "10G" -o "bam/${sample}.sorted.bam" "bam/${sample}.bam"
#sambamba markdup -r -t "$threads" "bam/${sample}.sorted.bam" "bam/${sample}.rmdup.bam"
#sambamba index "bam/${sample}.sorted.bam"

samtools collate -@ 4 -O -u "bam/${sample}.bam" | samtools fixmate -@ 4 -m -u - - | samtools sort -@ 4 -u - | samtools markdup -@ 4 - "bam/${sample}.markdup.bam"
samtools view -F 1024 -f 2 -b "bam/${sample}.markdup.bam" > "bam/${sample}.rmdup.bam"

#prepare for seacr: bamtobed/awk
sambamba sort -n -t "$threads" -m "10G" -o "bam/${sample}.rmdup.sorted.bam" "bam/${sample}.rmdup.bam"
sambamba index "bam/${sample}.rumdup.sorted.bam"

bedtools bamtobed -i "bam/${sample}.rmdup.sorted.bam" -bedpe > "processed/${sample}.bedpe"
awk '$1==$4 && $1!="." && $6-$2 < 1000 {print $0}' "processed/${sample}.bedpe" | cut -f 1-6 | sort -k1,1 -k2,2n -k3,3n > "processed/${sample}.filt.bedgraph"












