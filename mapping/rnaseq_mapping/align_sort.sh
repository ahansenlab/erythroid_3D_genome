#!/bin/bash

while getopts "f:n:s:o:" opt; do
  case $opt in
    f)
      fq_dir="$OPTARG"
      ;;
    n)
      fq_name="$OPTARG"
      ;;
    s)
      sample_name="$OPTARG"
      ;;
    o)
      outdir="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

# Shift off the parsed options so $1 points to the first non-flag argument

genome_fa="/mnt/md0/DataRepository/genomes/hg38_analysisSet/hg38.analysisSet.fa"
threads=8
genomedir="/mnt/md1/varshini/star_hg38"
fq1="${fq_dir}/${fq_name}/${fq_name}_1.fq.gz"
fq2="${fq_dir}/${fq_name}/${fq_name}_2.fq.gz"
gtf="/mnt/md0/varshini/Gene_Annotations/gencode.v45.annotation.gtf"
echo "aligning ${fq1} ${fq2}..."

# star-align
if [ ! -d "$genomedir" ]; then
STAR --genomeDir "$genomedir" \
--genomeFastaFiles "$genome_fa" \
--sjdbGTFfile "$gtf" \
--runMode genomeGenerate \
--runThreadN "$threads"
fi 

#STAR --genomeDir "$genomedir" \
#--runThreadN "$threads" \
#--readFilesIn "$fq1" "$fq2" \
#--outFileNamePrefix "${outdir}/${sample_name}" \
#--outSAMtype BAM Unsorted \
#--outSAMattributes Standard \
#--quantMode GeneCounts \
#--outTmpDir "${outdir}/${fq_name}_temp" \
#--readFilesCommand "gunzip -c"
echo "${outdir}/${sample_name}.sorted.bam" | cat -v

#sort 
sambamba sort --tmpdir "${outdir}/${fq_name}_temp_sambamba" -t "$threads" -m 30G -o "${outdir}/${sample_name}.sorted.bam" "${outdir}/${sample_name}Aligned.out.bam"

# deduplicate
if [ ! -f "${outdir}/${sample_name}.sorted.rmdup.bam" ]; then
sambamba markdup --tmpdir "${outdir}/${fq_name}_temp_sambamba_dup" --overflow-list-size 600000 -r -t "$threads" "${outdir}/${sample_name}.sorted.bam" "${outdir}/${sample_name}.sorted.rmdup.bam"
sambamba index -t "$threads" "${outdir}/${sample_name}.sorted.rmdup.bam"
else
echo "${outdir}/${sample_name}.sorted.rmdup.bam exists, continuing..."
fi

# index 
sambamba index -t "$threads" "${outdir}/${sample_name}.sorted.rmdup.bam"

echo "done! at ${outdir}/${sample_name}"
