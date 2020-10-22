#!/bin/bash

set -e
set -u
set -o pipefail

if [ "$1" == "-h" ]
then
	echo "Usage: `basename $0` work_path index threads"
  	echo "Default is At"
	echo 
  	
	echo "This script is for ATAC_align"
  	echo -e "After QC\nuse the bowti2 to align\nuse the samtools to sort\nues the sambamba to mark dup\nuse the bedtools to filter"
  
  	exit 0
fi

# set up the software environment
module load multiqc/1.8

# set the work path
work_path=$1

# set directory with bowtie genome index
index=${2-~/reference/index/bowtie2/TAIR10/Athaliana}

# set threads
threads=${3-40}

# check if the filter.bed exists
# We next use the filter.bed to filter reads from Mt & Ch
if [ -f ./filter.bed ] || [ -f ./filter_plus.bed ]
then
	echo "there is existing filter.bed"
	echo "Please make sure you uses the filter or filter_plus"
else
	echo "Sorry there is not filter.bed"
	exit 1
fi

# output agrs name
echo "your work_path: ${work_path}"
echo "your index: ${index}"


# make the sample_for_align.txt
realpath ${work_path}/result/01_cleandata/*.gz | sed 's/_[12]\..*//' | uniq > ${work_path}/tmp/sample_for_align.txt

# make all of the output directories
mkdir -p ${work_path}/result/02_alignment
mkdir -p ${work_path}/result/03_filter_alignment
mkdir -p ${work_path}/logs/bowtie2_alignment
mkdir -p ${work_path}/logs/sambamba_markdup

# set up file names
alignment_out=${work_path}/result/02_alignment
filter_alignment_out=${work_path}/result/03_filter_alignment
bowtie2_log=${work_path}/logs/bowtie2_alignment
sambamba_log=${work_path}/logs/sambamba_markdup


# alignment using bowtie2
echo bowtie2 begins

cat ${work_path}/tmp/sample_for_align.txt | while read id;
do
	base=$(basename ${id})
	echo The sample name is ${base}

    bowtie2 -p ${threads} -x ${index} \
    -1 ${id}_1.clean.fq.gz \
    -2 ${id}_2.clean.fq.gz 2> ${bowtie2_log}/${base}.log \
    | samtools sort -@ 20 -O bam -o ${alignment_out}/${base}.sorted.bam - 
    
    samtools index ${alignment_out}/${base}.sorted.bam

    echo The sample ${base} is done
done

echo bowtie2 ends


# markdup with sambamba
echo sambamba begins

for i in ${alignment_out}/*.sorted.bam
do
	base=$(basename ${i} .sorted.bam)
	echo The sample name is ${base}

	sambamba markdup -t 5 ${i} ${alignment_out}/${base}.sorted.markdup.bam \
	2> ${sambamba_log}/${base}.log

	echo The sample ${base} is done
done

echo sambamba ends


# filter duplcation、multi-mappers、low_quality reads with samtools
echo filter begins

for j in ${alignment_out}/*.sorted.markdup.bam
do
	base=$(basename ${j} .sorted.markdup.bam)
	echo The sample name is ${base}

	samtools view -@ 20 -bF 1804 -q 20 ${j} -o ${filter_alignment_out}/${base}.flt.bam 

	samtools index -@ 20 ${filter_alignment_out}/${base}.flt.bam
	echo The sample ${base} is done
done

echo filter ends

# filter ChrM ChrCh with bedtools
echo bedtools filter organelle begins

for i in ${filter_alignment_out}/*.flt.bam
do
	base=$(basename ${i} .flt.bam)
	echo The sample name is ${base}

    bedtools intersect -abam ${i} -b filter.bed -v > ${filter_alignment_out}/${base}.rm_organelle.bam
    samtools index ${filter_alignment_out}/${base}.rm_organelle.bam
    
    echo The sample ${base} is done
done

echo bedtools filter is done

# rm tmp *.sorted.markdup.bam to save more space
echo begin removing tmp markdup bam and flt bam
rm -f ${alignment_out}/*.sorted.markdup.bam
rm -f ${alignment_out}/*.sorted.markdup.bam.bai

rm -f ${filter_alignment_out}/*.flt.bam
rm -f ${filter_alignment_out}/*.flt.bam.bai
echo removment complete

# multiqc 
multiqc -o ${bowtie2_log} ${bowtie2_log}


