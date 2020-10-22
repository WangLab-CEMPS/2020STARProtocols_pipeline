#!/bin/bash

set -e
set -u
set -o pipefail

if [ "$1" == "-h" ]
then
  	echo "Usage: `basename $0` work_path"
  
  	echo
  	echo "This script is for ATAC_QC"
  	echo "After fastp, use the fastqc and multiqc to check the quality"
  
  	exit 0
fi

# set up the software environment
module load FastQC/0.11.7

# set the work path
work_path=$1

# make the tree file and test
mkdir -p ${work_path}/result ${work_path}/logs ${work_path}/tmp

file_number=`ls ${work_path} -F | grep "/" | wc -l`

if [  ${file_number} -eq 4 ]
then
  	echo "everything is OK. Let's do it"
else
  	echo "Sorry, something is wrong"
  	exit 1
fi

# output agrs name
echo "your work_path: ${work_path}"

# make the sample_for_QC.txt
realpath ${work_path}/rawdata/*.gz | sed 's/_[12]\..*//' | uniq > ${work_path}/tmp/sample_for_QC.txt

# make the output file
mkdir -p ${work_path}/result/01_cleandata
mkdir -p ${work_path}/logs/fastqc_v1
mkdir -p ${work_path}/logs/fastp

# set up file names
cleandata_out=${work_path}/result/01_cleandata
fastqc_v1_log=${work_path}/logs/fastqc_v1
fastp_log=${work_path}/logs/fastp

# fastp
echo fastp begins

cat ${work_path}/tmp/sample_for_QC.txt | while read id; 
do
    base=$(basename ${id})
    echo Sample name is ${base}

    fastp --thread 16 \
    -a CTGTCTCTTATACACATCT \
    -i ${id}_1.fq.gz \
    -I ${id}_2.fq.gz \
    -o ${cleandata_out}/${base}_1.clean.fq.gz \
    -O ${cleandata_out}/${base}_2.clean.fq.gz \
    -j ${fastp_log}/${base}.fastp.json \
    -h ${fastp_log}/${base}.fastp.html \
    2> ${fastp_log}/${base}.logs

    echo Sample ${base} is done
done

echo fastp ends

# fastqc
fastqc -q -t 20 ${cleandata_out}/*.fq.gz -o ${fastqc_v1_log}


# multiqc
module load multiqc/1.8
multiqc -o ${fastqc_v1_log} ${fastqc_v1_log}
multiqc -o ${fastp_log} ${fastp_log}

