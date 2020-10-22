#!/bin/bash

set -e
set -u
set -o pipefail

if [ "$1" == "-h" ]
then
    echo "Usage: `basename $0` work_path thread"
    echo 

    echo "This script is for ATAC_bamtobw"
    echo -e "After align & filter, use the deeptools to transfrom bam to bw"

    exit 0
fi

# set up the software environment
module load deeptools/2.0

# set the work path
work_path=$1

# set threads
threads=${2-20}

# make the output file
mkdir -p ${work_path}/result/05_bamtobw

# set up file names
input=${work_path}/result/03_filter_alignment
bw_out=${work_path}/result/05_bamtobw

# change bam to bw
ls ${input}/*.bam | while read id;
do
	base=$(basename ${id} .bam)

	bamCoverage -b ${id} \
	--binSize 10 --numberOfProcessors ${threads} \
	-o ${bw_out}/${base}.bw --normalizeUsing BPM
done

