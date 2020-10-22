#!/bin/bash

set -e
set -u
set -o pipefail

if [ "$1" == "-h" ]
then
    echo "Usage: `basename $0` Rep1.narrowPeak Rep2.narrowPeak Polled.narrowPeak NaiveOverlap_outputdir"

    exit 0
fi

# This script is from https://github.com/ENCODE-DCC/atac-seq-pipeline
# https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit

# This sctipr is ued to get the naiveOverlap peak
# bedtools v2.25.0

Rep1=$1
Rep2=$2
Polled=$3

output_dir=$4

base=$(basename ${Polled} .narrowPeak)

intersectBed -wo \
-a ${Polled} -b ${Rep1} | \
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | \
cut -f 1-10 | sort | uniq | \
intersectBed -wo \
-a stdin -b ${Rep2} | \
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | \
cut -f 1-10 | sort -k1,1 -k2,2n | uniq > ${output_dir}/${base}.NaiveOverlap.narrowPeak
