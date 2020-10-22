#!/bin/bash

set -e
set -u
set -o pipefail

if [ "$1" == "-h" ]
then
    echo "Usage: `basename $0` work_path bed_path thread" 
	echo "default is At"
	echo "please please pay attention the species !!!!!!!!!!!"
    echo 

    echo "This script is for Chip_plot"
    echo -e "After transfrom, use the deeptools for plot heatmap、profiler or something else"

    exit 0
fi

# set up the software environment
module load deeptools/2.0

# set the work path
work_path=$1

# set the threads
threads=${3-20}

# set the bed file
bed_file=${2-/home/sgd/reference/annoation/Athaliana/Araport11/Araport11_whole_gene_for_deeptools.bed}

# make the output file
mkdir -p ${work_path}/result/06_plot

# set up filenames
input=${work_path}/result/05_bamtobw
plot_out=${work_path}/result/06_plot

# computerMatirx(This should change)
## scale-regions
computeMatrix scale-regions \
-S ${input}/*.bw \
-R ${bed_file} \
--regionBodyLength 2000 \
--beforeRegionStartLength 3000 --afterRegionStartLength 3000 \
--skipZeros --numberOfProcessors ${threads} \
-o ${plot_out}/matrix_scale_region.scale.gz

## reference-point
computeMatrix reference-point \
-S ${input}/*.bw \
-R ${bed_file} \
--referencePoint TSS \
-b 3000 -a 3000 \
--skipZeros --numberOfProcessors ${threads} \
-o ${plot_out}/matrix_reference_point.reference.gz


# computer multiBigwigSummary
multiBigwigSummary bins -b ${input}/*.bw \
--numberOfProcessors ${threads} \
-o ${plot_out}/multibw_results.npz

# plot

## peak plot
plotProfile -m ${plot_out}/matrix_scale_region.scale.gz -out ${plot_out}/scale_region.pdf --perGroup
plotProfile -m ${plot_out}/matrix_scale_region.scale.gz -out ${plot_out}/scale_region_persample.pdf --numPlotsPerRow 4
plotProfile -m ${plot_out}/matrix_reference_point.reference.gz -out ${plot_out}/reference_point_region.pdf --perGroup
plotProfile -m ${plot_out}/matrix_reference_point.reference.gz -out ${plot_out}/reference_point_region_persample.pdf --numPlotsPerRow 4

rm ${plot_out}/matrix_reference_point.reference.gz -rf
rm ${plot_out}/matrix_scale_region.scale.gz -rf

## correlation plot

### plot scatterplot
plotCorrelation -in ${plot_out}/multibw_results.npz \
--corMethod spearman --skipZeros \
--whatToPlot scatterplot \
--plotTitle "Spearman Correlation" \
--removeOutliers \
--plotFile ${plot_out}/correlation_spearman_bwscore_scatterplot.pdf

### plot heatmap
plotCorrelation -in ${plot_out}/multibw_results.npz \
--corMethod spearman --skipZeros \
--whatToPlot heatmap \
--plotTitle "Spearman Correlation" \
--removeOutliers \
--plotNumbers \
--plotFile ${plot_out}/correlation_spearman_bwscore_heatmapplot.pdf


rm ${plot_out}/multibw_results.npz -rf
