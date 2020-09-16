#!/bin/sh

###############################
# This makes the swarm file for converting the FIMO output to bed files using the Fimo2bed_v1.1.R function.
# Written by Matthew J Brooks on Sept 11th, 2020
###############################

# Variables
wd=/data/brooksma/CutNRun/Mouse/Analysis/20200309_16-29/Meme/fimo/Analysis_200910/
mat=/data/brooksma/Index/TRANSFAC_Pro_data/TFP_2017.3_data/dat/matrixMAPgene_Mm_v2.txt
outfile=fimo2bed_v1_1.swarm


# Get paths/files needed
array=()
while IFS=  read -r -d $'\0'; do
    array+=("$REPLY")
done < <(find ${wd} -name "fimo.tsv" -print0)

# Write run statment for swarm file
echo "#swarm -f "${outfile}" -g 8 --time 1:00:00 --gres=lscratch:100 --logdir logs" > $outfile

# Swarm commands
for i in "${array[@]}";
do
echo "
Rscript Fimo2bed_v1_1.R \
${i} \
${mat}" >> $outfile
done;
