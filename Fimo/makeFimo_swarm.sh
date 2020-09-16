#!/bin/sh

###############################
# This makes the swarm file for FIMO analysis using the fimo_v2_2.sh function.
# Written by Matthew J Brooks on Sept 10th, 2020
###############################

# Variables
bed_dir=/data/brooksma/CutNRun/Mouse/Analysis/20200309_16-29/Beds/idr/Final
outfile=fimo_200910.swarm

# Write run statment for swarm file
echo "#swarm -f "${outfile}" -g 8 --time 4:00:00 --logdir logs" > $outfile

# Swarm commands 
for i in `ls ${bed_dir}/*narrowPeak`
do
    base=${i##*/}
    echo "bash fimo_v2_2.sh ${base} matrix_Mm.meme" >> $outfile
done;
