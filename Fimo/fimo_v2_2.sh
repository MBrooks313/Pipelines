#!/bin/sh

###############################
# This function run the fimo function in Meme Suite for a bed file and matrix named in the cml.
# This version v2.2 is adapted for HPC running SLURM.
# Written by Matthew J Brooks on Sept 10th, 2020
###############################

# Run statement
if [[ $# -lt 2 ]] ; then
    printf "\nUSAGE: fimo_v2_*.sh <peaks.bed> <motifs.meme>\n\n"
    printf "You need cml args for the bed and motif. - MJB\n\n"
    exit 1
fi

# Get cml variables
bed=$1
mot=$2


#-----------USER MODIFIED VARIABLES--------------------------------------------#
# Module versions for Biowulf
bedtools_ver=2.29.2
meme_ver=5.1.0

# Motif
motif_dir=/data/brooksma/Index/TRANSFAC_Pro_data/TFP_2017.3_data/dat/Meme
motif=${motif_dir}/$mot

# Bed file and location
bed_dir=/data/brooksma/CutNRun/Mouse/Analysis/20200309_16-29/Beds/idr/Final
base=${bed%%.*}
fasta=${base}.fa

# Reference fasta
ref=/data/brooksma/Index/Mouse/ENS/v98/Mus_musculus.GRCm38.dna.primary_assembly.fa

# Working directory
work_dir=/data/brooksma/CutNRun/Mouse/Analysis/20200309_16-29/Meme/fimo/Analysis_200910
#------------------------------------------------------------------------------#


# Date
now=$(date +'%y%m%d')

# Load modules
module load bedtools/$bedtools_ver || exit 1
module load meme/$meme_ver || exit 1


# Make fasta from bed file
if [ -f "${bed_dir}/$fasta" ]; then
    echo "Fasta exists moving on ..."
else
    echo "Sort and modifying bed file..."
    bedtools sort -i ${bed_dir}/${bed} | \
    tr -d 'chr' > ${bed_dir}/${base}_tmp.bed;
    echo "Getting the fasta ..."
    bedtools getfasta -fi $ref -bed ${bed_dir}/${base}_tmp.bed -fo ${bed_dir}/${fasta};
    rm ${bed_dir}/${base}_tmp.bed
fi


# USAGE: fimo [options] <motif file> <sequence file>
fimo -oc ${work_dir}/${base}_${now} $motif ${bed_dir}/${fasta}
