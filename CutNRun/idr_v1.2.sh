#!/bin/bash

#####################################
# This script performs IDR analysis.
# Based on Qunhua Li and Peter Bickel’s group
# This script generates a swarm file for running on an HPC running SLURM.
# Written by Matthew J Brooks on July 20th, 2020.
#####################################



# Set list variables
tf_list=(Crx Nrl)
tp_list=(E14 P2 P4 P10 P28)
rep_list=(0 1 2 3)

# Set analysis variables
wd=/data/brooksma/CutNRun/Mouse/src
IDR_THRESH=0.05
idr_ver=2.0.3
bed_dir=/data/brooksma/CutNRun/Mouse/Analysis/20200309_16-29/Beds

# Goto the directory with MACS peaks
cd $bed_dir;

# ----------------------#
# Write to swarm script

# submit script for swarm file
echo '#swarm -f idr_v1.1.swarm -g 8 -t 4 --time 8:00:00 --logdir logs' > ${wd}/idr_v1.1.swarm

# loop for each transcription factor
for i in ${tf_list[@]};
do

    # Loop for each time point
    for j in ${tp_list[@]};
    do

        # Get samples and IDs for tf-timepoint group
        samps=(`ls *${j}-*${i}*rmblklst.narrowPeak`)
        id=${i}_${j}

        # Print TF and TP info
        echo "#${id}" >> ${wd}/idr_v1.1.swarm

        ## Rep comparison 1
        rep=Rep${rep_list[0]}vRep${rep_list[1]}
        out_dir=${bed_dir}/idr/${id}/${rep}

        echo "#${rep}" >> ${wd}/idr_v1.1.swarm
        echo "
        mkdir -p $out_dir; \\
        module load idr/$idr_ver || exit 1; \\
        idr --samples ${bed_dir}/${samps[0]} ${bed_dir}/${samps[1]}  \\
        --input-file-type narrowPeak \\
        --output-file ${out_dir}/${id}.${rep}.idr.narrowPeak \\
        --output-file-type narrowPeak \\
        --log-output-file ${out_dir}/${id}.$rep.log \\
        --soft-idr-threshold $IDR_THRESH \\
        --plot \\
        --use-best-multisummit-IDR; \\
        IDR_THRESH_TRANSFORMED=$(awk -v p=${IDR_THRESH} 'BEGIN{print -log(p)/log(10)}'); \\
        awk 'BEGIN{OFS=\"\\t\"} \$12>='\"\${IDR_THRESH_TRANSFORMED}\"' {print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}' ${out_dir}/${id}.${rep}.idr.narrowPeak | sort | uniq | sort -k7n,7n > ${out_dir}/${id}.${rep}.idr05.narrowPeak; \\
        " >> ${wd}/idr_v1.1.swarm

        ## Rep comparison 2
        rep=Rep${rep_list[0]}vRep${rep_list[2]}
        out_dir=${bed_dir}/idr/${id}/${rep}

        echo "#${rep}" >> ${wd}/idr_v1.1.swarm
        echo "
        mkdir -p $out_dir; \\
        module load idr/$idr_ver || exit 1; \\
        idr --samples ${bed_dir}/${samps[0]} ${bed_dir}/${samps[2]}  \\
        --input-file-type narrowPeak \\
        --output-file ${out_dir}/${id}.${rep}.idr.narrowPeak \\
        --output-file-type narrowPeak \\
        --log-output-file ${out_dir}/${id}.$rep.log \\
        --soft-idr-threshold $IDR_THRESH \\
        --plot \\
        --use-best-multisummit-IDR; \\
        IDR_THRESH_TRANSFORMED=$(awk -v p=${IDR_THRESH} 'BEGIN{print -log(p)/log(10)}'); \\
        awk 'BEGIN{OFS=\"\\t\"} \$12>='\"\${IDR_THRESH_TRANSFORMED}\"' {print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}' ${out_dir}/${id}.${rep}.idr.narrowPeak | sort | uniq | sort -k7n,7n > ${out_dir}/${id}.${rep}.idr05.narrowPeak; \\
        " >> ${wd}/idr_v1.1.swarm

        ## Rep comparison 3
        rep=Rep${rep_list[0]}vRep${rep_list[3]}
        out_dir=${bed_dir}/idr/${id}/${rep}

        echo "#${rep}" >> ${wd}/idr_v1.1.swarm
        echo "
        mkdir -p $out_dir; \\
        module load idr/$idr_ver || exit 1; \\
        idr --samples ${bed_dir}/${samps[0]} ${bed_dir}/${samps[3]}  \\
        --input-file-type narrowPeak \\
        --output-file ${out_dir}/${id}.${rep}.idr.narrowPeak \\
        --output-file-type narrowPeak \\
        --log-output-file ${out_dir}/${id}.$rep.log \\
        --soft-idr-threshold $IDR_THRESH \\
        --plot \\
        --use-best-multisummit-IDR; \\
        IDR_THRESH_TRANSFORMED=$(awk -v p=${IDR_THRESH} 'BEGIN{print -log(p)/log(10)}'); \\
        awk 'BEGIN{OFS=\"\\t\"} \$12>='\"\${IDR_THRESH_TRANSFORMED}\"' {print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}' ${out_dir}/${id}.${rep}.idr.narrowPeak | sort | uniq | sort -k7n,7n > ${out_dir}/${id}.${rep}.idr05.narrowPeak; \\
        " >> ${wd}/idr_v1.1.swarm

        ## Rep comparison 4
        rep=Rep${rep_list[1]}vRep${rep_list[2]}
        out_dir=${bed_dir}/idr/${id}/${rep}

        echo "#${rep}" >> ${wd}/idr_v1.1.swarm
        echo "
        mkdir -p $out_dir; \\
        module load idr/$idr_ver || exit 1; \\
        idr --samples ${bed_dir}/${samps[1]} ${bed_dir}/${samps[2]}  \\
        --input-file-type narrowPeak \\
        --output-file ${out_dir}/${id}.${rep}.idr.narrowPeak \\
        --output-file-type narrowPeak \\
        --log-output-file ${out_dir}/${id}.$rep.log \\
        --soft-idr-threshold $IDR_THRESH \\
        --plot \\
        --use-best-multisummit-IDR; \\
        IDR_THRESH_TRANSFORMED=$(awk -v p=${IDR_THRESH} 'BEGIN{print -log(p)/log(10)}'); \\
        awk 'BEGIN{OFS=\"\\t\"} \$12>='\"\${IDR_THRESH_TRANSFORMED}\"' {print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}' ${out_dir}/${id}.${rep}.idr.narrowPeak | sort | uniq | sort -k7n,7n > ${out_dir}/${id}.${rep}.idr05.narrowPeak; \\
        " >> ${wd}/idr_v1.1.swarm

        ## Rep comparison 5
        rep=Rep${rep_list[1]}vRep${rep_list[3]}
        out_dir=${bed_dir}/idr/${id}/${rep}

        echo "#${rep}" >> ${wd}/idr_v1.1.swarm
        echo "
        mkdir -p $out_dir; \\
        module load idr/$idr_ver || exit 1; \\
        idr --samples ${bed_dir}/${samps[1]} ${bed_dir}/${samps[3]}  \\
        --input-file-type narrowPeak \\
        --output-file ${out_dir}/${id}.${rep}.idr.narrowPeak \\
        --output-file-type narrowPeak \\
        --log-output-file ${out_dir}/${id}.$rep.log \\
        --soft-idr-threshold $IDR_THRESH \\
        --plot \\
        --use-best-multisummit-IDR; \\
        IDR_THRESH_TRANSFORMED=$(awk -v p=${IDR_THRESH} 'BEGIN{print -log(p)/log(10)}'); \\
        awk 'BEGIN{OFS=\"\\t\"} \$12>='\"\${IDR_THRESH_TRANSFORMED}\"' {print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}' ${out_dir}/${id}.${rep}.idr.narrowPeak | sort | uniq | sort -k7n,7n > ${out_dir}/${id}.${rep}.idr05.narrowPeak; \\
        " >> ${wd}/idr_v1.1.swarm

        ## Rep comparison 6
        rep=Rep${rep_list[2]}vRep${rep_list[3]}
        out_dir=${bed_dir}/idr/${id}/${rep}

        echo "#${rep}" >> ${wd}/idr_v1.1.swarm
        echo "
        mkdir -p $out_dir; \\
        module load idr/$idr_ver || exit 1; \\
        idr --samples ${bed_dir}/${samps[2]} ${bed_dir}/${samps[3]}  \\
        --input-file-type narrowPeak \\
        --output-file ${out_dir}/${id}.${rep}.idr.narrowPeak \\
        --output-file-type narrowPeak \\
        --log-output-file ${out_dir}/${id}.$rep.log \\
        --soft-idr-threshold $IDR_THRESH \\
        --plot \\
        --use-best-multisummit-IDR; \\
        IDR_THRESH_TRANSFORMED=$(awk -v p=${IDR_THRESH} 'BEGIN{print -log(p)/log(10)}'); \\
        awk 'BEGIN{OFS=\"\\t\"} \$12>='\"\${IDR_THRESH_TRANSFORMED}\"' {print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}' ${out_dir}/${id}.${rep}.idr.narrowPeak | sort | uniq | sort -k7n,7n > ${out_dir}/${id}.${rep}.idr05.narrowPeak; \\
        " >> ${wd}/idr_v1.1.swarm

    done;
done;
