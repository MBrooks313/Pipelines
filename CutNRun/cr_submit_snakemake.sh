#!/bin/sh
#SBATCH --job-name="brooksma"
#SBATCH --mail-type=FAIL
#SBATCH --output=logs/snakemake.%j.o
#SBATCH --cpus-per-task=1
#SBATCH --mem=1g
#SBATCH --time=1-00:00:00

#####################################
# This script is the Cut&Run submit script for a snakemake pipeline.
# This pipeline was created by Matthew J Brooks in March 2020
# The pipeline was based off the Henikoff lab pipeline.
# This pipeline adapted to run on HPCs running SLURM
# This requires the snakefile CutNRun_v3.1.py and cr_config.json
#####################################

# Load module
module load python/3.7

# Export variables
NOW=$(date +"%Y%m%d_%H-%M")
export NGS_PIPELINE="/data/brooksma/CutNRun/Mouse/src"
export WORK_DIR="/data/brooksma/CutNRun/Mouse/Analysis/"${NOW}
SNAKEFILE=$NGS_PIPELINE/CutNRun_v3.1.py

# Make result directories and change into result directory
mkdir -p ${WORK_DIR}/logs
cd $WORK_DIR

# Snakemake command
echo "Get ready for snakemake..." >> logs/snakemake.%j.o
snakemake\
	--directory $WORK_DIR \
	--snakefile $SNAKEFILE \
	--jobname '{rulename}.{jobid}' \
	--rerun-incomplete \
	--nolock \
	--verbose \
	-k -p \
	-j 3000 \
	--stats cr_pipeline_${NOW}.stats \
	--cluster "sbatch --mail-type=FAIL -o logs/{params.rulename}.%j.o {params.batch}" \
	>& cr_pipeline_${NOW}.log

# Summary
#  snakemake --directory $WORK_DIR --snakefile $SNAKEFILE --configfile $SAM_CONFIG --summary

## DRY Run with Print out the shell commands that will be executed
#  snakemake --directory $WORK_DIR --snakefile $SNAKEFILE --configfile $SAM_CONFIG --dryrun -p -r
# snakemake --directory $WORK_DIR --snakefile $SNAKEFILE --dryrun -p -r

#DAG
 # snakemake --directory $WORK_DIR --snakefile $SNAKEFILE  --dag | dot -Tpng > dag.png

#Rulegraph
#  snakemake --directory $WORK_DIR --snakefile $SNAKEFILE  -n --forceall --rulegraph | dot -Tpng > rulegraph.png

# Mail Rulegraph and DAG to self
#  echo DAG |mutt -s "DAG" -a dag.png -a rulegraph.png -- brooksma@mail.nih.gov
