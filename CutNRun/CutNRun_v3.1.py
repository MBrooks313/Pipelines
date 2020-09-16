#######################################
# This is a Cut&Run analysis snakemake script.
# Written my Matthew J. Brooks on March 9th, 2020
# This runs on an HPC running SLURM
#######################################



#######################################
# Import config file and modules needed
#######################################

# Import modules
import glob
import os
import json

# Snakemake Base location
try:
	NGS_PIPELINE=os.environ['NGS_PIPELINE']
except KeyError:
	print("I can not locate your NGS_PIPELINE directory.")
	pass


# Import configs
configfile: NGS_PIPELINE + "/cr_config.json"
FASTQS = config["fastq_dir"]


########################################
# Import sample names from the FQ folder
########################################

SAMPLES = [os.path.basename(fname).split('.')[0] for fname in glob.glob(FASTQS + '/*.R1.fastq.gz')]


#############################################################
# List of directories needed and end point files for analysis
#############################################################

BAMS = expand("Bams/{sample}.bam", sample=SAMPLES)
RMBAMS = expand("Bams/{sample}_rmdup.bam", sample=SAMPLES)
RAWPEAKS = expand("Beds/{sample}_peaks.narrowPeak", sample=SAMPLES)
FINPEAKS = expand("Beds/{sample}_peaks_rmblklst.narrowPeak", sample=SAMPLES)
BIGWIG = expand("Beds/{sample}_peaks_bpm.bw", sample=SAMPLES)
QCDEEP = ["QCstats/Summary_Fingerprint.tsv"]

##############################
# Snakemake rules for analysis
##############################

localrules: all

rule all:
        input:  BAMS + RMBAMS + RAWPEAKS + FINPEAKS + BIGWIG + QCDEEP
        params:
                batch = config["job_all"]


rule trimAlign:
        """
        This trims fastqs, aligns, filters for MAPQ > 20, sorts, and indexes bam file.
        """
        input:
                R1 = FASTQS + "/{sample}.R1.fastq.gz",
                R2 = FASTQS + "/{sample}.R2.fastq.gz",
                adapter = config["trimmoadapt"]
        output:
                forward = temp("trim/{sample}.R1.fastq.gz"),
                reverse = temp("trim/{sample}.R2.fastq.gz"),
                for_un = temp("trim/unpaired_{sample}.R1.fastq.gz"),
                rev_un = temp("trim/unpaired_{sample}.R2.fastq.gz"),
                flag_raw = "QCstats/{sample}_Flagstat_raw.txt",
                flag_mapq = "QCstats/{sample}_Flagstat_mapq.txt",
                bam = "Bams/{sample}.bam",
                bai = "Bams/{sample}.bam.bai"
        log:    "logs/bowtie.{sample}.log"
        version: config["bowtie"]
        params:
                rulename = "trimAlign",
                batch = config["job_trimAlign"],
                trimmo = config["trimmo"],
                samtools = config["samtools"],
                ref = config["bt2idx"]
        shell:  """
        ######################
        #### Load modules ####
        module load trimmomatic/{params.trimmo} || exit 1
        module load bowtie/{version} || exit 1
        module load samtools/{params.samtools} || exit 1
        ####################
        #### Trim fastq ####
        java -jar $TRIMMOJAR PE -threads ${{SLURM_CPUS_ON_NODE}} \
        {input.R1} {input.R2} \
        {output.forward} {output.for_un} \
        {output.reverse} {output.rev_un} \
        ILLUMINACLIP:{input.adapter}:2:15:4:4:true TAILCROP:6 LEADING:20 TRAILING:20 MINLEN:25 2>> {log}
        ###################################################
        #### Align, filter for MAPQ > 20, and sort bam ####
        bowtie2 -p ${{SLURM_CPUS_ON_NODE}} -t \
        --end-to-end --dovetail -I 10 -X 700 \
        --very-sensitive --no-unal --no-mixed --no-discordant -q --phred33  \
        -x {params.ref} \
        -1 {output.forward} \
        -2 {output.reverse} 2>> {log} | \
        samtools view -bS - > /lscratch/${{SLURM_JOBID}}/raw.bam
        samtools view -bSq 20 /lscratch/${{SLURM_JOBID}}/raw.bam 2>> {log} | \
        samtools sort -@ ${{SLURM_CPUS_ON_NODE}} -o {output.bam} -
        ########################
        #### Index bam file ####
        samtools index {output.bam}
        ############################
        #### Flagstat bam files ####
        samtools flagstat /lscratch/${{SLURM_JOBID}}/raw.bam > {output.flag_raw}
        samtools flagstat {output.bam} > {output.flag_mapq}
        """


rule rmDup:
        """
        This removes optical and PCR duplicates.
        """
        input:  "Bams/{sample}.bam"
        output:
                bam = "Bams/{sample}_rmdup.bam",
                met = "QCstats/{sample}.rmdup_metrix.txt",
                flag = "QCstats/{sample}_Flagstat_picard.txt"
        log:    "logs/rmDup.{sample}.log"
        version: config["picard"]
        params:
                rulename = "rmDup",
                batch   = config["job_rmDup"],
                samtools = config["samtools"],
        shell:  """ \
        module load picard/{version} || exit 1
        module load samtools/{params.samtools} || exit 1
        java -XX:ParallelGCThreads=${{SLURM_CPUS_ON_NODE}} \
        -Xmx${{SLURM_MEM_PER_NODE}}m \
        -jar ${{PICARDJARPATH}}/picard.jar MarkDuplicates \
        I={input} \
        OUTPUT={output.bam} \
        M={output.met} \
        REMOVE_DUPLICATES=true \
        CREATE_INDEX=true \
        TMP_DIR=/lscratch/${{SLURM_JOBID}}
        samtools flagstat {output.bam} > {output.flag}
        """


rule macs:
        """
        This finds peaks in the bam file.
        """
        input:  "Bams/{sample}_rmdup.bam"
        output: "Beds/{sample}_peaks.narrowPeak",
        log:    "logs/macs.{sample}.log"
        version: config["macs2"]
        params:
                rulename = "macs",
                batch = config["job_macs"],
                base = "{sample}"
        shell:  """
        module load macs/{version}
        macs2 callpeak -t {input} \
        -p 0.01 -f BAMPE \
        -n {params.base} \
        --outdir Beds \
        --keep-dup all \
        --tempdir /lscratch/${{SLURM_JOB_ID}}
        """


rule blklist:
        """
        This removes peaks in blacklisted regions.
        """
        input:
                bed = "Beds/{sample}_peaks.narrowPeak",
                blk = config["blacklist"]
        output: "Beds/{sample}_peaks_rmblklst.narrowPeak"
        log:    "logs/blacklist.{sample}.log"
        version: config["bedtools"]
        params:
                rulename = "blklist",
                batch = config["job_blklist"],
        shell:  """
        module load bedtools/{version}
        echo {wildcards.sample} 2>> {log}
        bedtools intersect -a {input.bed} -b {input.blk} -wa | wc -l 2>> {log}
        bedtools intersect -a {input.bed} -b {input.blk} -v > {output}
        """


rule bigWig:
        """
        This creates a raw and normalized bigWig files.
        """
        input:  "Bams/{sample}_rmdup.bam"
        output:
                raw = "Beds/{sample}_peaks_raw.bw",
                bpm = "Beds/{sample}_peaks_bpm.bw",
        log:    "logs/bigWig.{sample}.log"
        version: config["deeptools"]
        params:
                rulename = "bigWig",
                batch = config["job_bigWig"]
        shell:  """
        module load deeptools/{version}
        ###########################
        #### Make bigWig files ####
        bamCoverage -b {input} -o {output.raw} 2>> {log}
        bamCoverage -b {input} -o {output.bpm} \
        --binSize 20 \
        --normalizeUsing BPM \
        --centerReads \
        -p ${{SLURM_CPUS_ON_NODE}} 2>> {log}
        """


rule qcBigWigSum:
        """
        This generates multi-bigWig summary using deeptools.
        """
        input:  BIGWIG
        output: "QCstats/bigWig_summary.npz"
        log:    "logs/qcBigWigSum.log"
        version: config["deeptools"]
        params:
                rulename = "qcBigWigSum",
                batch = config["job_qcBigWigSum"]
        shell:  """
        module load deeptools/{version}
        #############################
        #### Make bigWig summary ####
        multiBigwigSummary bins \
        --bwfiles {input} \
        -o {output} \
        --smartLabels \
        --binSize 1000 \
        -p ${{SLURM_CPUS_ON_NODE}} 2>> {log}
        """

rule qcDeepPlots:
        """
        This generates QC plots using deeptools.
        """
        input:
                sum = "QCstats/bigWig_summary.npz",
                bam = RMBAMS
        output:
                pca = "QCstats/Summary_PCA.pdf",
                heat = "QCstats/Summary_CorrHeatmap.pdf",
                finger = "QCstats/Summary_Fingerprint.pdf",
                pca_tsv = "QCstats/Summary_PCA.tsv",
                cor_tsv = "QCstats/Summary_CorrHeamtmap.tsv",
                finger_tsv = "QCstats/Summary_Fingerprint.tsv"
        log:    "logs/qcDeepPlots.log"
        version: config["deeptools"]
        params:
                rulename = "qcDeepPlots",
                batch = config["job_qcDeepPlots"]
        shell:  """
        module load deeptools/{version}
        #########################
        #### bigWig QC plots ####
        # PCA #
        plotPCA --corData {input.sum} \
        -o {output.pca} \
        --transpose \
        --plotFileFormat pdf \
        --plotHeight 10 \
        --plotWidth 16 \
        --outFileNameData {output.pca_tsv} 2>> {log}
        # Heatmap #
        plotCorrelation --corData {input.sum} \
        -o {output.heat} \
        --corMethod pearson \
        --whatToPlot heatmap \
        --plotFileFormat pdf \
        --colorMap RdYlBu \
        --plotNumbers \
        --outFileCorMatrix {output.cor_tsv} 2>> {log}
        # Fingerprint #
        plotFingerprint \
        -b {input.bam} \
        -o {output.finger} \
        --plotFileFormat pdf \
        --smartLabels \
        --skipZeros \
        --numberOfSamples 50000 \
        -p ${{SLURM_CPUS_ON_NODE}} \
        --outRawCounts {output.finger_tsv} 2>> {log}
        """


# rule :
#         input:
#         output:
#         log:
#         version:
#         params:
#                 rulename =
#                 batch =
#         shell:  """
#     """
