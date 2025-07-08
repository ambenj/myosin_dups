#!/usr/bin/env python
from os.path import join, abspath, expanduser

# Read in config file
OUTDIR = config['outdir']
REF_ALIGN = config['ref_align']
REF_SIM_PATH = config['ref_sim_path']
REF_SIM = config['ref_sim']
READ_LENGTH = config['read_length']
FRAG_MEAN = config['frag_mean']
FRAG_SD = config['frag_sd']
DEPTHS = config['depth']
REPS = list(range(int(config['reps'])))
REGIONS = config['regions']

# Convert out directory to absolute path
if OUTDIR[0] == '~':
    OUTDIR = expanduser(OUTDIR)
OUTDIR = abspath(OUTDIR)

# Get ref prefix
REF_PREFIX = os.path.splitext(os.path.basename(REF_ALIGN))[0]

#########################################################################

rule all:
    input:
#        expand(join(OUTDIR, "01_alignment/{ref_prefix}/{prefix}_sim{depth}X.rep{rep}.RG.sorted.bam"), ref_prefix=REF_PREFIX, prefix=REF_SIM, depth=DEPTHS, rep=REPS),
        expand(join(OUTDIR, "02_coverage/{ref_prefix}/{prefix}_sim{depth}X.rep{rep}.mapq3_roi_coverage.txt"), ref_prefix=REF_PREFIX, prefix=REF_SIM, depth=DEPTHS, rep=REPS),
        expand(join(OUTDIR, "02_coverage/{ref_prefix}/{prefix}_sim{depth}X.rep{rep}.mapq3_wg_coverage.txt"), ref_prefix=REF_PREFIX, prefix=REF_SIM, depth=DEPTHS, rep=REPS),


#########################################################################
########################### Simulate genomes #########################
#########################################################################


rule simulate:
    input:
        source_fasta = join(REF_SIM_PATH, "{prefix}"),
    output:
        fq1 = join(OUTDIR, "00_simulations/{prefix}_sim{depth}X.rep{rep}_1.fq.gz"),
        fq2 = join(OUTDIR, "00_simulations/{prefix}_sim{depth}X.rep{rep}_2.fq.gz")
    params:
        prefix = join(OUTDIR, "00_simulations/{prefix}_sim{depth}X.rep{rep}_"),
        fq1 = join(OUTDIR, "00_simulations/{prefix}_sim{depth}X.rep{rep}_1.fq"),
        fq2 = join(OUTDIR, "00_simulations/{prefix}_sim{depth}X.rep{rep}_2.fq"),
    resources:
        time = 3,
        mem = 32
    shell:"""
        # Load modules
        module load art/20160605

        # Set seed
        seed=$RANDOM
        echo "Random seed for illumina art is $seed"

        # Simulate reads
        art_illumina -ss HS20 -na -i {input.source_fasta} -p -m {FRAG_MEAN} -s {FRAG_SD} -l {READ_LENGTH} -f {wildcards.depth} -o {params.prefix} -rs $seed

        # gzip fastq
        gzip {params.fq1}
        gzip {params.fq2}
    """

rule align:
    input:
        ref = REF_ALIGN,
        fq1 = rules.simulate.output.fq1,
        fq2 = rules.simulate.output.fq2
    output:
        bam = temp(join(OUTDIR, "01_alignment/{ref_prefix}/{prefix}_sim{depth}X.rep{rep}.sorted.bam")),
        bam_rg = join(OUTDIR, "01_alignment/{ref_prefix}/{prefix}_sim{depth}X.rep{rep}.RG.sorted.bam")
    resources:
        time = 12,
        mem = 32
    shell: """
        # Load modules
        module load bwa/0.7.17
        module load samtools/1.19
        module load picard/3.1.1

        # Align reads
        bwa mem {input.ref} {input.fq1} {input.fq2} | samtools view -b - | samtools sort - > {output.bam}
        # Add read groups
        picard AddOrReplaceReadGroups -I {output.bam} -O {output.bam_rg} -SO coordinate -ID 1 -RGLB lib1 -RGPL illumina -RGPU unit1 -RGSM "{wildcards.prefix}.{wildcards.depth}X.rep{wildcards.rep}"
        # Index bam
        samtools index {output.bam_rg}
    """

rule roi_cov:
    input:
        bam = rules.align.output.bam_rg,
        regions = REGIONS,
    output:
        roi = join(OUTDIR, "02_coverage/{ref_prefix}/{prefix}_sim{depth}X.rep{rep}.mapq3_roi_coverage.txt")
    resources:
        time = 3,
        mem = 8
    shell:"""
        # Load modules
        module load samtools/1.19

        BAM={input.bam}
        REGIONS={input.regions}

        printf "desc\tregion\tchr\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\tbam\n" > {output.roi}

        # Get coverage for rois
        while IFS=$'\t' read -r -a region_list
        do
                region="${{region_list[0]}}:${{region_list[1]}}-${{region_list[2]}}"
                result=$(samtools coverage -H -q 3 --ff 260 -r $region $BAM)
                printf "${{region_list[3]}}\t${{region}}\t${{result}}\t${{BAM}}\n" >> {output.roi}

        done < $REGIONS
    """

rule wg_cov:
    input:
        bam = rules.align.output.bam_rg,
    output:
        wg = join(OUTDIR, "02_coverage/{ref_prefix}/{prefix}_sim{depth}X.rep{rep}.mapq3_wg_coverage.txt")
    resources:
        time = 3,
        mem = 16
    shell:"""
        # Load modules
        module load samtools/1.19

        printf "chr\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n" > {output.wg}

        # Get coverage for whole genome
        samtools coverage -H -q 3 --ff 260 {input.bam} >> {output.wg}

    """
