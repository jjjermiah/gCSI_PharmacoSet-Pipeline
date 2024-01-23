from os.path import join

# when running on an HPC cluster, star index requires a lot of memory
# resources required:
# 1. 20 cores
# 2. 100GB RAM
# minimum time required: 2 hours
# recommended time: 4 hours


refDir = "metadata/references"
ref_build = ["GRCh38", "GRCh37"]
ref_version = "44"
kallisto_version = "0.46.1"


rule star_index:
    input:
        gtf=lambda wildcards: gencodeAnnotation(
            dirPath=wildcards.refDir,
            ref_build=wildcards.ref_build,
            gencode_ver=wildcards.ref_version,
            species="human",
        ),
        fasta=lambda wildcards: gencodeGenome(
            dirPath=wildcards.refDir,
            ref_build=wildcards.ref_build,
            gencode_ver=wildcards.ref_version,
            species="human",
        ),
    output:
        directory(
            "{refDir}/{species}/{ref_build}_v{ref_version}/STAR.v{STAR_version}_index"
        ),
    params:
        sjdbOverhang=100,
    log:
        "logs/rnaseq/star_index/{refDir}{species}/{ref_build}_v{ref_version}/STAR.v{STAR_version}_index.log",
    threads: 24
    conda:
        lambda wildcards: f"../envs/star.v{wildcards.STAR_version}.yaml"
    shell:
        """
        STAR \
            --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang {params.sjdbOverhang} \
            --genomeDir {output} \
            > {log} 2>&1
        """


rule star_align_paired:
    input:
        fq1="rawdata/rnaseq/{sample}_1_1.rnaseq.fastq.gz",
        fq2="rawdata/rnaseq/{sample}_1_2.rnaseq.fastq.gz",
        # metadata/references/human/GRCh38_v44/STAR.v2.7.11_index
        index=directory(
            join(
                refDir,
                config["ref_species"],
                "{ref_build}_v{ref_version}",
                "STAR.v{STAR_VERSION}_index",
            )
        ),
    output:
        # note there are other possible outputs, but I am not using them for now
        # algned_sam = "procdata/rnaseq/star_v{STAR_VERSION}_{ref_build}.{ref_version}/{sample}_aligned.sam",
        sorted_bam="procdata/rnaseq/star_v{STAR_VERSION}_{ref_build}.{ref_version}/{sample}/Aligned.sortedByCoord.out.bam",
        transcript_bam="procdata/rnaseq/star_v{STAR_VERSION}_{ref_build}.{ref_version}/{sample}/Aligned.toTranscriptome.out.bam",
        log="logs/rnaseq/star_v{STAR_VERSION}_{ref_build}.{ref_version}/{sample}_aligned.log",
        sj="procdata/rnaseq/star_v{STAR_VERSION}_{ref_build}.{ref_version}/{sample}/SJ.out.tab",
    log:
        "logs/rnaseq/star_v{STAR_VERSION}_{ref_build}.{ref_version}/{sample}_aligned.log",
    threads: 20
    # params:
    #     extra = "--alignMatesGapMax 500000 \
    #     --outSAMunmapped Within \
    #     --outSAMtype BAM SortedByCoordinate"
    # script:
    #     "../scripts/rnaseq/star_align_paired.py"
    conda:
        lambda wildcards: f"../envs/star.v{wildcards.STAR_VERSION}.yaml"
    shell:
        """
        STAR \
        --runThreadN {threads} \
        --genomeDir {input.index} \
        --readFilesIn {input.fq1} {input.fq2} \
        --readFilesCommand gunzip -c \
        --outFileNamePrefix $(dirname {output.sorted_bam})/ \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outFilterMultimapNmax 1 \
        --outFilterMismatchNmax 3 \
        --outFilterMismatchNoverLmax 0.3 \
        --quantMode TranscriptomeSAM \
        --sjdbOverhang 100 \
        --alignIntronMax 500000 \
        --alignMatesGapMax 500000 \
        > {log} 2>&1
        """
