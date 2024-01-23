

rule rsem_index:
    input:
        reference_genome=lambda wc: gencodeGenome(
            dirPath="metadata/references",
            ref_build=wc.ref_build,
            gencode_ver=wc.ref_version,
        ),
        gtf=lambda wc: gencodeAnnotation(
            dirPath="metadata/references",
            ref_build=wc.ref_build,
            gencode_ver=wc.ref_version,
        ),
    output:
        seq=multiext(
            "{refDir}/{species}/{ref_build}_v{ref_version}/RSEM.v{RSEM_VERSION}_index/reference.",
            "seq",
            "grp",
            "ti",
        ),
    threads: 8
    log:
        "logs/{refDir}/rsem_index/{species}/{ref_build}_v{ref_version}/RSEM.v{RSEM_VERSION}_index.log",
    conda:
        lambda wc: f"../envs/rsem.v{wc.RSEM_VERSION}.yaml"
    shell:
        """
        rsem-prepare-reference \
            --num-threads {threads} \
            --gtf {input.gtf} \
            {input.reference_genome} \
            $(dirname {output.seq[0]})/reference \
            > {log} 2>&1
        """


rule rsem_calculateExpression:
    input:
        bam=lambda wc: f"procdata/rnaseq/star_v{config['star_version']}_{wc.ref_build}.{wc.ref_version}/{wc.sample}/Aligned.toTranscriptome.out.bam",
        reference="metadata/references/human/{ref_build}_v{ref_version}/RSEM.v{RSEM_VERSION}_index/reference.seq",
    output:
        # this file contains per-gene quantification data for the sample
        genes_results="results/rnaseq/rsem_v{RSEM_VERSION}_{ref_build}.{ref_version}/{sample}/a.genes.results",
        # this file contains per-transcript quantification data for the sample
        isoforms_results="results/rnaseq/rsem_v{RSEM_VERSION}_{ref_build}.{ref_version}/{sample}/a.isoforms.results",
    params:
        extra="--seed 42 --paired-end",
    threads: 30
    conda:
        lambda wc: f"../envs/rsem.v{wc.RSEM_VERSION}.yaml"
    log:
        "logs/rnaseq/rsem_v{RSEM_VERSION}_{ref_build}.{ref_version}/{sample}/rsem_calculateExpression.log",
    shell:
        """
        rsem-calculate-expression \
            --num-threads {threads} \
            {params.extra} \
            --alignments {input.bam} \
            $(dirname {input.reference})/reference \
            {output.genes_results} \
            > {log} 2>&1
        """
