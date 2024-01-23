from os.path import join


# have to use lambda wildcard functions here because cant pass both config["star_version"] and use
# the other wildcards in the same input line. Until I find a solution


rule feature_counts:
    input:
        sample=lambda wc: f"procdata/rnaseq/star_v{config['star_version']}_{wc.ref_build}.{wc.ref_version}/{wc.sample}/Aligned.sortedByCoord.out.bam",
        annotation=lambda wc: gencodeAnnotation(
            dirPath="metadata/references",
            ref_build=wc.ref_build,
            gencode_ver=wc.ref_version,
        ),
    output:
        multiext(
            "results/rnaseq/featureCounts_v{featureCounts_VERSION}_{ref_build}.{ref_version}/{sample}",
            ".featureCounts",
            ".featureCounts.summary",
            ".featureCounts.jcounts",
        ),
    threads: 8
    params:
        strand=0,  # 0: unstranded, 1: stranded, 2: reverse stranded
        extra="-O --fracOverlap 0.2 -J -p",
    log:
        "logs/rnaseq/featureCounts_v{featureCounts_VERSION}_{ref_build}.{ref_version}/{sample}.log",
    conda:
        lambda wc: f"../envs/featureCounts.v{wc.featureCounts_VERSION}.yaml"
    shell:
        """
        featureCounts \
            -T {threads} \
            -a {input.annotation} \
            -o {output[0]} \
            -s {params.strand} \
            {params.extra} \
            {input.sample} \
            &> {log}
        """
