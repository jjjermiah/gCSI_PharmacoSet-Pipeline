refDir = "references"
ref_build = ["GRCh38", "GRCh37"]
ref_version = [str(i) for i in range(33, 45)]  # ["44", "45"]
kallisto_version = ["0.46.1", "0.46.2", "0.48.0"]


rule quant_sample_586986:
    input:
        h5=expand(
            "results/rnaseq/kallisto_v{kallisto_version}_{ref_build}.{ref_version}/{sample}/abundance.h5",
            kallisto_version=kallisto_version,
            ref_build=ref_build,
            ref_version=ref_version,
            sample="586986",
        ),
        tsv=expand(
            "results/rnaseq/kallisto_v{kallisto_version}_{ref_build}.{ref_version}/{sample}/abundance.tsv",
            kallisto_version=kallisto_version,
            ref_build=ref_build,
            ref_version=ref_version,
            sample="586986",
        ),


rule build_all_kallisto_indices:
    input:
        expand(
            "{ref_dir}/human/{ref_build}_v{ref_version}/transcriptome-kallisto_v{kallisto_version}.idx",
            ref_dir=refDir,
            ref_build=ref_build,
            ref_version=ref_version,
            kallisto_version=kallisto_version,
        ),


rule kallisto_index:
    input:
        ref=lambda wildcards: gencodeTranscriptome(
            dirPath=wildcards.refDir,
            ref_build=wildcards.ref_build,
            gencode_ver=wildcards.ref_version,
            species=wildcards.species,
        ),
    output:
        index="{refDir}/{species}/{ref_build}_v{ref_version}/transcriptome-kallisto_v{kallisto_version}.idx",
    log:
        "logs/rnaseq/kallisto/{refDir}/{species}_{ref_build}.{ref_version}/kallisto-index_v{kallisto_version}.log",
    conda:
        lambda wildcards: f"../envs/kallisto.v{wildcards.kallisto_version}.yaml"
    envmodules:
        lambda wildcards: f"kallisto/{wildcards.kallisto_version}",
    resources:
        mem_mb=12000,
    threads: 2
    retries: 3
    shell:
        """
        kallisto index \
            -i {output.index} \
            {input.ref} 2>&1 | tee {log}
        """


rule kallisto_quant:
    input:
        fastq=lambda wc: [
            "rawdata/rnaseq/fastq/{}_1_{}.rnaseq.fastq.gz".format(wc.sample, pair)
            for pair in [1, 2]
        ],
        index=lambda wc: "references/human/{}_v{}/transcriptome-kallisto_v{}.idx".format(
            wc.ref_build, wc.ref_version, wc.kallisto_version
        ),
    output:
        files=multiext(
            "results/rnaseq/kallisto_v{kallisto_version}_{ref_build}.{ref_version}/{sample}/",
            "abundance.h5",
            "abundance.tsv",
            "run_info.json",
        ),
    threads: 8
    log:
        "logs/rnaseq/kallisto_v{kallisto_version}/{sample}_{ref_build}.{ref_version}_kallisto-quant.log",
    conda:
        lambda wildcards: f"../envs/kallisto.v{wildcards.kallisto_version}.yaml"
    envmodules:
        lambda wildcards: f"kallisto/{wildcards.kallisto_version}",
    shell:
        """
        kallisto quant \
            --threads {threads} \
            --index {input.index} \
            --output-dir $(dirname {output.h5}) \
            {input.fastq[0]} {input.fastq[1]} \
            2>&1 | tee {log}
        """
