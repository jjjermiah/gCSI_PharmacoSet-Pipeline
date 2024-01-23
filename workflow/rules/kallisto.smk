refDir = "references"
ref_build = ["GRCh38", "GRCh37"]
ref_version = ["44", "45"]
kallisto_version = ["0.46.1", "0.46.2", "0.48.0"]


rule build_all_kallisto_indices:
    input:
        expand(
            "{ref_dir}/human/{ref_build}_v{ref_version}/transcriptome-kallisto_v{kallisto_version}.idx", 
            ref_dir=refDir,
            ref_build=ref_build,
            ref_version=ref_version,
            kallisto_version=kallisto_version)

rule kallisto_index:
    input:
        ref = lambda wildcards: 
            gencodeTranscriptome(
                dirPath = wildcards.refDir, 
                ref_build = wildcards.ref_build, 
                gencode_ver = wildcards.ref_version, 
                species = wildcards.species,)
    output:
        index = "{refDir}/{species}/{ref_build}_v{ref_version}/transcriptome-kallisto_v{kallisto_version}.idx"
    log:
        "logs/rnaseq/kallisto/{refDir}/{species}_{ref_build}.{ref_version}/kallisto-index_v{kallisto_version}.log"
    conda:
        lambda wildcards: f"../envs/kallisto.v{wildcards.kallisto_version}.yaml"
    envmodules:
        lambda wildcards: f"kallisto/{wildcards.kallisto_version}"
    threads:
        2
    shell:
        """
        # OPTIONS=""
        # if [ "{wildcards.kallisto_version}" == "0.50.1" ]; then
        #     OPTIONS="--threads {threads}"
        # fi

        kallisto index \
            -i {output.index} \
            {input.ref} # > {log} 2>&1
        """

rule kallisto_quant:
    input:
        fastq =
            lambda wc: 
                ["rawdata/rnaseq/{}_1_{}.rnaseq.fastq.gz".format(wc.sample, pair) for pair in [1,2]],
        index =
            lambda wc: 
                "metadata/references/human/{}_v{}/transcriptome-kallisto_v{}.idx".format(
                    wc.ref_build, wc.ref_version, wc.kallisto_version)
    output:
        output_dir = 
            directory("results/rnaseq/kallisto_v{kallisto_version}_{ref_build}.{ref_version}/{sample}/")
    threads:
        10
    log:
        "logs/rnaseq/kallisto_v{kallisto_version}/{sample}_{ref_build}.{ref_version}_kallisto-quant.log"
    conda:
        lambda wildcards: f"../envs/kallisto.v{wildcards.kallisto_version}.yaml"
    envmodules:
        lambda wildcards: f"kallisto/{wildcards.kallisto_version}"    
    shell:
        """
        kallisto quant \
            --threads {threads} \
            --index {input.index} \
            --output-dir {output.output_dir} \
            {input.fastq[0]} {input.fastq[1]} \
            > {log} 2>&1
        """