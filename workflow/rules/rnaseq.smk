
# THESE NEED TO BE DEFINED!
# 1. samples
# 2. config["ref_build"]
# 3. config["ref_version"]
# 4. config["kallisto_version"]
# 5. config["featureCounts_VERSION"]
# 6. config["rsem_version"]
# 7. config["star_version"]


include: "downloadGencode.smk"
include: "kallisto.smk"
include: "star.smk"
include: "featureCounts.smk"
include: "rsem.smk"

# rule kallisto:
#     input:
#         expand(
#             "results/rnaseq/kallisto_v{kallisto_version}_{ref_build}.{ref_version}/{sample}/",
#             sample=samples,
#             ref_build=config["ref_build"],
#             ref_version=config["ref_version"],
#             kallisto_version=config["kallisto_version"]
#             ),


# rule runFeatureCounts:
#     input:
#         expand(
#             "results/rnaseq/featureCounts_v{featureCounts_VERSION}_{ref_build}.{ref_version}/{sample}.featureCounts",
#             sample=samples,
#             ref_build=config["ref_build"],
#             ref_version=config["ref_version"],
#             featureCounts_VERSION=config["featureCounts_VERSION"]
#             ),

# rule calcRSEM:
#     input:
#         expand(
#             "results/rnaseq/rsem_v{RSEM_VERSION}_{ref_build}.{ref_version}/{sample}/a.genes.results",
#             sample = samples,
#             ref_build=config["ref_build"],
#             ref_version=config["ref_version"],
#             RSEM_VERSION = config["rsem_version"],
#         )
        