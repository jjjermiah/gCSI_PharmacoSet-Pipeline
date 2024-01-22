from os.path import join
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

rule downloadTreatmentResponseANDMetadata:
    input: 
        metadata = HTTP.remote(config["treatmentResponse_and_Metadata"]["TSV_files"])
    output:
        GRMetrics = "rawdata/TreatmentResponse_and_Metadata/gCSI_GRmetrics_v1.3.tsv",
        GRValues = "rawdata/TreatmentResponse_and_Metadata/gCSI_GRvalues_v1.3.tsv",
    shell: 
        """
        tar -xzf {input.metadata} -C rawdata/TreatmentResponse_and_Metadata;
        find rawdata/TreatmentResponse_and_Metadata -name '*.tsv' -exec mv '{{}}' rawdata/TreatmentResponse_and_Metadata/;
        """

rule downloadExpressionMetadata:
    input:
        HTTP.remote(config["molecularProfiles"]["expression"]["metadata"]["url"])
    output:
        join(metadata,"EGAD000010000725/EGAD000010000725_map.txt")
    shell:
        "wget -O {output} {input}"


