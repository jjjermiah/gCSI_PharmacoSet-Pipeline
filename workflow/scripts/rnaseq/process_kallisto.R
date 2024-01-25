library(tximport)
library(GenomicRanges)

annotation.gtf <- "/home/bioinf/bhklab/jermiah/psets/PharmacoSet-Pipelines/gCSI/references/human/GRCh38_v45/annotation.gtf"
KALLISTO_DIR <- "/home/bioinf/bhklab/jermiah/psets/PharmacoSet-Pipelines/gCSI/procdata/rnaseq/kallisto_v0.46.1_GRCh38.45"

txdb <- GenomicFeatures::makeTxDbFromGFF(
    file = annotation.gtf,
    dataSource = "GENCODE GRCh38_v45"
)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")


samples <- list.files(KALLISTO_DIR)

abundance_files <- list.files(file.path(path, samples), full.names=T, pattern = "*abundance.tsv")
names(abundance_files) <- samples

rnaseq.genes <- tximport(h5_files, type = "kallisto", tx2gene = tx2gene, ignoreAfterBar = TRUE, ignoreTxVersion = FALSE)

# countsFromAbundance = c("no", "scaledTPM", "lengthScaledTPM", "dtuScaledTPM"),
countsFromAbundance <- "no"

rnaseq.transcripts <- tximport(h5_files, type = "kallisto", txOut = T, countsFromAbundance = countsFromAbundance)

rownames(rnaseq.transcripts$counts) <- sub("\\|.*", "", rownames(rnaseq.transcripts$counts))
rownames(rnaseq.transcripts$abundance) <- sub("\\|.*", "", rownames(rnaseq.transcripts$abundance))


assays <- list(
    kallisto.genes = rnaseq.genes$abundance,
    kallisto.genes_counts = rnaseq.genes$counts,
    kallisto.transcripts = rnaseq.transcripts$abundance,
    kallisto.transcripts_counts = rnaseq.transcripts$counts
)

