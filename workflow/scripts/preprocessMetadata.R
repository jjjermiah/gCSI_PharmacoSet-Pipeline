
## ------------------- Parse Snakemake Object ------------------- ##
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    save.image()
}

dt <- data.table::fread(INPUT[[1]])
# names(dt)
#  [1] "DrugName"            "CellLineName"        "PrimaryTissue"      
#  [4] "ExperimentNumber"    "TrtDuration"         "doublingtime"       
#  [7] "log10Concentration"  "GR"                  "relative_cell_count"
# [10] "activities.mean"     "activities.std"      "Norm_CellLineName"  
# [13] "Norm_DrugName"

dt <- unique(dt[, .(DrugName, CellLineName, PrimaryTissue, Norm_CellLineName, Norm_DrugName)])

# Do the following : numSamplesPerTreatment <- dt[, .N, by = .(Norm_DrugName, DrugName)][order(N, decreasing = TRUE)]
# but rename the N row as numSamplesPerTreatment
numSamplesPerTreatment <- dt[, .N, by = .(Norm_DrugName, DrugName)][order(N, decreasing = TRUE)]
names(numSamplesPerTreatment)[3] <- "numSamplesPerTreatment"

numTreatmentsPerSample <- dt[, .N, by = .(Norm_CellLineName, CellLineName, PrimaryTissue)][order(N, decreasing = TRUE)]
names(numTreatmentsPerSample)[4] <- "numTreatmentsPerSample"

data.table::fwrite(numSamplesPerTreatment, OUTPUT[["treatmentMetadata"]], sep = "\t", quote = FALSE, row.names = FALSE)
data.table::fwrite(numTreatmentsPerSample, OUTPUT[["sampleMetadata"]], sep = "\t", quote = FALSE, row.names = FALSE)

