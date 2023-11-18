
## ------------------- Parse Snakemake Object ------------------- ##
if(exists("snakemake")){
    INPUT <- snakemake@input
    OUTPUT <- snakemake@output
    WILDCARDS <- snakemake@wildcards
    THREADS <- snakemake@threads
    save.image()
}

GRMetrics <- data.table::fread(INPUT[["GRMetrics"]])
# names(GRMetrics)
# names(GRMetrics)
#  [1] "DrugName"              "CellLineName"          "PrimaryTissue"        
#  [4] "ExperimentNumber"      "TrtDuration"           "doublingtime"         
#  [7] "maxlog10Concentration" "GR_AOC"                "meanviability"        
# [10] "GRmax"                 "Emax"                  "GRinf"                
# [13] "GEC50"                 "GR50"                  "R_square_GR"          
# [16] "pval_GR"               "flat_fit_GR"           "h_GR"                 
# [19] "N_conc"                "log10_conc_step"       "Norm_CellLineName"    
# [22] "Norm_DrugName"         "GR_05uM_fit"      

GRValues <- data.table::fread(INPUT[["GRValues"]])
# names(GRValues)
#  [1] "DrugName"            "CellLineName"        "PrimaryTissue"      
#  [4] "ExperimentNumber"    "TrtDuration"         "doublingtime"       
#  [7] "log10Concentration"  "GR"                  "relative_cell_count"
# [10] "activities.mean"     "activities.std"      "Norm_CellLineName"  
# [13] "Norm_DrugName"


raw_dt <- GRValues[, .(Norm_DrugName, Norm_CellLineName, ExperimentNumber, log10Concentration, GR, relative_cell_count, activities.mean, activities.std)]

# Steps from existing pset pipeline:
#' Create Experiment ID: 
#'      The script creates a new column "expid" which is a combination of 
#'      "cellid", "drugid", and "ExperimentNumber", separated by an underscore.

#' Calculate Maximum Concentration: 
#'      The script calculates the maximum concentration by counting the number of
#'      rows per "expid" and taking the maximum of these counts.

#' Create Unique IDs: 
#'      The script creates a vector of unique experiment IDs.

#' Initialize Matrices: 
#'      The script initializes two matrices, "doses_final" and "viability_final", 
#'      with dimensions equal to the number of unique experiment IDs and the 
#'      maximum concentration. The matrices are filled with NA values.

#' Order Data: 
#'      The script orders the data by "expid" and "log10Concentration" using the 
#'      setorder function.

#' Split Data: The script splits the data by "expid" using the split function, 
#'      resulting in a list of data.tables, each corresponding to a unique experiment.
