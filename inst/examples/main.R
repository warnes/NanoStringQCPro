
library(NanoStringQCPro)

dataDir <- system.file("extdata", package="NanoStringQCPro")
rccDir <- file.path(dataDir, "RCC")

example_rccSet <- newRccSet(
     rccFiles               = file.path(rccDir, dir(rccDir))
    ,rlf                    = file.path(dataDir, "RLF", "NQCP_example.rlf")
    ,cdrDesignData          = file.path(dataDir, "CDR", "CDR-DesignData.csv")
    ,extraPdata             = file.path(dataDir, "extraPdata", "SampleType.txt")
    ,blankLabel             = "blank"
    ,addEgAnnotations       = TRUE
    ,experimentData.name    = "Dorothee Nickles"
    ,experimentData.lab     = "Richard Bourgon"
    ,experimentData.contact = "nickles.dorothee@gene.com"
    ,experimentData.title   = "NanoStringQCPro example dataset"
    ,experimentData.abstract= "Example data for the NanoStringQCPro package"
)

save(list="example_rccSet", file="../../data/example_rccSet.rdata")

norm_example_rccSet <- preprocRccSet(example_rccSet)

qc_example_rccSet <- makeQCReport(norm_example_rccSet,
                                  outputBaseName="example_QC_report",
                                  heatmaps=TRUE,
                                  verbose=TRUE)

