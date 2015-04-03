library(NanoStringQCPro)

dataDir <- system.file("extdata", package="NanoStringQCPro")

example_rccSet <- newRccSet(
     rccDir                 = file.path(dataDir, "RCC")
    ,rlfPath                = file.path(dataDir, "RLF", "NQCP_example.rlf")
    ,cdrDesignData   = file.path(dataDir, "CDR", "CDR-DesignData.csv")
    ,extraPdata        = file.path(dataDir, "extraPdata", "SampleType.txt")
    ,blankLabel             = "blank"
    ,addEgAnnotations       = TRUE
    ,experimentData.name    = "Dorothee Nickles"
    ,experimentData.lab     = "Richard Bourgon"
    ,experimentData.contact = "nickles.dorothee@gene.com"
    ,experimentData.title   = "NanoStringQCPro example dataset"
    ,experimentData.abstract= "Example data for the NanoStringQCPro package"
)

save(list="example_rccSet", file="../../data/example_rccSet.rdata")

norm_example_rccSet <- preprocRccSet(
     rccSet         = example_rccSet
    ,method    = "median"
)

qc_example_rccSet <- makeQCReport(
     rccSet             = norm_example_rccSet
    ,outputBaseName     = "example_QC_report"
)
