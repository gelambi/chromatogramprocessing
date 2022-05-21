rm(list=ls())

### Deconvolution and alignment of GCMS amide samples using eRah###

### EFFECT OF BAT HANDLING ###
----#####

## DATA FORMAT: First, I used MSConvert to transform .D files to .mzXML. Parameters suggested in GNPS (https://ccms-ucsd.github.io/GNPSDocumentation/fileconversion/):

# Click Browse and select file(s) for conversion. Then click Add to add them to the MSConvert workflow.

# Choose an Output Directory

# Under Options, choose mzXML for output format, 32-bit for binary encoding precision and uncheck Use zlib compression.

# Under filters, choose Peak Picking with Vendor checked, in order to centroid the data. Indicate MS-Levels 1-2. Click Add to add the filter.

library(erah)
library(mzR) # it reads the mzXML files 

createdt("erahfolderbat") # This step creates the two .csv files 

ex_bat <- newExp(instrumental="erahfolderbat/erahfolderbat_inst.csv",
             phenotype = "erahfolderbat/erahfolderbat_pheno.csv") # create a new experiment

metaData(ex_bat) # check that everything looks good
phenoData(ex_bat)

# Deconvolution of GC-MS data

ex.dec.par <- setDecPar(min.peak.width = 1) # I am not sure about the 1, I tried 0.5 but the deconvolution failed

ex_bat_dec <- deconvolveComp(ex_bat, ex.dec.par)

# Save
save(ex_bat_dec, file = "deconvolution_batamides.rda")
# Load
load("deconvolution_batamides.rda")

# Alignment of GC-MS deconvolved compounds

ex.al.par <- setAlPar(min.spectra.cor = 0.90,
                      max.time.dist = 3,
                      mz.range = 70:600)

ex_align <- alignComp(ex_bat_dec, alParameters = ex.al.par)

ex_misscomp <- recMissComp(ex_align, min.samples = 3, free.model = FALSE) # This failed too

aligment <- alignList(ex_align, by.area=TRUE)
aligment

write.csv(aligment, file = "aligment_batdata.csv") # I got too many peaks

### EFFECT OF ANT HANDLING ###
----#####

createdt("erahfolderants") # This step creates the two .csv files 

ex_ants <- newExp(instrumental="erahfolderants/erahfolderants_inst.csv",
             phenotype = "erahfolderants/erahfolderants_pheno.csv") # create a new experiment

metaData(ex_ants) # check that everything looks good
phenoData(ex_ants)

# Deconvolution 

ex.dec.par <- setDecPar(min.peak.width = 1,
                        avoid.processing.mz = c(73:75,147:149))

ex <- deconvolveComp(ex_ants, ex.dec.par)

# Save
save(ex_ants, file = "testPCOS.rda")
# Load
load("testPCOS.rda")

# Alignment

ex.al.par <- setAlPar(min.spectra.cor = 0.90,
                      max.time.dist = 3,
                      mz.range = 70:600)

ex <- alignComp(ex_ants, alParameters = ex.al.par)
ex <- recMissComp(ex_ants, min.samples = 3)

aligment <- alignList(ex_ants, by.area=TRUE)
write.csv(aligment, file = "aligment_antdata.csv")