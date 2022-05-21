

library(erah)
library(mzR) # it reads the mzXML files 

createdt("erahalkenylphenols") # This step creates the two .csv files 

ex_alk <- newExp(instrumental="erahalkenylphenols/erahalkenylphenols_inst.csv",
                 phenotype = "erahalkenylphenols/erahalkenylphenols_pheno.csv") # create a new experiment

metaData(ex_alk) # check that everything looks good
phenoData(ex_alk)

# Deconvolution 

ex.dec.par <- setDecPar(min.peak.width = 1)

ex_bat_dec <- deconvolveComp(ex_alk, ex.dec.par)

# Save
save(ex_bat_dec, file = "deconvolution_alk.rda")
# Load
load("deconvolution_alk.rda")

# Alignment

ex.al.par <- setAlPar(min.spectra.cor = 0.90,
                      max.time.dist = 3,
                      mz.range = 70:600)

ex_alk <- alignComp(ex_alk, alParameters = ex.al.par)

ex_alk <- recMissComp(ex_alk, min.samples = 3, free.model = FALSE)

alignment <- alignList(ex_alk, by.area=TRUE)
alignment

write.csv(alignment, file = "alignment_alkenylphenols.csv")