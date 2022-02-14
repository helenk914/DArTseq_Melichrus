### Looped script, to generate subsets of mygl0, save subsets.
###Filter subsets, save filtered subsets.

##Create subsets of mygl0 and save them.EDIT for different subsetting levels.
for (subsets in c("AAA","BBA","BAA", "AAB")){
  Subsetgl <- mygl0[match(metadata_ingroup$id[which(metadata_ingroup$Tertiary_PCA_group==subsets)],mygl0$ind.names),]
  Subsetname = paste0("mygl0", subsets)
  save(Subsetgl,file = file.path("temp",paste0(Subsetname, ".rda")))
  ## Average repeatability
  Subsetgl1 <- gl.filter.RepAvg(Subsetgl,threshold = 0.95)
  ## Filter out large amounts of missing data by locus
  Subsetgl1 <- gl.filter.callrate(Subsetgl1,method = "loc",threshold = 0.95)
  ##  Thin to one snp per tag. Usually best to choose the one with the highest minor allele frequency.
  Subsetgl1 <- gl.filter.secondaries(Subsetgl1,method = "best")
  ## Filter out individuals with too much missing data
  Subsetgl1 <- gl.filter.callrate(Subsetgl1,method = "ind",threshold = 0.9)
  ##  Filter out loci with very low minor allele frequencies
  ##  exclude step for pop stats as this biases statistics like heterozygosities
  ##  and analyses based on allele frequency spectra.
  Subsetgl1 <- gl.filter.maf(Subsetgl1,threshold = 0.05)
  nLoc(Subsetgl1);nInd(Subsetgl1)
  # nLoc(mygl0);nInd(mygl0)
  ## Recalculate all of the metrics for your filtered genlight
  Subsetgl1 <- gl.recalc.metrics(Subsetgl1)
  ##Save filtered subset
  Subsetname1 = paste0("SubFiltgl", subsets)
  save(Subsetgl1,file = file.path("temp",paste0(Subsetname1, ".rda")))
}

################################################################################
## Check output genlights have enough SNPs, low missing data, and haven't lost
## many individual samples

# load(file = file.path("temp",paste0("mygl0C.rda")))
# load(file = file.path("temp",paste0("SubFiltglC.rda")))
# Subsetgl
# Subsetgl1
