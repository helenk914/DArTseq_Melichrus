################################################################################
## Run 'Set up script for PCOA' first to set up and load mygl (don't need to filter)
################################################################################
## Looped script to;
## Generate subsets of mygl based on columns in metadata, save subsets
## Filter subsets, save filtered subsets.

##Pre-subsetting filtering of mygl
## Remove accidental duplicates in dataset (HTK96c.1, HTK96e.1 and accidental outgroup ARC870b)
mygl <- gl.drop.ind(mygl, ind.list=c("HTK96c.1", "HTK96e.1","ARC870b"), recalc = TRUE, mono.rm = TRUE)

metadata <- metadata0[match(mygl@ind.names,metadata0$id),]

##Create subsets of mygl and save them.EDIT for different subsetting levels.
for (subsets in c("AAAA")){
  Subsetgl <- mygl[match(metadata$id[which(metadata$Quarternary_PCA_group==subsets)],mygl$ind.names),]
  Subsetname = paste0("mygl", subsets)
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
  # nLoc(mygl);nInd(mygl)
  ## Recalculate all of the metrics for your filtered genlight
  Subsetgl1 <- gl.recalc.metrics(Subsetgl1)
  ##Save filtered subset
  Subsetname1 = paste0("SubFiltgl", subsets)
  save(Subsetgl1,file = file.path("temp",paste0(Subsetname1, ".rda")))
}

################################################################################
## Check output genlights have enough SNPs, low missing data, and haven't lost
## many individual samples

# load(file = file.path("temp",paste0("myglAAAA.rda")))
# load(file = file.path("temp",paste0("SubFiltglAAAA.rda")))
# Subsetgl
# Subsetgl1

