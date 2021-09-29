
###Modified bits of this script to test out whether it is better to use the partition of ingroup (+silly accidental inclusion of one WA Melichrus sample) or
### to use the original partition and just remove the outgroup. Result: The whole dataset partition without removing the outgroup had the most SNPs and the best RepAvg. 
##But the ingroup partition(after removing WA specimen) has 82287 SNPs and repavg is 100% for 78.9% of loci. The whole dataset-outgroup produced 78474 SNPs with 76% of loci having a 100%repavg. 
##So basically there isn't alot of difference, but it looks like it is better to use the ingroup partition. 


#Set up script for PCOAs
##install.packages("BiocManager")
##BiocManager::install("SNPRelate")
##install.packages("ggrepel")
library(dartR)
library(StAMPP)
library(gplots) 
library(plotly)
library(hierfstat)
library(poppr)
library(adegenet)
library(ggrepel)

##################################################################
## Setup
##################################################################
## Specify main parameters
getwd()
# setwd("C:/Users/hkenned6/Documents/Melichrus/Molecular/DArTseq_Melichrus/DArTseq_Melichrus_ingroup")
overwrite <- TRUE
## Directories: choose an indir and an analysis name
indir <- "data" 
analysis_name <- "MELIngroup"
runfilt <- FALSE

## Calculate other parameters
dartname <- list.files(indir,"SNP_2.csv")
if(length(dartname)==0)dartname <- list.files(indir,"SNP_1.csv")
dartfile <- file.path(indir,dartname)
metadataname <- list.files(indir,"metadata.csv")
metadatafile <- file.path(indir,metadataname)
outdir <- file.path("output",analysis_name)
if(!dir.exists("output")) dir.create("output")
if(!dir.exists("temp")) dir.create("temp")
if(dir.exists(outdir)){
  if(overwrite){
    print("Overwriting existing directory")
  }else print("Directory already exists")
} else dir.create(outdir)

##################################################################
## Save gl file
##################################################################
## Read data file. Note that metadata are not read in here, as
## the order must be the same, but the dartfile is not in a
## sensible order.
library(dartR)
# mygl <- gl.read.dart(dartfile)
# save (mygl, file= "MELIngroup.rdata")
## Inspect gl
#mygl

##################################################################
## Open saved gl file 
################################################################
load("MELIngroup.rdata")
##To Inspect genlight object and list names
mygl
#indNames(mygl)

## Remove any accidental duplicates in dataset (HTK96c.1 and HTK96e.1)
mygl0 <- gl.drop.ind(mygl, ind.list=c("HTK96c.1", "HTK96e.1", "ARC870b"), recalc = TRUE, mono.rm = TRUE)
mygl
###################################################################
## Read metadata file and sort according to sample order in dartfile
###################################################################
metadata0 <- read.csv(metadatafile, stringsAsFactors = F)
metadata <- metadata0[match(mygl@ind.names,metadata0$id),]

## Add sample metadata to mygl
mygl@pop <- factor(metadata$pop)
##Check that poulations and names are all there
#mygl@pop
#mygl@ind.names
##Return or set population assignment for individuals
#pop(mygl)
##List the population labels##
#levels(pop(mygl))

##################################################################
## Explore data quality
##################################################################
## Check the number of individual samples
#nInd(mygl)
## Labels for individuals
#indNames(mygl)
##  Base frequencies across all sites, plus proportions of snps by type:
##  transitions (a<>g or c<>t) and transversions (all other snps)
#gl.report.bases(mygl)
##  Missing data: proportion of SNPs with X missing genotypes 
gl.report.callrate(mygl)
###  Missing data: proportion of individuals with X missing genotypes
gl.report.callrate(mygl0,method = "ind")
# ##  Distances between 'tags' - identify oversplitting of loci 
# gl.report.hamming(mygl) # requires too much memory
# gl.report.hamming(mygl[,1:10],rs = 5)
# ##  Average repeatability of each snp
mygl0
gl.report.RepAvg(mygl)
# ##  Secondary snps occur when there are multiple snps in a single tag.
# ##  They will be in linkage disequilibrium, so we will thin the data to 
# ##  one snp per tag. Some may still be linked, but not many. 
# gl.report.secondaries(mygl)
# ##  Monomorphic sites should be filtered out. DArT will have done this
# ##  before the data are provided to us, but if we subsample it, we 
# ##  need to do that again.
# gl.report.monomorphs(mygl)
# ##  Check LD - a few for illustration
# # gl.report.ld(mygl) # takes a while without parallelisation
# gl.report.ld(mygl[,1:10])
# ##  HWE testing needs a population specification
# mygl@pop=factor(rep(1,nInd(mygl)))
# gl.report.hwe(mygl) # produces a large dataframe
#gl.report.hwe(mygl[,1:10])
##  Look at distribution of observed heterozygosity, as it can indicate
##  paralogous loci or clonal reproduction. 
# hist(gl.Ho(mygl))

############################################################
## Filter data - using the same object to limit memory usage
############################################################
if(runfilt){
  
  # ## Filtering for PCA #######################################

  ## Average repeatability
  mygl1 <- gl.filter.RepAvg(mygl1,threshold = 0.95)
  ## Remove monomorphic loci
  mygl1 <- gl.filter.monomorphs(mygl1)
  ## Filter out large amounts of missing data by locus
  mygl1 <- gl.filter.callrate(mygl1,method = "loc",threshold = 0.97)
  ##  Thin to one snp per tag. Usually best to choose the one with the highest minor allele frequency.
  mygl1 <- gl.filter.secondaries(mygl1,method = "best")
  ## Filter out individuals with too much missing data
  mygl1 <- gl.filter.callrate(mygl1,method = "ind",threshold = 0.91,recalc = TRUE, mono.rm = TRUE)
  ##  Filter out loci with very low minor allele frequencies
  ##  exclude step for pop stats as this biases statistics like heterozygosities
  ##  and analyses based on allele frequency spectra.
  mygl1 <- gl.filter.maf(mygl1,threshold = 0.05)
  nLoc(mygl1);nInd(mygl1)
  nLoc(mygl);nInd(mygl)
  mygl1
  ## Recalculate all of the metrics for your filtered genlight
  mygl1 <- gl.recalc.metrics(mygl1)
  
  #############################################################
  ## Save filtered gl and associated metadata
  #############################################################
  
  ## Save the population IDs as a vector
  pop0 <- mygl1@pop
  #pop0
  ## Save filtered data as a compressed R data object.
  save(mygl1,file = file.path("temp",paste0(analysis_name,"_filtered1.rda")))
  #mygl1
}else{
  load(file.path("temp",paste0(analysis_name,"_filtered1.rda")))
}


##match metadata to filtered dataset
metadata1 <- metadata0[match(mygl1@ind.names,metadata0$id),]
metadata1
