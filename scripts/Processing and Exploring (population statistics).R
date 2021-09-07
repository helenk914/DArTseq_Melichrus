##################################################################
##################################################################
##R.Andrew Original script
##Modified by H.Kennedy 2021 
##Processing, exploration of DArTseq data by population statistics
##################################################################
##################################################################
## Requires *SNP_2.csv file supplied by DArT and matching sample
## *metadata.csv in a folder specific to the DArT run (within 
## the data directory). 
##################################################################
##Products are; mygl, unfiltered genlight object (load), mygl2, filtered gl
## for population statistics. 
##################################################################
## Install and Load packages
##################################################################
##Check where packages are installed 
## .libPaths()
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
##setwd("C:/Users/hkenned6/Documents/Melichrus/Molecular/DArTseq_Melichrus/")
overwrite <- TRUE
## Directories: choose an indir and an analysis name
indir <- "data" 
analysis_name <- "MELIngroup"

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
## Load data
##################################################################
## Read data file. Note that metadata are not read in here, as
##  the order must be the same, but the dartfile is not in a
##  sensible order.
##library(dartR)
##mygl <- gl.read.dart(dartfile)
##save (mygl, file= "MEL.rdata")
##Inspect gl
##mygl

##################################################################
## Open saved gl file 
##################################################################
load("MELIngroup.rdata")
##To Inspect genlight object and list names
mygl
#indNames(mygl)

###################################################################
## Read metadata file and sort according to sample order in dartfile
###################################################################
metadata0 <- read.csv(metadatafile, stringsAsFactors = F)
metadata <- metadata0[match(mygl@ind.names,metadata0$id),]

## Add sample metadata to mygl
mygl@pop <- factor(metadata$pop)
##Check that poulations and names are all there
mygl@pop
mygl@ind.names
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
#gl.report.callrate(mygl)
###  Missing data: proportion of individuals with X missing genotypes
# gl.report.callrate(mygl,method = "ind")
# ##  Distances between 'tags' - identify oversplitting of loci 
# gl.report.hamming(mygl) # requires too much memory
# gl.report.hamming(mygl[,1:10],rs = 5)
# ##  Average repeatability of each snp
# gl.report.RepAvg(mygl)
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
## hist(gl.Ho(mygl))
##  If the species reproduces sexually, it is unlikely that observed 
##  heterozygosity will be > 0.5 - how unlikely will depend on the sample 
##  size, so choose a threshold based on your gl.report.heterozygosity.
## gl.report.heterozygosity(mygl)
## hist(gl.Ho(mygl))

##############################################################
## Subset data to include only populations with 5 or more individuals
##############################################################

## Examine spread in number of individuals per population (sample size can effect heterozygosity)
as.data.frame(table(pop(mygl))) %>% arrange(Freq)
barplot(table(pop(mygl)), las=2)
## Use, if barplot doesn't show all populations. 
#library(tidyverse)
NumberInd %>% arrange(Freq)
names(NumberInd)[1] = 'pop'

## Subset mygl to include only pops with 5 or more individuals ##

myglSP <- mygl[mygl@pop %in% levels(pop(mygl))[which(table(pop(mygl)) > 4)],]

############################################################
## Filter data - using the same object to limit memory usage
############################################################
## Filtering for population stats ##########################
## Remove any accidental duplicates in dataset (HTK96c.1 and HTK96e.1)
mygl2<- gl.drop.ind(mygl, ind.list=c("HTK96c.1", "HTK96e.1"))
## Average repeatability
mygl2 <- gl.filter.RepAvg(mygl2,threshold = 0.95)
## Filter out large amounts of missing data by locus
mygl2 <- gl.filter.callrate(mygl2,method = "loc",threshold = 0.95)
##  Thin to one snp per tag. Usually best to choose the one with the highest minor allele frequency. 
mygl2 <- gl.filter.secondaries(mygl2,method = "best")
## Filter out individuals with too much missing data
mygl2 <- gl.filter.callrate(mygl2,method = "ind",threshold = 0.87)
nLoc(mygl);nInd(mygl)
nLoc(mygl2);nInd(mygl2)
indNames(mygl2)
## Extra filtering ##########################################
## May be too intensive for PC, or more useful for PCA etc.

##  Filtering out over-split loci (likely already done by DArT) (need to use UNE computation)
##mygl1 <- gl.filter.hamming(mygl1,threshold = 0.2,rs = 5)
##  Filter out loci with very low minor allele frequencies
##  exclude step for pop stats as this biases statistics like heterozygosities
##  and analyses based on allele frequency spectra.
# mygl2 <- gl.filter.maf(mygl1,threshold = 0.05)


## Recalculate all of the metrics for your filtered genlight 
mygl2 <- gl.recalc.metrics(mygl2)

## It's often useful to save the population IDs as a vector
pop0 <- mygl2@pop

## To save time filtering in the future, save your filtered data as a
## compressed R data object.
save(mygl2,file = file.path("temp",paste0(analysis_name,"_filtered2.rda")))
#mygl2

##match metadata to filtered dataset
metadata2 <- metadata0[match(mygl2@ind.names,metadata0$id),]
metadata2

##################################################################
## Population-level statistics
##################################################################

## Compare heterozygosity among populations
ho_pop <- gl.report.heterozygosity(mygl2)
pdf(file.path(outdir,paste0(analysis_name,"_ho_pop.pdf")));gl.report.heterozygosity(mygl2);dev.off()
write.table(ho_pop,file = file.path(outdir,paste0(analysis_name,"_ho_pop.txt")))
#gl.report.maf(mygl) # takes a lot of time

## Examine private alleles by pairs of populations
pa_pw_pop <- gl.report.pa(mygl2)
write.table(pa_pw_pop,file = file.path(outdir,paste0(analysis_name,"_pa_pw_pop.txt")))

## Genetic distances calculated in dartR
mydist1 <- gl.dist.pop(mygl2,method = "pcfixed"); mydist1
mydist2 <- gl.dist.pop(mygl2,method = "pa"); mydist2
mydist3 <- gl.dist.pop(mygl2,method = "gower",upper = T); mydist3 # Why the warning??
mydist4 <- gl.dist.pop(mygl2,method = "euclidean",upper = T); mydist4
write.table(mydist1,file = file.path(outdir,paste0(analysis_name,"_pop_dist_pcfixed.txt")),sep="\t")
write.table(mydist2,file = file.path(outdir,paste0(analysis_name,"_pop_dist_pa.txt")),sep="\t")
write.table(as.matrix(mydist3),file = file.path(outdir,paste0(analysis_name,"_pop_dist_gower.txt")),sep="\t")
write.table(as.matrix(mydist4),file = file.path(outdir,paste0(analysis_name,"_pop_dist_euclidean.txt")),sep="\t")

## Nei's (1972) genetic distance
mygeno <- stamppConvert(mygl2,type = "genlight")
mydist5 <- stamppNeisD(mygeno,pop=T)
write.table(mydist5,file = file.path(outdir,paste0(analysis_name,"_pop_dist_Nei1972.txt")),sep="\t")

## F statistics with STAMPP
pwfst <-stamppFst(mygl2, nboots=1000, percent=95, nclusters=1)
pwfst$Fsts
pwfst$Pvalues
write.table(pwfst$Fsts,file = file.path(outdir,paste0(analysis_name,"_pop_fst.txt")),sep="\t")

## Heterozygosity and F statistics with hierfstat
sumstats <- basic.stats(cbind(mygl2@pop,as.data.frame(mygl2)))
sumstats$overall
pop_summary <- data.frame(N=colMeans(sumstats$n.ind.samp,na.rm=T),
                          Ho=colMeans(sumstats$Ho,na.rm=T),
                          He=colMeans(sumstats$Hs,na.rm=T),
                          Fis=colMeans(sumstats$Fis,na.rm=T))
write.table(pop_summary,file = file.path(outdir,paste0(analysis_name,"_pop_summary.txt")),sep="\t")
pop_summary

## Take a look at how Ho and He are distributed
hist(sumstats$perloc$Ho)
plot(sumstats$perloc$Hs,sumstats$perloc$Ho)
hist(sumstats$perloc$Fis)

##################################################################
## Population-level statistics
##################################################################
load(file = file.path("temp",paste0(analysis_name,"_filtered2.rda")))

# Compare heterozygosity among populations
mygl2Ho <- gl.report.heterozygosity(mygl2)
basic.stats(mygl2)
hist(gl.Ho(mygl2))
ho_pop <- gl.report.heterozygosity(mygl2)
# pdf(file.path(outdir,paste0(analysis_name,"_ho_pop.pdf")));gl.report.heterozygosity(mygl2);dev.off()
# write.table(ho_pop,file = file.path(outdir,paste0(analysis_name,"_ho_pop.txt")))
# # gl.report.maf(mygl) # takes a lot of time
###M.sp.mareeba have much higher heterozygosity than other populations, investigate Fis. 
##Basic.stats 
# 
# ## Examine private alleles by pairs of populations
# pa_pw_pop <- gl.report.pa.pop(mygl2)
# pa_pw_pop
# write.table(pa_pw_pop,file = file.path(outdir,paste0(analysis_name,"_pa_pw_pop.txt")))
# 
# ## Genetic distances calculated in dartR
# mydist1 <- gl.dist.pop(mygl2,method = "pcfixed"); mydist1
# mydist2 <- gl.dist.pop(mygl2,method = "pa"); mydist2
# mydist3 <- gl.dist.pop(mygl2,method = "gower",upper = T); mydist3 # Why the warning??
# mydist4 <- gl.dist.pop(mygl2,method = "euclidean",upper = T); mydist4
# write.table(mydist1,file = file.path(outdir,paste0(analysis_name,"_pop_dist_pcfixed.txt")),sep="\t")
# write.table(mydist2,file = file.path(outdir,paste0(analysis_name,"_pop_dist_pa.txt")),sep="\t")
# write.table(as.matrix(mydist3),file = file.path(outdir,paste0(analysis_name,"_pop_dist_gower.txt")),sep="\t")
# write.table(as.matrix(mydist4),file = file.path(outdir,paste0(analysis_name,"_pop_dist_euclidean.txt")),sep="\t")
# 
# ## Nei's (1972) genetic distance
# mygeno <- stamppConvert(mygl1,type = "genlight")
# mydist5 <- stamppNeisD(mygeno,pop=T)
# write.table(mydist5,file = file.path(outdir,paste0(analysis_name,"_pop_dist_Nei1972.txt")),sep="\t")
# 
# ## F statistics with STAMPP
# pwfst <-stamppFst(mygl2, nboots=1000, percent=95, nclusters=1)
# pwfst$Fsts
# pwfst$Pvalues
# write.table(pwfst$Fsts,file = file.path(outdir,paste0(analysis_name,"_pop_fst.txt")),sep="\t")
# 
# ## Heterozygosity and F statistics with hierfstat
# sumstats <- basic.stats(cbind(mygl2@pop,as.data.frame(mygl2)))
# sumstats$overall
# pop_summary <- data.frame(N=colMeans(sumstats$n.ind.samp,na.rm=T),
#                           Ho=colMeans(sumstats$Ho,na.rm=T),
#                           He=colMeans(sumstats$Hs,na.rm=T),
#                           Fis=colMeans(sumstats$Fis,na.rm=T))
# write.table(pop_summary,file = file.path(outdir,paste0(analysis_name,"_pop_summary.txt")),sep="\t")
# pop_summary
# 
# ## Take a look at how Ho and He are distributed
# hist(sumstats$perloc$Ho)
# plot(sumstats$perloc$Hs,sumstats$perloc$Ho)

##Fixed difference analysis (script from DArTR manual)-measure of contemporary gene flow, by pair-wise comparison of loci for fixed alleles
fd <- gl.collapse.recursive(gl, t=0)
gl <- fd.sig$gl phy <- gl2phylip(gl, outfile="turtle.phy", bstrap=1000)

#Isolation by distance plot
