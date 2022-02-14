##################################################################
##################################################################
## Processing and exploration of DArTseq data
##################################################################
##################################################################
## In-progress script for preliminary processing and exploration,
## not final - adjust as appropriate for your data. 
##
## Requires *SNP_2.csv file supplied by DArT and matching sample
## *metadata.csv in a folder specific to the DArT run (within 
## the data directory). 
##
##################################################################
## Load packages
##################################################################

library(dartR)
library(StAMPP)
library(gplots) 
library(plotly)
library(hierfstat)
library(poppr)

##################################################################
## Setup
##################################################################
## Specify main parameters
setwd("C:/Users/hkenned6/Documents/Dartseq")
overwrite <- TRUE
## Directories: choose an indir and an analysis name
# indir <- "data/OrderAppendix_3_DPro19-4334" ## PgRels
# analysis_name <- "PgRels"
indir <- "data/Report-DMano19-4611/" ## PgRels
analysis_name <- "NoisyMiner"

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
if(grepl("SNP_2",dartfile)){
  mygl <- gl.read.dart.2row(datafile = dartfile,
                            topskip = 6,nmetavar = 26)
} else {
  mygl <- gl.read.dart(dartfile)
}

## Read metadata file and sort according to sample order in dartfile
metadata0 <- read.csv(metadatafile, stringsAsFactors = F)
metadata <- metadata0[match(mygl@ind.names,metadata0$id),]

## Add sample metadata to mygl
mygl@pop <- factor(metadata$pop)
table(mygl@pop)

##################################################################
## Filter data
##################################################################
## Explore data quality
##  Base frequencies across all sites, plus proportions of snps by type:
##  transitions (a<>g or c<>t) and transversions (all other snps)
gl.report.bases(mygl)
##  Missing data: proportion of snps with X missing genotypes
gl.report.callrate(mygl)
##  Missing data: proportion of individuals with X missing genotypes
gl.report.callrate(mygl,method = "ind")
##  Distances between 'tags' - identify oversplitting of loci
# gl.report.hamming(mygl) # requires too much memory
gl.report.hamming(mygl[,1:10],rs = 5)
##  Average repeatability of each snp
gl.report.repavg(mygl)
##  Secondary snps occur when there are multiple snps in a single tag.
##  They will be in linkage disequilibrium, so we will thin the data to 
##  one snp per tag. Some may still be linked, but not many. 
gl.report.secondaries(mygl)
##  Monomorphic sites should be filtered out. DArT will have done this
##  before the data are provided to us, but if we subsample it, we 
##  need to do that again.
gl.report.monomorphs(mygl)
##  Check LD - a few for illustration
# gl.report.ld(mygl) # takes a while without parallelisation
gl.report.ld(mygl[,1:10])
##  HWE testing needs a population specification
# mygl@pop=factor(rep(1,nInd(mygl)))
# gl.report.hwe(mygl) # produces a large dataframe
gl.report.hwe(mygl[,1:10])
##  Look at distribution of observed heterozygosity, as it can indicate
##  paralogous loci or clonal reproduction. 
hist(gl.Ho(mygl))

## Filter data - using the same object to limit memory usage
mygl1 <- gl.filter.repavg(mygl,threshold = 0.95)
nLoc(mygl1);nInd(mygl1)
mygl1 <- gl.filter.callrate(mygl1,method = "loc",threshold = 0.95)
mygl1 <- gl.filter.callrate(mygl1,method = "ind",threshold = 0.9)
nLoc(mygl1);nInd(mygl1)
##  Filtering out over-split loci is a good idea, but too time consuming
##  for us today.
# mygl1 <- gl.filter.hamming(mygl1,threshold = 0.2,rs = 5)
# nLoc(mygl1);nInd(mygl1)
##  If your species reproduces sexually, it is unlikely that observed 
##  heterozygosity will be > 0.5 - how unlikely will depend on the sample 
##  size, so choose a threshold based on your gl.report.heterozygosity.
mygl1 <- mygl1[,gl.Ho(mygl1)<0.8]
nLoc(mygl1);nInd(mygl1)
##  In most circumstances, you'll want to thin to one snp per tag.
##  It will usually make sense to select the one with the highest minor allele frequency 
mygl1 <- gl.filter.secondaries(mygl1,method = "best")
nLoc(mygl1);nInd(mygl1)
##  It often makes sense to filter out loci with very low minor allele 
##  frequencies, but bear in mind this biases statistics like heterozygosities
##  and analyses based on allele frequency spectra.
mygl1 <- gl.filter.maf(mygl1,threshold = 0.05)
nLoc(mygl1);nInd(mygl1)

## Recalculate all of the metrics for your filtered genlight
mygl1 <- gl.recalc.metrics(mygl1)
## It's often useful to save the population IDs as a vectore
pop0 <- mygl1@pop
## To save time filtering in the future, save your filtered data as a
## compressed R data object.
save(mygl1,file = file.path("temp",paste0(analysis_name,"_filtered.rda")))

##################################################################
## Population-level statistics
##################################################################

## Compare heterozygosity among populations
ho_pop <- gl.report.heterozygosity(mygl1)
pdf(file.path(outdir,paste0(analysis_name,"_ho_pop.pdf")));gl.report.heterozygosity(mygl1);dev.off()
write.table(ho_pop,file = file.path(outdir,paste0(analysis_name,"_ho_pop.txt")))
# gl.report.maf(mygl) # takes a lot of time

## Examine private alleles by pairs of populations
pa_pw_pop <- gl.report.pa.pop(mygl1)
pa_pw_pop
write.table(pa_pw_pop,file = file.path(outdir,paste0(analysis_name,"_pa_pw_pop.txt")))

## Genetic distances calculated in dartR
mydist1 <- gl.dist.pop(mygl1,method = "pcfixed"); mydist1
mydist2 <- gl.dist.pop(mygl1,method = "pa"); mydist2
mydist3 <- gl.dist.pop(mygl1,method = "gower",upper = T); mydist3 # Why the warning??
mydist4 <- gl.dist.pop(mygl1,method = "euclidean",upper = T); mydist4
write.table(mydist1,file = file.path(outdir,paste0(analysis_name,"_pop_dist_pcfixed.txt")),sep="\t")
write.table(mydist2,file = file.path(outdir,paste0(analysis_name,"_pop_dist_pa.txt")),sep="\t")
write.table(as.matrix(mydist3),file = file.path(outdir,paste0(analysis_name,"_pop_dist_gower.txt")),sep="\t")
write.table(as.matrix(mydist4),file = file.path(outdir,paste0(analysis_name,"_pop_dist_euclidean.txt")),sep="\t")

## Nei's (1972) genetic distance
mygeno <- stamppConvert(mygl1,type = "genlight")
mydist5 <- stamppNeisD(mygeno,pop=T)
write.table(mydist5,file = file.path(outdir,paste0(analysis_name,"_pop_dist_Nei1972.txt")),sep="\t")

## F statistics with STAMPP
pwfst <-stamppFst(mygl1, nboots=1000, percent=95, nclusters=1)
pwfst$Fsts
pwfst$Pvalues
write.table(pwfst$Fsts,file = file.path(outdir,paste0(analysis_name,"_pop_fst.txt")),sep="\t")

## Heterozygosity and F statistics with hierfstat
sumstats <- basic.stats(cbind(mygl1@pop,as.data.frame(mygl1)))
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

##################################################################
## Ordination
##################################################################
load(file = file.path("temp",paste0(analysis_name,"_filtered.rda")))

## Principal coordinates analysis
mypc <- gl.pcoa(mygl1,nfactors = 6)
## Scree plot
barplot(mypc$eig/sum(mypc$eig)*100)
gl.pcoa.scree(mypc)

gl.pcoa.plot(mypc,mygl1,labels="ind",xaxis=1,yaxis=2)
plot1 <- gl.pcoa.plot(mypc,mygl1,labels="pop",xaxis=1,yaxis=2)
plot2 <- gl.pcoa.plot(mypc,mygl1,labels="pop",xaxis=1,yaxis=3)
plot3 <- gl.pcoa.plot(mypc,mygl1,labels="pop",xaxis=1,yaxis=4)
plot4 <- gl.pcoa.plot(mypc,mygl1,labels="pop",xaxis=1,yaxis=5)
plot5 <- gl.pcoa.plot(mypc,mygl1,labels="pop",xaxis=1,yaxis=6)
plot6 <- gl.pcoa.plot(mypc,mygl1,labels="pop",xaxis=2,yaxis=3)
pdf(file.path(outdir,paste0(analysis_name,"_pcoa_12.pdf")));plot1;dev.off()
pdf(file.path(outdir,paste0(analysis_name,"_pcoa_13.pdf")));plot2;dev.off()
pdf(file.path(outdir,paste0(analysis_name,"_pcoa_14.pdf")));plot3;dev.off()
pdf(file.path(outdir,paste0(analysis_name,"_pcoa_15.pdf")));plot4;dev.off()
pdf(file.path(outdir,paste0(analysis_name,"_pcoa_16.pdf")));plot5;dev.off()
pdf(file.path(outdir,paste0(analysis_name,"_pcoa_23.pdf")));plot6;dev.off()

gl.pcoa.plot(mypc,mygl1,labels="ind",xaxis=1,yaxis=2)
ggplotly(tooltip=c("pop","ind")) # shows pop twice or not at all

## Try making a 3d plot - may not work if necessary packages aren't installed
gl.pcoa.plot.3d(mypc, mygl1)

##################################################################
## Neighbour-joining tree
##################################################################
## By population
gl.tree.nj(mygl1, type="fan")
gl.tree.nj(mygl1, type="phylogram")





