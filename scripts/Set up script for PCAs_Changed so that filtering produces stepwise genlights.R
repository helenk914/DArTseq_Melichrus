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
analysis_name <- "MELICHRUS"
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
mygl <- gl.read.dart(dartfile)
save (mygl, file= "MELICHRUS.rdata")
## Inspect gl
mygl

##################################################################
## Open saved gl file 
################################################################
load("MELICHRUS.rdata")
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
gl.report.callrate(mygl,method = "ind")
# ##  Distances between 'tags' - identify oversplitting of loci 
# gl.report.hamming(mygl) # requires too much memory
# gl.report.hamming(mygl[,1:10],rs = 5)
# ##  Average repeatability of each snp
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
  ## Remove any accidental duplicates in dataset (HTK96c.1 and HTK96e.1)
  mygl2 <- gl.drop.ind(mygl, ind.list=c("HTK96c.1", "HTK96e.1", "ARC870b"), recalc = TRUE, mono.rm = TRUE)
  ## Average repeatability
  mygl3 <- gl.filter.RepAvg(mygl2,threshold = .98)
  mygl3
  ## Remove monomorphic loci
  mygl4 <- gl.filter.monomorphs(mygl3)
  mygl4
  ## Filter out large amounts of missing data by locus
  mygl5 <- gl.filter.callrate(mygl4,method = "loc",threshold = 0.95)
  mygl5
  ##  Thin to one snp per tag. Usually best to choose the one with the highest minor allele frequency.
  mygl6 <- gl.filter.secondaries(mygl5,method = "best")
  mygl6
  ## Filter out individuals with too much missing data
  mygl7 <- gl.filter.callrate(mygl6,method = "ind",threshold = 0.9,recalc = TRUE, mono.rm = TRUE)
  mygl7
  ##  Filter out loci with very low minor allele frequencies
  ##  exclude step for pop stats as this biases statistics like heterozygosities
  ##  and analyses based on allele frequency spectra.
  mygl8 <- gl.filter.maf(mygl7,threshold = 0.05)
  mygl8
  mygl9 <- gl.filter.monomorphs(mygl8)
  mygl9
  ## Recalculate all of the metrics for your filtered genlight
  mygl10 <- gl.recalc.metrics(mygl9)
  mygl10
  mygl1 <- mygl10
  
  nLoc(mygl1);nInd(mygl1)
  nLoc(mygl);nInd(mygl)
  
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
metadata1 <- metadata0[match(mygl6@ind.names,metadata0$id),]
metadata1

##################################################################
## Ordination
##################################################################
load(file = file.path("temp",paste0(analysis_name,"_filtered1.rda")))
mygl6
## Principal coordinates analysis ################################
mypc <- gl.pcoa(mygl6,nfactors = 6)
## Scree plot (represents how the axes represent variation)
# barplot(mypc$eig/sum(mypc$eig)*100)
# gl.pcoa.scree(mypc)

## Plot PCA ####################################################
## labels = none, as too many labels cause viewing error
gl.pcoa.plot(mypc,mygl6,labels="none",xaxis=1,yaxis=2)

gl.pcoa.plot(mypc,mygl6,labels="ind",xaxis=1,yaxis=2)
## labels= ind, gives an error, but can be viewed with ggplotly interactive
ggplotly(tooltip=c("pop","ind")) # shows pop twice or not at all

gl.pcoa.plot(mypc,mygl6,labels="ind",xaxis=1,yaxis=2)
ggplotly(tooltip=c("pop","ind")) # shows pop twice or not at all

gl.pcoa.plot(mypc,mygl1,labels="ind",xaxis=2,yaxis=3)
ggplotly(tooltip=c("pop","ind")) # shows pop twice or not at all

## Save PCA ####################################################

##Name, edit and save plot ######################################
##Colour code by specific metadata column, bw theme, axes proportional to % explained, axes labelled
GenuswidePCA1 <- ggplot(as.data.frame(mypc$scores), aes(PC1,PC2,col=metadata1$sens.lat.complex))+
  geom_point()+theme_bw()+xlab(paste0("PC1 (",round(100*mypc$eig[1]/sum(mypc$eig),1),"%)"))+
  ylab(paste0("PC2 (",round(100*mypc$eig[2]/sum(mypc$eig),1),"%)"))+
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")

GenuswidePCA2 <-ggplot(as.data.frame(mypc$scores), aes(PC1,PC3,col=metadata1$scientific_name_OTU))+
  geom_point()+theme_bw()+xlab(paste0("PC1 (",round(100*mypc$eig[1]/sum(mypc$eig),1),"%)"))+
  ylab(paste0("PC3 (",round(100*mypc$eig[3]/sum(mypc$eig),1),"%)"))+
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")

GenuswidePCA3 <-ggplot(as.data.frame(mypc$scores), aes(PC2,PC3,col=metadata1$scientific_name_OTU))+
  geom_point()+theme_bw()+xlab(paste0("PC2 (",round(100*mypc$eig[2]/sum(mypc$eig),1),"%)"))+
  ylab(paste0("PC3 (",round(100*mypc$eig[3]/sum(mypc$eig),1),"%)"))+
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")

## To control group colour add this ##############################
#scale_colour_manual(values=mypalette,name = 'Species complex')
##################################################################

## Save as png (easy to copy and paste) and/or pdf (better resolution)
png(file.path(outdir,paste0(analysis_name,"_pcoa_GenuswidePCA1.png")));print(GenuswidePCA1);dev.off()
png(file.path(outdir,paste0(analysis_name,"_pcoa_GenuswidePCA2.png")));print(GenuswidePCA2);dev.off()
png(file.path(outdir,paste0(analysis_name,"_pcoa_GenuswidePCA3.png")));print(GenuswidePCA3);dev.off()
dev.off()
pdf(file.path(outdir,paste0(analysis_name,"_pcoa_GenuswidePCA1.pdf")));print(GenuswidePCA1);dev.off()
pdf(file.path(outdir,paste0(analysis_name,"_pcoa_GenuswidePCA2.pdf")));print(GenuswidePCA2);dev.off()
pdf(file.path(outdir,paste0(analysis_name,"_pcoa_GenuswidePCA3.pdf")));print(GenuswidePCA3);dev.off()

##################################################################
## 3D ordination plot
##################################################################
## Use if the 3rd axes explain > 10
gl.pcoa.plot.3d(mypc, mygl6)

##################################################################
## Neighbour-joining tree
##################################################################
## By population
mygl1nj <- gl.tree.nj(mygl1, type="fan")
plot.phylo
gl.tree.nj(mygl1, type="phylogram")
##Need to show a longer label...but how?
