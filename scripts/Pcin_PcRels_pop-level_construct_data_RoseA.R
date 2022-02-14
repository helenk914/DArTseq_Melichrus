######################################################################
## Population structure of native mints (Prostanthera) in R
######################################################################
## This script is modified from the Module 1.3 activity for use in 
## addressing the first question in the Prostanthera cineolifera project.
## A new script... and possibly a new dataset... may be needed
## for the other questions.

######################################################################
## Setup
######################################################################


## Check version of R
R.version.string

## Install packages - commented out because they only need to be installed once. 
#     The first time, select the first command, but omit the hash (#) character, 
#     and click the "Run" button above. Repeat with other line.
# install.packages("adegenet",dep=T)
# install.packages("hierfstat",dep=T)

## Load packages
library(adegenet)
library(hierfstat)
library(dartR)

######################################################################
## Set wd and load metadata
######################################################################

## Check working directory
getwd()
#setwd("C:/Users/randre20/Dropbox/Ruth Palsson/R-conStruct")

## Read metadata
metadata0=read.csv("data/metadata.csv",
                   na.strings = c(""," ","#N/A"))
metadata=metadata0
str(metadata)

######################################################################
## Produce data for conStruct - M.urceolatus s.s. group (AA)
######################################################################

load("temp/SubFiltglAA.rda",verbose = T)
urc0=gl2gi(Subsetgl1)

## Sort populations 
# oldpop=pop(urc)
# newpop=factor(pop(urc),levels=sort(as.numeric(levels(pop(urc)))))
# pop(urc)=newpop

## Inspect data structure
urc0
levels(urc0@pop)
head(locNames(urc0))
head(indNames(urc0))

## inspect sample sizes
table(pop(urc0))

## plot sample sizes
barplot(table(pop(urc0)), col=funky(39),las=3,
        xlab="Population",ylab="Sample size")

## Make new genind with all populations of species of interest
# poplist=c("cin","lan","ova","haw","oln")
# include=which(indNames(urc0) %in% metadata$id[which(metadata$species_short %in% poplist)])
#urc=urc0[include,]
urc=urc0
table(urc$pop)

## calculate allele frequencies for individuals, then loci
## for filtering out invariant loci
af.ind <- makefreq(urc)
col.index <- c(1:(ncol(af.ind)/2))*2
af.ind.loc <- af.ind[,col.index]
dim(af.ind.loc)
hist(af.ind)
hist(af.ind.loc)
#making a dataframe with allele freq by locus for each pop
af.loc <- colMeans(af.ind.loc,na.rm=TRUE)
str(af.loc)
loc.to.keep <- af.loc > 0.05 & af.loc < 0.95
loc.to.keep[1:20]
table(loc.to.keep,useNA = "always")

## Filter invariant loci
urc.filt <- urc[,loc=which(loc.to.keep)]
urc.filt$pop
nLoc(urc.filt)

## Convert to genpop and calculate allele frequencies
urc.pop <- genind2genpop(urc.filt,pop=)
urc.freq <- makefreq(urc.pop)[,seq(1,2*nLoc(urc.pop)-1,2)]
dim(urc.freq)
urc.freq[,1:10]

# ## Exclude columns with NaN
# col.to.exclude <- sapply(1:ncol(urc.freq),function(x)({
#   any(is.nan(urc.freq[,x]))
# }))
# table(col.to.exclude)

## Filter allele frequencies
# allele.freqs <- urc.freq[,!col.to.exclude]
# allele.frequencies <- allele.freqs[order(rownames(allele.freqs)),]
allele.frequencies <- urc.freq[order(rownames(urc.freq)),]
dim(allele.frequencies)

## Extract coordinates
coords <- metadata[match(unique(urc$pop),metadata$pop),
                   c("decimal_longitude","decimal_latitude")]
coords=as.matrix(coords)

## Geographic distance matrix
geoDist <- dist(coords,diag=T,upper=T)
geoDist <- as.matrix(geoDist)

## Check the formats
str(allele.frequencies)
str(coords)
str(geoDist)

## Bundle together in a list
urc.data <- list(allele.frequencies=allele.frequencies,
                  coords=coords,
                  geoDist=geoDist)
str(urc.data)

## Save to a file that can be loaded from a construct script.
save(urc.data,file="temp/urcAA_construct_pop_data.rda")

