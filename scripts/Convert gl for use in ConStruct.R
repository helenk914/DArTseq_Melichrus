######################################################################
## Preparing SubsetAA(M.urceolatus s.s. group) data for ConStruct analysis
######################################################################

######################################################################
## Setup
######################################################################
## Check version of R
# R.version.string
## Install packages 
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
#setwd ()

## Read metadata (read in as metadata0 
## so you can return to manipulating original if needed)

metadata0=read.csv("data/metadata.csv",
                   na.strings = c(""," ","#N/A"))
metadata=metadata0
str(metadata)

######################################################################
## Produce data for conStruct - M.urceolatus s.s. group (AA)
## See script 'Subsetting and filtering groups_secondary' for creation of filtered subset
######################################################################
## Load filtered dataset (check that filtering is optimized for ConStruct)

load("temp/SubFiltglAA.rda",verbose = T)
## Convert from genlight to gen ind format
urc0=gl2gi(gl.filter.monomorphs(Subsetgl1))

## Sort populations 
# oldpop=pop(urc)
# newpop=factor(pop(urc),levels=sort(as.numeric(levels(pop(urc)))))
# pop(urc)=newpop

metadata1 <- metadata[match(indNames(urc0),metadata$id),]
urc0@pop <- factor(metadata1$pop)


## Inspect data structure
urc0
levels(urc0@pop)
head(locNames(urc0))
head(indNames(urc0))

## inspect sample sizes
table(pop(urc0))

urc0$pop
## plot sample sizes
barplot(table(pop(urc0)), col=funky(39),las=3,
        xlab="Population",ylab="Sample size")

##Change name to urc
urc=urc0
table(urc$pop)
urc
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

dim(af.ind.loc)
length(loc.to.keep)
## Filter invariant loci
urc.filt <- urc[,loc=which(loc.to.keep)]
urc.filt$pop
nLoc(urc.filt)

## Convert to genpop and calculate allele frequencies
urc.pop <- genind2genpop(urc.filt)
urc.freq <- makefreq(urc.pop)[,seq(1,2*nLoc(urc.pop)-1,2)]
dim(urc.freq)
urc.freq[,1:10]

# ## Exclude columns with NaN 
col.to.exclude <- sapply(1:ncol(urc.freq),function(x)({
  any(is.nan(urc.freq[,x]))
}))
table(col.to.exclude)
row.to.exclude <- sapply(1:nrow(urc.freq),function(x)({
  any(is.nan(urc.freq[x,]))
}))
nmissing <- rowSums(is.nan(urc.freq))

nmissing

## Filter allele frequencies ######### dont run both 108 and 109 (try both ways, filter out NAN and not)
allele.freqs <- urc.freq[,!col.to.exclude]
allele.frequencies <- allele.freqs[order(rownames(allele.freqs)),]
#allele.frequencies <- urc.freq[order(rownames(urc.freq)),]
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
save(urc.data,file="temp/urcAA.3_construct_pop_data.rda")

