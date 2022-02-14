################################################################################
## Run 'Set up script for PCOA' first to set up and load mygl (don't need to filter)
## Run 'Subsetting and filtering group...' to produce and save filtered subsets
################################################################################
## Run PCOA for filtered subsets 
## Inspect PCOA
## Adjust graphics and save output for PCOA's

################################################################################
## !!!EDIT PCA name vector or you will save over previous graphs!!!
################################################################################

################################################################################
## Ordination
################################################################################
## Load specific filtered subset "SubFiltgl(AA)" etc. Load metadata and apply to gl. 

load(file = file.path("temp",paste0("SubFiltglAA.rda")))
Subsetgl1

#Read in metadata
metaSubsetgl1 <- metadata0[match(Subsetgl1@ind.names,metadata0$id),]

# Add metadata to gl
Subsetgl1@pop <- factor(metaSubsetgl1$pop)

##Check that populations and names are all there
Subsetgl1@pop
mygl@ind.names

################################################################################
## Principal coordinates analysis 

Subsetpc <- gl.pcoa(Subsetgl1,nfactors = 6)
## Scree plot (represents how the axes represent variation)
# barplot(mypc$eig/sum(mypc$eig)*100)
# gl.pcoa.scree(mypc)

################################################################################
## Plot PCA 
## labels = none, as too many labels cause viewing error

gl.pcoa.plot(Subsetpc,Subsetgl1,labels="none",xaxis=1,yaxis=2)

gl.pcoa.plot(Subsetpc,Subsetgl1,labels="ind",xaxis=1,yaxis=2)
## labels= ind, gives an error, but can be viewed with ggplotly interactive
ggplotly(tooltip=c("pop","ind")) # shows pop twice or not at all

################################################################################
## Save PCA plots 

## Define palette
mypalette <- c("#0ff1ce","#696969","#bada55","#7fe5f0","#ff0000",
               "#ffd700","#407294","#cbcba9","#065535","#5ac18e",
               "#f7347a","#576675","#ffe4e1","#008080","#ffa500",
               "#8a2be2", "#800080","#F5F57A","#F37B7A","#E69F00","#F5F57A",
               "#F37B7A", "#E69F00", "#CC79A7", 
               "#7CB4F7", "#78F5C9", "#D55E00", "#41780f",
               "#d8a6d0","#ec8b3c","#441848","#efe50d")
palette_colourblind <- c("#F5F57A", "#F37B7A", "#E69F00", "#CC79A7", 
                          "#7CB4F7", "#78F5C9", "#D55E00", "#41780f",
                         "#d8a6d0","#ec8b3c","#441848","#efe50d")

## Name, edit and save plot !!EDIT PCA NAME!! 
## Colour code by specific metadata column, bw theme,
## axes proportional to % explained, axes labelled
PCA1 <- ggplot(as.data.frame(Subsetpc$scores), aes(PC1,PC2,col=metaSubsetgl1$scientific_name_OTU))+
  geom_point()+theme_bw()+
  scale_colour_manual(values=mypalette,name = 'OTU')+
  xlab(paste0("PC1 (",round(100*Subsetpc$eig[1]/sum(Subsetpc$eig),1),"%)"))+
  ylab(paste0("PC2 (",round(100*Subsetpc$eig[2]/sum(Subsetpc$eig),1),"%)"))+
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")

PCA2 <-ggplot(as.data.frame(Subsetpc$scores), aes(PC1,PC3,col=metaSubsetgl1$scientific_name_OTU))+
  geom_point()+theme_bw()+scale_colour_manual(values=mypalette,name = 'OTU')+
  xlab(paste0("PC1 (",round(100*Subsetpc$eig[1]/sum(Subsetpc$eig),1),"%)"))+
  ylab(paste0("PC3 (",round(100*Subsetpc$eig[3]/sum(Subsetpc$eig),1),"%)"))+
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")

PCA3 <-ggplot(as.data.frame(Subsetpc$scores), aes(PC2,PC3,col=metaSubsetgl1$scientific_name_OTU))+
  geom_point()+theme_bw()+scale_colour_manual(values=mypalette,name = 'OTU')+
  xlab(paste0("PC2 (",round(100*Subsetpc$eig[2]/sum(Subsetpc$eig),1),"%)"))+
  ylab(paste0("PC3 (",round(100*Subsetpc$eig[3]/sum(Subsetpc$eig),1),"%)"))+
  coord_fixed(ratio = 1, xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")

## Save as png (easy to copy and paste) and/or pdf (better resolution)
#!!
PCAname = paste0("SubsetpcAA")
#!!
png(file.path(outdir,paste0(PCAname, "1.png")));print(PCA1);dev.off()
png(file.path(outdir,paste0(PCAname, "2.png")));print(PCA2);dev.off()
png(file.path(outdir,paste0(PCAname, "3.png")));print(PCA3);dev.off()
pdf(file.path(outdir,paste0(PCAname, "1.pdf")));print(PCA1);dev.off()
pdf(file.path(outdir,paste0(PCAname, "2.pdf")));print(PCA2);dev.off()
pdf(file.path(outdir,paste0(PCAname, "3.pdf")));print(PCA3);dev.off()
dev.off()

