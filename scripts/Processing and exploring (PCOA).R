##################################################################
##################################################################
##R.Andrew Original script
##Modified by H.Kennedy 2021 
##Processing, exploration of DArTseq data by PCOA
##################################################################
##################################################################
## Requires *SNP_2.csv file supplied by DArT and matching sample
## *metadata.csv in a folder specific to the DArT run (within 
## the data directory). 
##################################################################
##Products are; mygl, unfiltered genlight object, mygl1, filtered gl
## for PCOA
## Populations defined under 'pop' in metadata
#help(package = 'base')
##################################################################
##Run 'set up script for PCOA first'
##################################################################

##################################################################
## Ordination
##################################################################
#load(file = file.path("temp",paste0(analysis_name,"_filtered1.rda")))
#mygl1
##################################################################
## Principal coordinates analysis 
##################################################################
mypc <- gl.pcoa(mygl1,nfactors = 6)
## Scree plot (represents how the axes represent variation)
# barplot(mypc$eig/sum(mypc$eig)*100)
# gl.pcoa.scree(mypc)

## Plot PCA ####################################################
## labels = none, as too many labels cause viewing error
gl.pcoa.plot(mypc,mygl1,labels="none",xaxis=1,yaxis=2)

gl.pcoa.plot(mypc,mygl1,labels="ind",xaxis=1,yaxis=2)
## labels= ind, gives an error, but can be viewed with ggplotly interactive
ggplotly(tooltip=c("pop","ind")) # shows pop twice or not at all

gl.pcoa.plot(mypc,mygl1,labels="ind",xaxis=1,yaxis=3)
ggplotly(tooltip=c("pop","ind")) # shows pop twice or not at all

gl.pcoa.plot(mypc,mygl1,labels="ind",xaxis=2,yaxis=3)
ggplotly(tooltip=c("pop","ind")) # shows pop twice or not at all

#################################################################
## Save PCOA Graphs
#################################################################

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
##gl.pcoa.plot.3d(mypc, mygl1)

##################################################################
## Neighbour-joining tree
##################################################################
## By population
mygl1nj <- gl.tree.nj(mygl1, type="fan")
plot.phylo
gl.tree.nj(mygl1, type="phylogram")
##Need to show a longer label...but how?


