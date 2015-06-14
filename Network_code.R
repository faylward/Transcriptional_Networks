
########################################################################################################################################################################################### 
## This code was run using R version 3.0.2 (2013-09-25) -- "Frisbee Sailing" on an Ubuntu 14.04.1 LTS. It is recommended at least 4G RAM is available for running this code. 
## This code was used for transcriptional network analyses in Aylward et al., PNAS, 2015. The R packages WGCNA (Langfelder and Horvath, BMC Bioinformatics, 2008) and igraph
## (Csardi and Nepusz, Interjournal, 2006) are required.
## It is assumed the user has a count table (timepoints or treatments as columns, genes/transcripts as rows) to use as input. As an example, the data used in Aylward et al., PNAS, 2015
## provided. To use custom count data the user will need to modify certain parameters to fit the data- for example the soft-thresholding power used in weighted networks will vary for
## dataset. For more information this see Zhang and Horvath, Stat. Appl. Genet. Mol. Biol., 2005.
## For questions regarding this code contact Frank Aylward <frank.o.aylward@gmail.com> 
###########################################################################################################################################################################################

# Load R Packages
library(WGCNA)
library(igraph)

#### Enable multithreading for WGCNA
enableWGCNAThreads()

############ Load raw count data and normalize by total counts in a given timepoint [column]
numsamples <- 35
canon.counts <- read.table(file="Data/canon.combined.counts.10_2014.txt", sep="\t", header=TRUE, row.names=1)
numgenes <- dim(canon.counts)[1]
Data.Counts <- canon.counts[1:numgenes,1:numsamples]
Total.Counts <- as.numeric(apply(canon.counts[2:numgenes,1:numsamples], 2, sum))
Norm=scale(Data.Counts, scale=Total.Counts, center=FALSE)
Norm_Final <- t(Norm)
#time <- as.numeric(canon.counts[1,])

#Pick Soft Threshold
powers =c(c(1:10),seq(from = 12, to=20,by=2))
sft=pickSoftThreshold(Norm_Final, powerVector=powers, verbose=5)

#Document Scale-free topology fits and power-connectivity relationship, output results as a figure
jpeg(file="Data/CANON_SFT_Fit_DESeq_Final.jpg", quality=100, res=400, height=5, width=9, units="in")
par(mfrow =c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signedR^2",type="n",main =paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main =paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
dev.off()

##Power of 5 chosen for the CANON 2012 ESP drifter dataset. The general rule of thumb is to choose the smallest value that gives a scale free index > 0.8.
pow <- 5
net = blockwiseModules(Norm_Final,power= pow,TOMType ="signed", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, verbose = 3) 

#Get unique colors to associate with each module identified. Custom colors were chosen here, but for datasets with variable numbers of modules the user may prefer to use the "labels2colors" function with an unspecified "colorSeq" parameter, as random colors will then be chosen. 
mergedColors = labels2colors(net$colors, colorSeq=c("steelblue1", "purple", "orange", "springgreen4", "red", "magenta", "steelblue4", "brown", "cyan", "turquoise3", "violetred4", "tomato1", "violet", "tan4", "turquoise4", "wheat3", "slategray2", "yellow", "tan1", "sienna2", "darkturquoise", "saddlebrown", "orangered", "firebrick1", "orchid4", "royalblue3", "tomato4", "olivedrab1", "olivedrab4", "springgreen3", "skyblue", "chocolate3", "darkcyan", "aquamarine", "coral", "darkgoldenrod", "blueviolet"))  

#Combine module information (ortholog cluster ID, module number, module color) and output table
#mod <- cbind(colnames(Norm_Final), net$colors, mergedColors); colnames(mod) <- c("Cluster", "Module", "Module_Color")
#write.table(mod, file="Data/canon.modules.10_2014", quote=F, sep="\t", row.names=F)

#Plot dendrogram
#jpeg(file="Data/combined.jpg", quality=75, res=150, height=15, width=30, units="in")
#plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],dendroLabels = FALSE, hang = 0.05,addGuide = TRUE, guideHang = 0.05)
#dev.off()

#calculate and output Topological Overlap Matrix (TOM, for later analyses)
#TOM = TOMsimilarityFromExpr(Norm_Final, power= pow)
colnames(TOM) = colnames(Norm_Final)
write.table(TOM, file="Data/canon.TOM", sep="\t")

#For your viewing pleasure
sort(table(mergedColors), decreasing=T)


################################################################################################
### Identify Eigengenes in Dataset
################################################################################################

#modules <- row.names(data.frame(sort(table(mergedColors), decreasing=T)[1:38]))
# Define modules (by color) for which eigengenes are to be calculated (color codes are identical to those used in script 1). 

modules <- c("steelblue1", "purple", "orange", "springgreen4", "red", "magenta", "steelblue4", "brown", "cyan", "grey40", "gold", "tomato1", "violet", "tan4", "turquoise4", "wheat3", "slategray2", "yellow", "tan1", "sienna2", "darkturquoise", "saddlebrown", "orangered", "firebrick1", "orchid4", "royalblue3", "tomato4", "olivedrab1", "olivedrab4", "springgreen3", "skyblue", "chocolate3", "darkcyan", "aquamarine", "coral", "darkgoldenrod", "blueviolet")

eig <- moduleEigengenes(Norm_Final, as.character(Module_Colors[,2]), excludeGrey=TRUE)
data <- eig$eigengenes
data1 <- eig$averageExpr
#times = as.numeric(x[1,])
timelist = list(time)
#data1 <- cbind(data, times)
colnames(data) <- gsub("ME", "", colnames(data))
#data <- data.frame(scale(data))
### Set up output for eigengene data; the user will need to adjust the plotting parameters
jpeg(file="Data/Eigengenes.jpg", quality=100, res=300, height=8, width=7, units="in")
par(mfrow=c(7,4), mar=c(0.5, 0.5, 1, 0.5), oma=c(2,0,1,0))
for(i in 29:37) {
#for(i in 1:(dim(data)[2]-1)) {
	color <- modules[i]
	#print(color)
	#ag <- aggregate(data[,color], timelist, mean)
	#axis = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24)
	#axis = c(1, 3, 5, 6, 7, 9, 11, 12, 13, 14, 15, 17, 19, 21, 22, 23)
	axis <- c(1:35)
	#row.names(ag) <- axis
	title <- paste(c("Module ", i, "; ", color, "; ", "n=", mod[modules[i],]), collapse="")
	title <- paste(c("Module ", i, "; ", "n=", mod[modules[i],]), collapse="")
	print(title)
	plot(axis, data[,modules[i]], col=color, type="l", ylab=NA, lwd=2, main=title, cex.main=1.2, axes=F)
	rect(xleft = c(0, 5.5, 11.5, 16.5, 23.5, 31.5), xright=c(2.5, 8.5, 13.5, 20.5, 27.5, 34.5), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(1, 1, 1, 1, 1, 1), col='grey90', border=NA)
	par(new=T)
	plot(axis, data[,modules[i]], col='steelblue', type="l", ylab=NA, lwd=2, main=title, cex.main=1.2, axes=F)
		if(i %in% c(34:37)) {axis(1, at = c(0.1, 5, 10, 15, 20, 25, 30, 35), labels = c(0, 5, 10, 15, 20, 25, 30, 35), cex.axis=1) }
	lines(seq(0, 0, length.out=35), pch=22, lty=5)
	box(lty=1, col="grey60")
	#par(new=T)
	#plot(data[,i], col=color, type="l", ylab=NA, lwd=1)
	#par(new=T)
	#plot(times, data[,i], col=color, ylab=NA, ylim=c(-0.2, 0.9))
	#par(new=T)
}
dev.off()



