
##################################################################################################################################################
## The following are a series of R commands to take a matrix of count data (timepoints or treatments as columns, genes/transcripts as rows) and calculate the appropriate value for soft-thresholding (which the user must choose based on the results of the scale-free topology and connectivity diagrams). Modules and module membership for each transcript are then calculated. 
#This code was used to analyze the CANON 2012 ESP Drifter dataset of 8848 ortholog clusters from 10 taxa measured over 35 timepoints. 
##################################################################################################################################################
############################## Load Necessary Packages
library(WGCNA)
library(igraph)

############################ Enable multithreading for WGCNA
enableWGCNAThreads()

####################################################################################################### CANON:
############ Load raw count data and normalize by total counts in a given timepoint [column]
numsamples <- 35
canon.counts <- read.table(file="R_network_pipeline/CANON/canon.combined.counts.10_2014.txt", sep="\t", header=TRUE, row.names=1)
numgenes <- dim(canon.counts)[1]
Data.Counts <- canon.counts[1:numgenes,1:numsamples]
Total.Counts <- as.numeric(apply(canon.counts[2:numgenes,1:numsamples], 2, sum))
Norm=scale(Data.Counts, scale=Total.Counts, center=FALSE)
Norm_Final <- t(Norm)
#time <- as.numeric(canon.counts[1,])

#################### Do upper quartile normalization ####################
uq <- c()
for ( i in 1:dim(Data.Counts)[2]) {col <- Data.Counts[,i]; col[col==0] <- NA; uq1 <- quantile(col, probs=0.95, na.rm = TRUE); uq[i] <- uq1 }
Norm=scale(Data.Counts, scale=uq, center=FALSE)
Norm_Final <- t(Norm)
###########################################################################

########### DESEQ
colData <- data.frame(condition=colnames(Norm))
dds <- DESeqDataSetFromMatrix(as.matrix(Data.Counts), colData, formula(~ condition))
var <- varianceStabilizingTransformation(dds, blind=T)
vsd <- t(assay(var))
##########################


#Pick Soft Threshold
powers =c(c(1:10),seq(from = 12, to=20,by=2))
sft=pickSoftThreshold(vsd, powerVector=powers, verbose=5)

#Document Scale-free topology fits and power-connectivity relationship, output results as a figure
jpeg(file="CANON_SFT_Fit_DESeq_Final.jpg", quality=100, res=400, height=5, width=9, units="in")
par(mfrow =c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signedR^2",type="n",main =paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main =paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
dev.off()

##Power of 5 chosen for the CANON 2012 ESP drifter dataset
pow <- 5
net = blockwiseModules(Norm_Final,power= pow,TOMType ="signed", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, verbose = 3) 

#Get unique colors to associate with each module identified. Custom colors were chosen here, but for datasets with variable numbers of modules the user may prefer to use the "labels2colors" function with an unspecified "colorSeq" parameter, as random colors will then be chosen. 
mergedColors = labels2colors(net$colors, colorSeq=c("steelblue1", "purple", "orange", "springgreen4", "red", "magenta", "steelblue4", "brown", "cyan", "turquoise3", "violetred4", "tomato1", "violet", "tan4", "turquoise4", "wheat3", "slategray2", "yellow", "tan1", "sienna2", "darkturquoise", "saddlebrown", "orangered", "firebrick1", "orchid4", "royalblue3", "tomato4", "olivedrab1", "olivedrab4", "springgreen3", "skyblue", "chocolate3", "darkcyan", "aquamarine", "coral", "darkgoldenrod", "blueviolet"))  

#Combine module information (ortholog cluster ID, module number, module color) and output table
mod <- cbind(colnames(Norm_Final), net$colors, mergedColors); colnames(mod) <- c("Cluster", "Module", "Module_Color")
write.table(mod, file="R_network_pipeline/CANON/canon.modules.10_2014", quote=F, sep="\t", row.names=F)

#Plot dendrogram
jpeg(file="combined.jpg", quality=75, res=150, height=15, width=30, units="in")
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],dendroLabels = FALSE, hang = 0.05,addGuide = TRUE, guideHang = 0.05)
dev.off()

#calculate and output Topological Overlap Matrix (TOM, for later analyses)
TOM = TOMsimilarityFromExpr(Norm_Final,power= pow)
colnames(TOM) = colnames(Norm_Final)
write.table(TOM, file="R_network_pipeline/CANON/canon.TOM", sep="\t")

#For your viewing pleasure
sort(table(mergedColors), decreasing=T)

