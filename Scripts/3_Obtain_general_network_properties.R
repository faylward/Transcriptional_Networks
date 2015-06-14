#########################################################################################################
############################ Obtain General Network Properties [Takes a Super Long Time!] ###############
#########################################################################################################
library(WGCNA)
library(igraph)

############################ Enable multithreading for WGCNA
enableWGCNAThreads()

#######################################################################################################Alltaxa:
############ Load raw count data and normalize by total counts in a given timepoint [column]
numsamples <- 35
Alltaxa.counts <- read.table(file="R_network_pipeline/CANON/canon.combined.counts.10_2014.txt", sep="\t", header=TRUE, row.names=1)
numgenes <- dim(Alltaxa.counts)[1]
Data.Counts <- Alltaxa.counts[1:numgenes,1:numsamples]
Total.Counts <- as.numeric(apply(Alltaxa.counts[1:numgenes,1:numsamples], 2, sum))
Norm=scale(Data.Counts, scale=Total.Counts, center=FALSE)
Norm_Final <- t(Norm)
#time <- as.numeric(Alltaxa.counts[1,])

# Load data generated in script 1
TOM <- as.matrix(read.table(file="R_network_pipeline/CANON/canon.TOM", sep="\t", header=TRUE, row.names=1, quote=''))

# This is the long step- generate all statistics (connectivity, etc)
conf <- conformityBasedNetworkConcepts(TOM)
k = conf$fundamentalNCs$Connectivity
stats <- cbind(conf$fundamentalNCs$Connectivity, conf$fundamentalNCs$ScaledConnectivity, conf$fundamentalNCs$MAR, conf$fundamentalNCs$ClusterCoef, conf$approximateConformityBasedNCs$Connectivity.CF.App, conf$approximateConformityBasedNCs$MAR.CF.App, conf$approximateConformityBasedNCs$ClusterCoef.CF.App)
colnames(stats)=c("Connectivity", "Scaled Connectivity", "MAR", "ClusteringCoef", "ConnectivityCF", "MAR.CF.App", "ClusteringCoef.CF.App")
row.names(stats)=colnames(Norm_Final)

#output data
write.table(stats, file="canon_Network_att.TC.final", sep="\t", quote=FALSE)
data=c(conf$fundamentalNCs$Density, conf$fundamentalNCs$Centralization, conf$fundamentalNCs$Heterogeneity, conf$approximateConformityBasedNCs$Density.CF.App, conf$approximateConformityBasedNCs$Centralization.CF.App, conf$approximateConformityBasedNCs$Heterogeneity.CF.App)
names <- c("Density", "Centralization", "Heterogeneity", "Density.CF.App", "Centralization.CF.App", "Heterogeneity.CF.App")
dataset <- cbind(names, data)
write.table(dataset, file="canon_att1.TC.final",  row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

###############################################################################################################

########### DESEQ
colData <- data.frame(condition=colnames(Norm))
dds <- DESeqDataSetFromMatrix(as.matrix(Data.Counts), colData, formula(~ condition))
var <- varianceStabilizingTransformation(dds, blind=T)
vsd <- t(assay(var))
TOM = TOMsimilarityFromExpr(vsd, power= 6)
colnames(TOM) = colnames(vsd)

conf <- conformityBasedNetworkConcepts(TOM)
k = conf$fundamentalNCs$Connectivity
stats <- cbind(conf$fundamentalNCs$Connectivity, conf$fundamentalNCs$ScaledConnectivity, conf$fundamentalNCs$MAR, conf$fundamentalNCs$ClusterCoef, conf$approximateConformityBasedNCs$Connectivity.CF.App, conf$approximateConformityBasedNCs$MAR.CF.App, conf$approximateConformityBasedNCs$ClusterCoef.CF.App)
colnames(stats)=c("Connectivity", "Scaled Connectivity", "MAR", "ClusteringCoef", "ConnectivityCF", "MAR.CF.App", "ClusteringCoef.CF.App")
row.names(stats)=colnames(Norm_Final)
write.table(stats, file="canon_Network_att.deseq", sep="\t", quote=FALSE)

data=c(conf$fundamentalNCs$Density, conf$fundamentalNCs$Centralization, conf$fundamentalNCs$Heterogeneity, conf$approximateConformityBasedNCs$Density.CF.App, conf$approximateConformityBasedNCs$Centralization.CF.App, conf$approximateConformityBasedNCs$Heterogeneity.CF.App)
names <- c("Density", "Centralization", "Heterogeneity", "Density.CF.App", "Centralization.CF.App", "Heterogeneity.CF.App")
dataset <- cbind(names, data)
write.table(dataset, file="canon_att1.deseq",  row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

