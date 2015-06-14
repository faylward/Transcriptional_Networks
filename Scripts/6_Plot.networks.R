
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
Total.Counts <- as.numeric(apply(canon.counts[1:numgenes,1:numsamples], 2, sum))
Norm=scale(Data.Counts, scale=Total.Counts, center=FALSE)
Norm_Final <- t(Norm)

TOM <- read.table(file="R_network_pipeline/CANON/canon.TC.TOM_11_2014", sep="\t", header=TRUE, row.names=1, nrows=8849, quote="")

########### Load module .module file generated in script 1. 
Module_Colors <- read.table(file="R_network_pipeline/CANON/canon.modules.TC.10_2014", sep="\t", header=TRUE, row.names=1)
mod <- data.frame(sort(table(Module_Colors[,2]), decreasing=T))
mergedColors <- Module_Colors[,2]

########################################################################################
############################# Construct Networks using igraph package ##################
########################################################################################

# Get connectivity data
canon <- read.table(file="R_network_pipeline/CANON/canon.final.combined.list", sep="\t", header=TRUE, row.names=1, quote="")
#connected <- row.names(canon)[which(canon$Scaled.Connectivity > 0.099)]

modules <- c("steelblue1", "purple", "orange", "springgreen4", "red", "magenta", "steelblue4", "brown", "cyan", "turquoise3", "violetred4", "tomato1", "violet", "tan4", "turquoise4", "wheat3", "slategray2", "yellow", "tan1", "sienna2", "darkturquoise", "saddlebrown", "orangered", "firebrick1", "orchid4", "royalblue3", "tomato4", "olivedrab1", "olivedrab4", "springgreen3", "skyblue", "chocolate3", "darkcyan", "aquamarine", "coral", "darkgoldenrod", "blueviolet")

modules =c("steelblue1", "purple", "orange", "springgreen4", "red", "magenta", "steelblue4", "brown", "cyan",  "turquoise3", "violetred4", "tomato1", "violet")
modules =c("steelblue1", "purple", "orange", "springgreen4", "red", "magenta", "steelblue4", "brown", "cyan",  "turquoise3", "gold", "tomato1", "violet")
modules =c("steelblue1", "purple", "orange", "springgreen4", "red", "magenta", "steelblue4", "brown", "cyan",  "turquoise3")
#modules = "gold"
#modules =c("steelblue1", "purple", "orange")

colnames(TOM) <- gsub("X.", "", colnames(TOM)); colnames(TOM) <- gsub(".", "", colnames(TOM), fixed=TRUE)
row.names(TOM) <- colnames(TOM)
nodes <- row.names(canon)[which(canon$Module_Color %in% modules & canon$Scaled.Connectivity.TC > 0.0)]
modTOM <- TOM[nodes, nodes]

### Choose appropriate threshold for network representation
quant = quantile(apply(TOM, 2, as.numeric), 0.993)

cyt = exportNetworkToCytoscape(modTOM, weighted= TRUE, threshold = quant, nodeNames = nodes, nodeAttr = as.character(Module_Colors[nodes,2]))
#cyt = exportNetworkToCytoscape(modTOM,weighted= TRUE, threshold = quant, nodeNames = modProbes, nodeAttr = mergedColors[inModule])

newedges=(cyt$edgeData[,3])
edges <- cbind(as.character(cyt$edgeData[,1]), as.character(cyt$edgeData[,2]), cyt$edgeData[,3])
write.table(edges, file="R_network_pipeline/CANON/canon.edges",  row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)

######### Make Fruchterman-Reingold layout for graph
igr=read.graph("R_network_pipeline/CANON/canon.edges", format="ncol", weights="yes")
set.seed(123)
l <- layout.fruchterman.reingold(igr,niter=500, area=vcount(igr)^3.2)

#########################################################################################
############# Color for cluster-based network ###########################################
#########################################################################################
nodes <- cbind(as.character(cyt$nodeData[,1]), as.character(cyt$nodeData[,3]))
write.table(nodes, file="R_network_pipeline/CANON/canon.nodes", row.names=FALSE, col.names=c("node", "attr"), sep="\t", quote=FALSE)

x=read.table(file="R_network_pipeline/CANON/canon.nodes", sep="\t", header=TRUE)
V(igr)$attr=as.character(x$attr[match(V(igr)$name, x$node)])
#V(igr)$color="grey75"
V(igr)$color=V(igr)$attr

mod_names <- c("Mod 1", "Mod 2", "Mod 3", "Mod 4", "Mod 5", "Mod 6", "Mod 7", "Mod 8", "Mod 9", "Mod 10", "Mod 11", "Mod 12", "Mod 13")
#tiff(file="clusters.tiff", quality=75, res=100, height=30, width=30, units="in")
jpeg(file="canon_module_Network.nokcutoff2.jpg", quality=100, res=600, height=8, width=8, units="in")
#postscript(file="canon_module_Network.3.eps", height=2, width=2)
plot.igraph(igr,vertex.label=NA,layout=l, vertex.size=1.3, edge.color="grey85", vertex.frame.color=NA, edge.width=0.3)
#legend("topleft", legend=mod_names, col=modules, pch=15, cex=0.7, bty="n")
dev.off()

#########################################################################################
######### Color for KEGG ################################################################
#########################################################################################

Kegg_Legend = c("Translation, Chaperones", "Oxidative Phosphorylation", "Amino Acid Met.", "Carbohydrate Met", "Nucleotide Met.", "Cofactor/Vitamin Met.", "Transporters", "DNA Processing/Division", "Protein Export", "Aminoacyl tRNA Biosyn.", "Transcription", "Flagella", "Photosynthesis")
Kegg_Colors = c("firebrick1", "red", "green4", "darkorange1", "gold", "yellow", "blue", "deepskyblue3", "darkolivegreen1", "darkorchid1", "chartreuse", "darkgoldenrod", "cyan")

x=read.table(file="Network_attr_Files/Path/alltaxa.path.nodes", sep="\t", header=TRUE)
V(igr)$attr=as.character(x$attr[match(V(igr)$name, x$node)])
#V(igr)$color="grey75"
V(igr)$color=V(igr)$attr

jpeg(file="Network_attr_Files/General/All_taxa.KEGG.jpg", quality=75, res=150, height=30, width=30, units="in")
plot.igraph(igr,vertex.label=NA, layout=l, vertex.size=1.3, edge.color="grey90", edge.width=0.3, vertex.frame.color=NA)
legend("topleft", legend=Kegg_Legend, col=Kegg_Colors, pch=15, cex=3, bty="n")
dev.off()

#########################################################################################
######### Color for taxa ################################################################
#########################################################################################
library(RColorBrewer)
mypal <- c("green4", brewer.pal(9, "Set1")); mypal[7] = "gold"; mypal[4] <- "#54D126"
mypal <- c("green4", "red", "blue", "yellow", "purple", "steelblue4", "orange", "magenta", "darkolivegreen2", "cyan")
node <- as.character(cyt$nodeData[,1])

taxa <- c("Ostreo_Cluster", "SAR11_Cluster", "SAR86_Cluster", "SAR116_Cluster",  "Roseo_Cluster", "SAR406_Cluster", "SAR92", "Flavo_Cluster", "Arctic_Cluster", "eury_Cluster")
attr <- list(1:length(node))
for(i in 1:length(node)) {
	for(j in 1:length(taxa)) {
		if(grepl(taxa[j], node[i]) == 1) {
			attr[i] <- mypal[j]
		}
	}
}
nodes <- cbind(node, attr)
x <- data.frame(nodes)
#x=read.table(file="Network_attr_Files/Path/nodes.taxa.based", sep="\t", header=TRUE)
V(igr)$attr=as.character(x$attr[match(V(igr)$name, x$node)])
#V(igr)$color="grey75"
V(igr)$color=V(igr)$attr

mod_names <- c("Ostreococcus", "SAR11", "SAR86", "SAR116", "Roseobacter", "SAR406", "SAR92", "Flavobacteria", "ARCTIC96-BD19", "GII Euryarchaeota")
color <- c("red", "green4", "blue", "yellow", "purple", "steelblue4", "orange", "darkolivegreen2", "magenta", "cyan")
#tiff(file="clusters.tiff", quality=75, res=100, height=30, width=30, units="in")
jpeg(file="Canon_Taxa_based_network.no.k.cutoff2.jpg", quality=100, res=600, height=8, width=6, units="in")
plot.igraph(igr,vertex.label=NA,layout=l, vertex.size=1.4, edge.color="grey85", edge.width=0.3, vertex.frame.color=NA)
legend("topleft", legend=mod_names, col=mypal, modules, pch=15, cex=0.5, bty="n")
dev.off()

#########################################################################################
######### Color for KEGG Fewer Pathways #################################################
#########################################################################################
Kegg_Legend = c("Translation, Chaperones", "Oxidative Phosphorylation", "Amino Acid Met.", "Carbohydrate Met", "Nucleotide Met.", "Cofactor/Vitamin Met.", "Transporters", "DNA Processing/Division", "Flagella", "Photosynthesis")
Kegg_Colors = c("brown4", "red", "green4", "darkorange1", "gold", "yellow", "blue", "deepskyblue3", "darkgoldenrod", "cyan")

x=read.table(file="Network_attr_Files/Path/alltaxa.path2.nodes", sep="\t", header=TRUE)
V(igr)$attr=as.character(x$attr[match(V(igr)$name, x$node)])
#V(igr)$color="grey75"
V(igr)$color=V(igr)$attr

jpeg(file="Network_attr_Files/General/All_taxa.KEGG.jpg", quality=75, res=150, height=30, width=30, units="in")
plot.igraph(igr,vertex.label=NA, layout=l, vertex.size=1.3, edge.color="grey90", edge.width=0.3, vertex.frame.color=NA)
legend("topleft", legend=Kegg_Legend, col=Kegg_Colors, pch=15, cex=3, bty="n")
dev.off()



