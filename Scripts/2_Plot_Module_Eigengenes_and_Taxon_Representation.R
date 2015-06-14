
#######################################################################################
################################ Calculate Eigengenes #################################
#######################################################################################

##################################################################################################################################################
############################## Load Necessary Packages
library(WGCNA)
library(igraph)

############################ Enable multithreading for WGCNA
enableWGCNAThreads()

#############################################################################################################################
############ Load raw count data and normalize by total counts in a given timepoint [column] [this is identical to script 1].
numsamples <- 35
canon.counts <- read.table(file="R_network_pipeline/CANON/canon.combined.counts.10_2014.txt", sep="\t", header=TRUE, row.names=1)
numgenes <- dim(canon.counts)[1]
Data.Counts <- canon.counts[1:numgenes,1:numsamples]
Total.Counts <- as.numeric(apply(canon.counts[1:numgenes,1:numsamples], 2, sum))
Norm=scale(Data.Counts, scale=Total.Counts, center=FALSE)
Norm_Final <- t(Norm)

########### Load module .module file generated in script 1. 
Module_Colors <- read.table(file="R_network_pipeline/CANON/canon.modules.TC.10_2014", sep="\t", header=TRUE, row.names=1)
mod <- data.frame(sort(table(Module_Colors[,2]), decreasing=T))

#modules <- row.names(data.frame(sort(table(mergedColors), decreasing=T)[1:38]))
# Define modules (by color) for which eigengenes are to be calculated (color codes are identical to those used in script 1). 
modules <- c("steelblue1", "purple", "orange", "springgreen4", "red", "magenta", "steelblue4", "brown", "cyan", "turquoise3", "violetred4", "tomato1", "violet", "tan4", "turquoise4", "wheat3", "slategray2", "yellow", "tan1", "sienna2", "darkturquoise", "saddlebrown", "orangered", "firebrick1", "orchid4", "royalblue3", "tomato4", "olivedrab1", "olivedrab4", "springgreen3", "skyblue", "chocolate3", "darkcyan", "aquamarine", "coral", "darkgoldenrod", "blueviolet")
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
jpeg(file="canon.eigengenes.3.jpg", quality=100, res=300, height=8, width=7, units="in")
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


############################################################################################################
########################## Calculate Taxa Representation Across the Modules Identified #####################
############################################################################################################

############### Set color palette using RColorBrewer
library(RColorBrewer)
mypal <- c("green4", brewer.pal(9, "Set1"))

###################################################
num <- 34
result <- c()
names <- c("Ostreococcus", "SAR11", "SAR116", "SAR86", "GII Euryarchaeota", "SAR92", "Flavobacteria", "Roseobacter", "ARCTIC96-BD19", "SAR406")
for(i in 0:num) {

	if(i == 0) {
	set <- subset(Module_Colors, grepl("grey", Module_Colors[,2]))
	subset <- as.character(row.names(set))
	}
	else {
	subset <- row.names(Module_Colors)[which(Module_Colors[,1] == i)]
	}
	string <- strsplit(subset, "Cluster")
	list <- c("Ostreo_"=0, "SAR11_" = 0, "SAR116_" = 0, "SAR86_"=0, "eury_"=0, "SAR92_"=0, "Flavo_"=0, "Roseo_"=0, "Arctic_"=0, "SAR406_"=0)
		for(j in 1:length(string)) {
			list[string[[j]][1]] = list[string[[j]][1]] + 1
		}
	final <- as.numeric(list)
	result <- data.frame(rbind(result, final))
}
colnames(result) <- names
row.names(result) <- c("Unassigned", 1:34)
sum <- apply(result, 1, sum)
result[,11] <- sum
order <- result[2:35,][order(result[,11], decreasing=TRUE),]

library(ggplot2)
library(reshape2)
dat1 <- melt(t(result[1:10]))
colnames(dat1) <- c("Taxon", "Module", "Transcripts")
jpeg(file="CANON_Modules.jpg", quality=100, res=500, height=8, width=8, units="in")
ggplot(dat1, aes(Module, Transcripts, fill=Taxon), guide=guide_legend(title="size")) + geom_bar(stat="identity") + coord_flip() + scale_fill_hue(l=50, c=60) + theme(axis.text.y  = element_text(colour="black", size=10, face="bold"), axis.text.x = element_text(colour="black", face="bold"), legend.title = element_text(size=12, face="bold"), legend.text = element_text(size = 10)) + scale_fill_manual(values=mypal)
dev.off()


