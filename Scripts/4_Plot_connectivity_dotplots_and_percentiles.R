

################################################################################################################################################
########### Plot Connectivity Percentiles and Rank-ordered Dotplots for CANON and BioLincs ESP Drifter Datasets ################################
################################################################################################################################################
######### The input file for these analyses is a connectivity file containing i) all nodes in a network with taxon-specific IDs (SAR11_Cluster_1, etc), and ii) the associated connectivity scores for each cluster. 
#This script uses the R packages "RColorBrewer" and "ggplot2".
 
############### Set color palette using RColorBrewer
library(RColorBrewer)
mypal <- c("green4", brewer.pal(9, "Set1")); mypal[7] = "gold"; mypal[4] <- "#54D126"

#######################################################################
############## Load BioLincs data and plot connectivity dotplot #######
#######################################################################

k <- read.table(file="R_network_pipeline/BioLincs/blc_Network_att.new", header=TRUE, row.names=1, sep="\t")
taxa <- c("Pro_Cluster", "SAR11_Cluster", "SAR116_Cluster", "SAR86_Cluster", "Roseo_Cluster", "SAR406_Cluster", "SAR324_Cluster")
jpeg(file="BLC.connectivity.TC.jpg", quality=100, res=300, height=7, width=10, units="in")
for(i in 1:length(taxa)) {
	clusters <- subset(k, grepl(taxa[i], row.names(k)))
	clusters <- clusters[order(clusters$Connectivity, decreasing=T),]
	plot(clusters$Scaled.Connectivity[1:1000], col=mypal[i], ylim = c(0, 1), xlim=c(1, 700), ylab="Scaled Connectivity", xlab="Connectivity-Ordered Transcripts", pch=16, cex.lab=1.2, font.lab=2, axes=F)
	par(new=T)
}
axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1), labels=c('0', '0.2', '0.4', '0.6', '0.8', '1'), cex.axis=1.3, las=2, font.axis='2', cex.axis=1.2)
axis(1, at=c(0, 100, 200, 300, 400, 500, 600, 700), labels=c('0', '100', '200', '300', '400', '500', '600','700'), cex.axis=1.2, font.axis='2')
labels <- c("Prochlorococcus", "SAR11", "SAR116", "SAR86", "Roseobacter", "SAR406", "SAR324")
legend("topright", legend=labels, col=mypal[1:7], pch=15, cex=1.3, bty="n")
box(lwd=2)
dev.off()

#########################################################################
############# Plot percentile bar graph  ################################
#########################################################################

name <- c("Prochlorococcus", "SAR11", "SAR116", "SAR86", "Roseobacter", "SAR406", "SAR324")
Result <- data.frame()
sequ <- seq(0, 1, by=0.05)
#sequ <- seq(0.9, 1, by=0.005)
	for(i in 2:length(sequ)) {
	quant_right <- quantile(k$Connectivity, sequ[i])
	quant_left <- quantile(k$Connectivity, sequ[i-1])
	names <- as.character(row.names(k)[which(k$Connectivity > quant_left & k$Connectivity < quant_right )])
	string <- strsplit(names, "Cluster")
	list <- c("Pro_"=0, "SAR11_" = 0, "SAR116_" = 0, "SAR86_"=0, "Roseo_"=0, "SAR406_"=0, "SAR324_" = 0)
		for(j in 1:length(string)) {
			list[string[[j]][1]] = list[string[[j]][1]] + 1
		}
	final <- as.numeric(list)
	Result <- data.frame(rbind(Result, final))
	}
colnames(Result) <- name
row.names(Result) <- 100*(sequ[2:(length(sequ))] - 0.025)
library(ggplot2)
library(reshape2)
dat1 <- melt(t(Result))
colnames(dat1) <- c("Taxon", "Percentile", "Proportion")
#cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
jpeg(file="BLC_quantiles_90_100.jpg", quality=100, res=300, height=4, width=8, units="in")
ggplot(dat1, aes(Percentile, Proportion, fill=Taxon), guide=guide_legend(title="size")) + geom_bar(stat="identity", position="fill") + coord_flip() + scale_fill_hue(l=50, c=60) + theme(legend.title = element_text(size=12, face="bold"), legend.text = element_text(size = 10)) + scale_fill_manual(values=mypal)
dev.off()

##########################################################################################
########################### Load CANON data and plot connectivity ########################
##########################################################################################
######## For connectivity dotplot

k <- read.table(file="R_network_pipeline/CANON/canon_Network_att.new", header=TRUE, row.names=1, sep="\t")
taxa <- c("Ostreo_Cluster", "SAR11_Cluster", "SAR116_Cluster", "SAR86_Cluster", "Roseo_Cluster", "SAR406_Cluster", "eury_Cluster", "Flavo_Cluster", "Arctic_Cluster", "SAR92")
jpeg(file="CANON.connectivity.TC.jpg", quality=100, res=300, height=7, width=10, units="in")
for(i in 1:length(taxa)) {
	
	clusters <- subset(k, grepl(taxa[i], row.names(k)))
	clusters <- clusters[order(clusters$Connectivity, decreasing=T),]
	clus <- row.names(clusters)[1:100]
	plot(clusters$Scaled.Connectivity[1:1000], col=mypal[i], ylim = c(0, 1), xlim=c(1, 700), ylab="Scaled Connectivity", xlab="Connectivity-Ordered Transcripts", pch=16, cex.lab=1.2, font.lab=2, axes=F)
	par(new=T)
}
axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1), labels=c('0', '0.2', '0.4', '0.6', '0.8', '1'), cex.axis=1.3, las=2, font.axis='2', cex.axis=1.2)
axis(1, at=c(0, 100, 200, 300, 400, 500, 600, 700), labels=c('0', '100', '200', '300', '400', '500', '600','700'), cex.axis=1.2, font.axis='2')
labels <- c("Ostreococcus", "SAR11", "SAR116", "SAR86", "Roseobacter", "SAR406", "GII Euryarchaeota", "Flavobacteria", "ARCTIC96-BD19", "SAR92")
legend("topright", legend=labels, col=mypal[1:10], pch=15, cex=1.3, bty="n")
box(lwd=2)
dev.off()

##################################################################################
############################# For CANON Percentile Barplot #######################
##################################################################################

name <- c("Ostreococcus", "SAR11", "SAR116", "SAR86", "GII Euryarchaeota", "SAR92", "Flavobacteria", "Roseobacter", "ARCTIC96-BD19", "SAR406")
Result <- data.frame()
#sequ <- seq(0, 1, by=0.05)
sequ <- seq(0.9, 1, by=0.005)
	for(i in 2:length(sequ)) {
	quant_right <- quantile(k$Connectivity, sequ[i])
	quant_left <- quantile(k$Connectivity, sequ[i-1])
	names <- as.character(row.names(k)[which(k$Connectivity > quant_left & k$Connectivity < quant_right )])
	string <- strsplit(names, "_")
	#print(quant)
	list <- c("Ostreo"=0, "SAR11" = 0, "SAR116" = 0, "SAR86"=0, "eury"=0, "SAR92"=0, "Flavo"=0, "Roseo"=0, "Arctic"=0, "SAR406"=0)
		for(j in 1:length(string)) {
			list[string[[j]][1]] = list[string[[j]][1]] + 1
		}
	final <- as.numeric(list)
	Result <- data.frame(rbind(Result, final))
	}
colnames(Result) <- name
row.names(Result) <- 100 * (sequ[2:(length(sequ))] - 0.0025)
library(ggplot2)
library(reshape2)
dat1 <- melt(t(Result))
colnames(dat1) <- c("Taxon", "Percentile", "Proportion")
#cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
jpeg(file="CANON_quantiles_1_100.jpg", quality=100, res=300, height=4, width=8, units="in")
ggplot(dat1, aes(Percentile, Proportion, fill=Taxon), guide=guide_legend(title="size")) + geom_bar(stat="identity", position="fill") + coord_flip() + scale_fill_hue(l=50, c=60) + theme(legend.title = element_text(size=12, face="bold"), legend.text = element_text(size = 10)) + scale_fill_manual(values=mypal)
dev.off()

###################################################################
########################### End ###################################
###################################################################

x <- read.table(file="R_network_pipeline/CANON/blc.all.ortho.clusters.list", sep="\t", quote="", header=TRUE, row.names=1)
sum <- read.table(file="R_network_pipeline/BioLincs/blc.combined.summary", sep="\t", quote="", header=TRUE, row.names=1)


canon <- read.table(file="R_network_pipeline/CANON/Combined_Canon", header=TRUE, row.names=1, sep="\t")
k <- read.table(file="R_network_pipeline/CANON/canon_Network_att", header=TRUE, row.names=1, sep="\t")

canon.cor <- cor(t(canon))
k.order <- k[order(k$Connectivity, decreasing=T),]
pro <- row.names(k.order)[1:8848]


##################################################
allclust <- list()
Result <- c()
for (i in 1:8848) {
	clusters <- colnames(blc.cor)[which(blc.cor[,pro[i]] > 0.75)]
	string <- strsplit(clusters, "Cluster")
	list <- c("Pro_"=0, "SAR11_" = 0, "SAR116_" = 0, "SAR86_"=0, "Roseo_"=0, "SAR406_"=0, "SAR324_" = 0)
		for(j in 1:length(string)) {
			list[string[[j]][1]] = list[string[[j]][1]] + 1
		}
	final <- as.numeric(list)
	Result <- data.frame(rbind(Result, final))
	}
row.names(Result) <- pro
colnames(Result) <- list


blc.cor <- cor(t(blc))
k.order <- k[order(k$Connectivity, decreasing=T),]
pro <- row.names(k.order)[1:100]


#### For adding up weighted networks
canon.cor2 <- (canon.cor)^5

pro <- row.names(k.order)[1:8848]
allclust <- list()
Result <- c()
name <- c("Ostreococcus", "SAR11", "SAR116", "SAR86", "GII Euryarchaeota", "SAR92", "Flavobacteria", "Roseobacter", "ARCTIC96-BD19", "SAR406")
for (i in 1:length(pro)) {
	clusters <- colnames(canon.cor2)
	list <- c("Ostreo_"=0, "SAR11_" = 0, "SAR116_" = 0, "SAR86_"=0, "eury_"=0, "SAR92_"=0, "Flavo_"=0, "Roseo_"=0, "Arctic_"=0, "SAR406_"=0)
		for(j in 1:length(clusters)) {
			id <- strsplit(clusters[j], "Cluster")
			list[id[[1]][1]] = list[id[[1]][1]] + canon.cor2[pro[i],clusters[j]]
		}
	final <- as.numeric(list)
	Result <- data.frame(rbind(Result, final))
	}
row.names(Result) <- pro
colnames(Result) <- name



sum <- apply(Result, 1, sum)
scale <- t(scale(t(Result), scale=sum, center=F))



	allclust[(length(allclust)+1):length(clusters)] <- clusters
	#allclust <- c(allclust, clusters)
	#l <- length(unique(allclust))
	#print(l)
}

############################################### Faster way for weighted network

Result <- c()
list <- c("Ostreo_"=0, "SAR11_" = 0, "SAR116_" = 0, "SAR86_"=0, "eury_"=0, "SAR92_"=0, "Flavo_"=0, "Roseo_"=0, "Arctic_"=0, "SAR406_"=0)
list1 <- c("Ostreo_", "SAR11_", "SAR116_", "SAR86_", "eury_", "SAR92_", "Flavo_", "Roseo_", "Arctic_", "SAR406_")
for (i in 1:length(list)) {
	clusters <- subset(canon.cor, grepl(list1[i], row.names(canon.cor)))
	sum <- apply(clusters, 2, sum)
	list <- c("Ostreo_"=0, "SAR11_" = 0, "SAR116_" = 0, "SAR86_"=0, "eury_"=0, "SAR92_"=0, "Flavo_"=0, "Roseo_"=0, "Arctic_"=0, "SAR406_"=0)
		for(j in 1:8848) {
			id <- strsplit(row.names(canon.cor)[j], "Cluster")
			#if(sum[j] > 50) { 
			list[id[[1]][1]] = list[id[[1]][1]] + sum[j]
			#}
		}
	final <- as.numeric(list)
	Result <- data.frame(rbind(Result, final))
	
	}
colnames(Result) <- list1
row.names(Result) <- list1
write.table(Result, file="output1", sep="\t")


############################################### Faster way for unweighted network

e <- seq(from=0, to=0, length.out=8868)

Result <- c()
list <- c("Ostreo_"=0, "SAR11_" = 0, "SAR116_" = 0, "SAR86_"=0, "eury_"=0, "SAR92_"=0, "Flavo_"=0, "Roseo_"=0, "Arctic_"=0, "SAR406_"=0)
list1 <- c("Ostreo_", "SAR11_", "SAR116_", "SAR86_", "eury_", "SAR92_", "Flavo_", "Roseo_", "Arctic_", "SAR406_")
for (i in 1:length(list)) {
	clusters <- subset(canon.cor2, grepl(list1[i], row.names(canon.cor2)))
	sum <- apply(clusters, 2, sum)
	list <- c("Ostreo_"=0, "SAR11_" = 0, "SAR116_" = 0, "SAR86_"=0, "eury_"=0, "SAR92_"=0, "Flavo_"=0, "Roseo_"=0, "Arctic_"=0, "SAR406_"=0)
		for(j in 1:8848) {
			id <- strsplit(row.names(canon.cor2)[j], "Cluster")
			#if(sum[j] > 50) { 
			list[id[[1]][1]] = list[id[[1]][1]] + sum[j]
			#}
		}
	final <- as.numeric(list)
	Result <- data.frame(rbind(Result, final))
	
	}
colnames(Result) <- list1
row.names(Result) <- list1
write.table(Result, file="output1", sep="\t")

