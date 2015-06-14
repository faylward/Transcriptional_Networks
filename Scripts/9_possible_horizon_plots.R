
#######################################################################################
################################ Plot Modules? #################################
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
modules <- c("steelblue1", "purple", "orange", "springgreen4", "red", "magenta", "steelblue4", "brown", "cyan", "turquoise3", "gold", "tomato1", "violet", "tan4", "turquoise4", "wheat3", "slategray2", "yellow", "tan1", "sienna2", "darkturquoise", "saddlebrown", "orangered", "firebrick1", "orchid4", "royalblue3", "tomato4", "olivedrab1", "olivedrab4", "springgreen3", "skyblue", "chocolate3", "darkcyan", "aquamarine", "coral", "darkgoldenrod", "blueviolet")
eig <- moduleEigengenes(Norm_Final, as.character(Module_Colors[,2]), excludeGrey=TRUE)
data <- eig$eigengenes
data1 <- eig$averageExpr
data <- as.matrix(data1)
timelist <- as.factor(as.character(c('02', '06', '08', 12, 16, 20, 24, '04', '08', 12, 16, 24, '04', '08', 14, 18, 22, 24, '02', '06', 10, 14, 18, 22, 24, '02', '04', '08', 10, 14, 18, 22, '02', '06', '08')))
times <- c(2, 6, 8, 12, 16, 20, 24, 4, 8, 12, 16, 24, 4, 8, 14, 18, 22, 24, 2, 6, 10, 14, 18, 22, 24, 2, 4, 8, 10, 14, 18, 22, 2, 6, 8)
ag <- aggregate(data, list(timelist), mean)
colnames(ag) <- gsub("AE", "", colnames(ag))
colnames(data1) <- gsub("AE", "", colnames(data1))
#ag[,1] <- c(10, 12, 14, 16, 18, '02', 20, 22, 24, '04', '06', '08')
order <- ag[order(ag[,1], decreasing=F),]

###########################################################################
###################### Sparklines #########################################
###########################################################################
	x <- ag[,2:35]
jpeg(file="canon.spark.jpg", quality=100, res=600, height=18, width=3, units="in")
    par(mfrow=c(34,1), mar=c(0.2,2.2,0.2,0), oma=c(0.2,0.2,0.2,0.2))

    for (i in 34:1){ # setup for statement to loops over all elements in a list or vector
        plot(x[,modules[i]], col=modules[i],lwd=2, axes=F,ylab='i',xlab="",main="",type="n")
        rect(xleft = c(0, 9.5), xright=c(3.5, 12), ybottom=c(-3, -3), ytop=c(3, 3), col='grey90', border=NA)
        lines(x[,modules[i]], col=modules[i],lwd=3)
	mtext(i, las=2, side=2, font=2, line=1)
        #plot(x[,modules[i]], col=modules[i],lwd=2, axes=F,ylab='i',xlab="",main="",type="l")
        #axis(2,yaxp=c(min(x[,i]),max(x[,i]),2), cex.axis=1.1,las=1, at=c(round(min(x[,i]),3),round(max(x[,i]),3)))
        #axis(2,yaxp=c(min(x[,i]),max(x[,i]),2),col="white",tcl=0,labels=FALSE)
        #ymin<-min(x[,i]); tmin<-which.min(x[,i]);ymax<-max(x[,i]);tmax<-which.max(x[,i]); # see the code from Jason below for what these do 
        #points(x=c(tmin,tmax),y=c(ymin,ymax),pch=19,col=c("red","blue"),cex=1) # add coloured points at max and min
        lines(seq(0, 0, length.out=12), pch=22, lty=5, lwd=1.5)
        }
        axis(1,pos=c(-5)) # places horizontal axis at the bottom of it all.
dev.off()
#########################################################
###### End Sparklines ###################################
#########################################################

diel <- read.table(file="R_network_pipeline/CANON/All.canon.hra.out", sep="\t", header=TRUE, row.names=1)
diel.transcript <- (diel$Perm.FDR <= 0.1 & diel$Regression.FDR <= 0.1)
diel.2 <- cbind(diel, diel.transcript)

canon <- read.table(file="R_network_pipeline/CANON/canon.combined.summary2", sep="\t", header=TRUE, row.names=1)
final.canon <- merge(diel.2, canon, by="row.names", all=T)
row.names(final.canon) <- final.canon$Row.names

# SAR11 
sar11.ribo <- row.names(final.canon)[which(final.canon$Organism == "SAR11_" & final.canon$diel.transcript == "TRUE" & final.canon$Final_Annote == "Translation")]
sar11.ribos <- apply(final.canon[sar11.ribo, 37:71], 2, mean)
sar11 <- final.canon[sar11.ribo, 37:71]
ag <- aggregate(t(sar11), list(timelist), mean)
row.names(ag) <- as.character(ag$Group.1)
sar11.pag <- apply(ag[,2:dim(ag)[2]], 1, mean)

# SAR116 
sar116.ribo <- row.names(final.canon)[which(final.canon$Organism == "SAR116_" & final.canon$diel.transcript == "TRUE" & final.canon$Final_Annote == "Translation")]
sar116.ribos <- apply(final.canon[sar116.ribo, 37:71], 2, mean)
sar116 <- final.canon[sar116.ribo, 37:71]
ag <- aggregate(t(sar116), list(timelist), mean)
row.names(ag) <- as.character(ag$Group.1)
sar116.pag <- apply(ag[,2:dim(ag)[2]], 1, mean)

# SAR86
sar86.ribo <- row.names(final.canon)[which(final.canon$Organism == "SAR86" & final.canon$diel.transcript == "TRUE" & final.canon$Final_Annote == "Translation")]
sar86.ribos <- apply(final.canon[sar86.ribo, 37:71], 2, mean)
sar86 <- final.canon[sar86.ribo, 37:71]
ag <- aggregate(t(sar86), list(timelist), mean)
row.names(ag) <- as.character(ag$Group.1)
sar86.pag <- apply(ag[,2:dim(ag)[2]], 1, mean)
###

plot(sar11.pag, type="l", col="red")
par(new=T)
plot(sar116.pag, type="l", col="green4")
par(new=T)
plot(sar86.pag, type="l", col="blue")
##########

## Try density plots
# SAR11 
sar11.ribo <- row.names(final.canon)[which(final.canon$Organism == "SAR11_" & final.canon$diel.transcript == "TRUE" & final.canon$Final_Annote == "Translation")]
sar11.ribos <- final.canon[sar11.ribo,]$Peak.Time
sar11.d <- density(sar11.ribos)

# SAR116
sar116.ribo <- row.names(final.canon)[which(final.canon$Organism == "SAR116_" & final.canon$diel.transcript == "TRUE" & final.canon$Final_Annote == "Translation")]
sar116.ribos <- final.canon[sar116.ribo,]$Peak.Time
sar116.d <- density(sar116.ribos)
plot(d, type="l")

# SAR86 
sar86.ribo <- row.names(final.canon)[which(final.canon$Organism == "SAR86" & final.canon$diel.transcript == "TRUE" & final.canon$Final_Annote == "Translation")]
sar86.ribos <- final.canon[sar86.ribo,]$Peak.Time
sar86.d <- density(sar86.ribos)

plot(sar11.d, type="l", col="red")
par(new=T)
plot(sar116.d, type="l", col="green4")
par(new=T)
plot(sar86.d, type="l", col="blue")



#######################################################################################################################

#####
# With dotplots on top
#data <- data.frame(scale(data))
### Set up output for eigengene data; the user will need to adjust the plotting parameters
jpeg(file="canon.eigengenes.3.jpg", quality=100, res=600, height=10, width=7, units="in")
par(mfrow=c(5,2), mar=c(2, 1, 2, 1))
for(i in 1:10) {
#for(i in 1:(dim(data)[2]-1)) {
	color <- modules[i]
	#print(color)
	#ag <- aggregate(order[,color], timelist, mean)
	#axis = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24)
	#axis = c(1, 3, 5, 6, 7, 9, 11, 12, 13, 14, 15, 17, 19, 21, 22, 23)
	axis <- c(1:24)
	#row.names(ag) <- axis
	title <- paste(c("Module ", i, "; ", color, "; ", "n=", mod[modules[i],]), collapse="")
	print(title)
	plot(seq(from=2, to=24, by=2), ag[,modules[i]], col=color, type="l", ylab=NA, lwd=2, main=title, cex.main=1.5, axes=F, ylim=c(-1.1, 1.1))
	rect(xleft = c(0, 19.5), xright=c(6.5, 24), ybottom=c(-1.5, -1.5), ytop=c(1.5, 1.5), col='grey90', border=NA)
	par(new=T)
	plot(seq(from=2, to=24, by=2), ag[,modules[i]], col=color, type="l", ylab=NA, lwd=2, main=title, cex.main=1.5, axes=F, ylim=c(-1.1, 1.1))
	#plot(c(1:12), ag[,modules[i]])
	par(new=T)
	plot(times, data1[,modules[i]], col=color, type='p', ylim=c(-1.1, 1.1), axes=F, pch=16, cex=0.8)
	#axis(1, at = c(0.1, 6, 12, 18, 24), labels = c(0, 6, 12, 18, 24), cex.axis=1.4)
	lines(seq(0, 0, length.out=24), pch=22, lty=5)
	#box(lty=1)
	#par(new=T)
	#plot(data[,i], col=color, type="l", ylab=NA, lwd=1)
	#par(new=T)
	#plot(times, data[,i], col=color, ylab=NA, ylim=c(-0.2, 0.9))
	#par(new=T)
}
dev.off()


# WithOUT dotplots on top
#data <- data.frame(scale(data))
### Set up output for eigengene data; the user will need to adjust the plotting parameters
jpeg(file="canon.eigengenes.3.jpg", quality=100, res=600, height=10, width=7, units="in")
par(mfrow=c(5,2), mar=c(2, 1, 2, 1))
for(i in 1:10) {
#for(i in 1:(dim(data)[2]-1)) {
	color <- modules[i]
	#print(color)
	#ag <- aggregate(order[,color], timelist, mean)
	#axis = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24)
	#axis = c(1, 3, 5, 6, 7, 9, 11, 12, 13, 14, 15, 17, 19, 21, 22, 23)
	axis <- c(1:24)
	#row.names(ag) <- axis
	title <- paste(c("Module ", i, "; ", color, "; ", "n=", mod[modules[i],]), collapse="")
	print(title)
	plot(seq(from=2, to=24, by=2), ag[,modules[i]], col=color, type="l", ylab=NA, lwd=2, main=title, cex.main=1.5, axes=F)
	rect(xleft = c(0, 19.5), xright=c(6.5, 24), ybottom=c(-1.5, -1.5), ytop=c(1.5, 1.5), col='grey90', border=NA)
	par(new=T)
	plot(seq(from=2, to=24, by=2), ag[,modules[i]], col=color, type="l", ylab=NA, lwd=2, main=title, cex.main=1.5, axes=F)
	#plot(c(1:12), ag[,modules[i]])
	#par(new=T)
	#plot(times, data1[,modules[i]], col=color, type='p', ylim=c(-1.1, 1.1), axes=F, pch=16, cex=0.8)
	#axis(1, at = c(0.1, 6, 12, 18, 24), labels = c(0, 6, 12, 18, 24), cex.axis=1.4)
	lines(seq(0, 0, length.out=24), pch=22, lty=5)
	#box(lty=1)
	#par(new=T)
	#plot(data[,i], col=color, type="l", ylab=NA, lwd=1)
	#par(new=T)
	#plot(times, data[,i], col=color, ylab=NA, ylim=c(-0.2, 0.9))
	#par(new=T)
}
dev.off()


fit <- loess(ag[,modules[i]] ~ seq(from=2, to=24, by=2))
plot(c(1:12), ag[,modules[i]])
lines(predict(fit), lwd=2)


############### With curves added

#data <- data.frame(scale(data))
### Set up output for eigengene data; the user will need to adjust the plotting parameters
jpeg(file="canon.eigengenes.3.jpg", quality=100, res=300, height=10, width=7, units="in")
par(mfrow=c(5,2), mar=c(2, 1, 2, 1))
for(i in 1:10) {
#for(i in 1:(dim(data)[2]-1)) {
	color <- modules[i]
	#print(color)
	#ag <- aggregate(order[,color], timelist, mean)
	#axis = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24)
	#axis = c(1, 3, 5, 6, 7, 9, 11, 12, 13, 14, 15, 17, 19, 21, 22, 23)
	axis <- c(1:24)
	#row.names(ag) <- axis
	title <- paste(c("Module ", i, "; ", color, "; ", "n=", mod[modules[i],]), collapse="")
	print(title)
	plot(c(1:12), ag[,modules[i]], col=color, type="p", ylab=NA, lwd=2, main=title, cex.main=1.5, axes=F)
	rect(xleft = c(0, 9.75), xright=c(3.25, 12), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(1, 1, 1, 1, 1, 1), col='grey90', border=NA)
	par(new=T)
	plot(c(1:12), ag[,modules[i]], col=color, type="p", ylab=NA, lwd=2, main=title, cex.main=1.5, axes=F)
	fit <- loess(ag[,modules[i]] ~ seq(from=2, to=24, by=2), span=0.7)
	#plot(c(1:12), ag[,modules[i]])
	axis(1, at = c(0.1, 6, 12, 18, 24), labels = c(0, 6, 12, 18, 24), cex.axis=1.4)
	#lines(seq(0, 0, length.out=24), pch=22, lty=5)
	lines(predict(fit), lwd=2, col=color)
	box(lty=1)
	#par(new=T)
	#plot(data[,i], col=color, type="l", ylab=NA, lwd=1)
	#par(new=T)
	#plot(times, data[,i], col=color, ylab=NA, ylim=c(-0.2, 0.9))
	#par(new=T)
}
dev.off()


############ Horizon plot

require(lattice)
require(latticeExtra)
require(reshape2)
require(quantmod)
 
tckrs <- c("^W0DOW","^GSPC","^RUT","^E1DOW","^P1DOW","^DJUBS")
 
getSymbols(tckrs,from="2011-12-31")
 
#combine prices together
prices <- na.omit(merge(W0DOW[,4],GSPC[,4],RUT[,4],E1DOW[,4],P1DOW[,4],DJUBS[,4]))
#get change since beginning of period
change <- prices/matrix(rep(prices[1,],NROW(prices)),nrow=NROW(prices),ncol=NCOL(prices),byrow=TRUE) -1
colnames(change) <- tckrs
 
#using the example as presented in horizonplot documentation
horizonplot(change,layout=c(1,NCOL(change)),
scale=0.05,
par.settings=theEconomist.theme(box="transparent"),
#if you want y labels in the graph uncomment
# panel = function (x,y,...) {
# panel.horizonplot(x,y,...)
# panel.text(x=x[1],y=0,label=colnames(change)[panel.number()],pos=3)
# },
strip.left = FALSE,
scales = list(y = list(draw = FALSE,relation = "same",alternating=FALSE)),
main="World Indexes Change Since 2011",
xlab=NULL,
ylab = list(rev(colnames(change)), rot = 0, cex = 0.8)) +
#add some separation between graphs with small white band
layer(panel.xblocks(height=0.001,col="white",...))

############################################################################################################
########################## Calculate Taxa Representation Across the Modules Identified #####################
############################################################################################################

############### Set color palette using RColorBrewer
library(RColorBrewer)
mypal <- c("green4", brewer.pal(9, "Set1"))

###################################################
num <- 37
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
row.names(result) <- c("Unassigned", 1:37)
sum <- apply(result, 1, sum)
result[,11] <- sum
order <- result[2:38,][order(result[,11], decreasing=TRUE),]

library(ggplot2)
library(reshape2)
dat1 <- melt(t(result[1:10]))
colnames(dat1) <- c("Taxon", "Module", "Transcripts")
jpeg(file="CANON_Modules.jpg", quality=100, res=500, height=8, width=8, units="in")
ggplot(dat1, aes(Module, Transcripts, fill=Taxon), guide=guide_legend(title="size")) + geom_bar(stat="identity") + coord_flip() + scale_fill_hue(l=50, c=60) + theme(axis.text.y  = element_text(colour="black", size=10, face="bold"), axis.text.x = element_text(colour="black", face="bold"), legend.title = element_text(size=12, face="bold"), legend.text = element_text(size = 10)) + scale_fill_manual(values=mypal)
dev.off()


