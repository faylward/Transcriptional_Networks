
##############################################################
############### Plot abundance of microbial groups ##########
##############################################################
combined <- read.table(file="R_network_pipeline/CANON/canon.final.combined.list3", sep="\t", quote="", header=TRUE, row.names=1)
#string <- strsplit(row.names(combined), "Cluster")
#taxa <- c(1:length(row.names(combined)))
#for(i in 1:length(taxa)) {taxa[i] <- string[[i]][1]}
#combined2 <- cbind(combined, taxa)
#combined2[,1:35][which(combined2$taxa == "Ostreo_")]

##########################################################################
#################### NEW PLOT WITH PROPER DISTNACES INBETWEEN ############
##########################################################################

library(RColorBrewer)
mypal <- c("green4", brewer.pal(9, "Set1")); mypal[7] = "gold"

counts1 <- as.numeric(c("259558", "161837", "14930", "153188", "191217", "141794", "134940", "11495", "710590", "21206", "304290", "40380", "10199", "304664", "200826", "49166", "54780", "125220", "99397", "299937", "191111", "224567", "231658", "70186", "125668", "101283", "23796", "136668", "122119", "237315", "163874", "168570", "136397", "92250", "94468"))

t=c(2, 6, 8, 12, 16, 20, 24, 28, 32, 36, 40, 48, 52, 56, 62, 66, 70, 72, 74, 78, 82, 86, 90, 94, 96, 98, 100, 104, 106, 110, 114, 118, 122, 126, 128)
list <- c("Ostreococcus", "SAR11_", "SAR116_", "SAR86", "Euryarchaota", "Roseobacter", "SAR92", "SAR406", "Flavobacteria", "ARCTIC")
list <- c("Ostreococcus", "SAR11_", "SAR116_", "SAR86", "Roseobacter", "SAR406", "Euryarchaota", "SAR92", "Flavobacteria", "ARCTIC")
result <- c()
for(i in 1:length(list)) {counts <- apply(combined[(combined$Organism == list[i]),31:65], 2, sum) ; result <- data.frame(rbind(result, counts))}
row.names(result) <- list
res <- 100*scale(result, scale=counts1, center=F)
colnames(res) <- as.numeric(t)
#dat1 <- melt(t(res)); colnames(dat1) <- c("time", "taxa", "abund")
#ggplot(data=dat1, aes(x=time, y=abund, group=taxa, colour=taxa)) + geom_line(size=1.5)+ scale_colour_manual(values=mypal)

################ Plot with regular R functions
names <- c("Ostreococcus", "SAR11", "SAR116", "SAR86", "GII Euryarchaota", "SAR92", "Flavobacteria", "Roseobacter", "ARCTIC96-BD19", "SAR406")
names <- c("Ostreococcus", "SAR11", "SAR116", "SAR86", "Roseobacter",  "SAR406", "MGII Euryarchaeota", "SAR92", "Flavobacteria", "ARCTIC96-BD19")
res2 <- t(res)
other <- c(1, 2, 3, 4, 5, 6, 8, 9, 10)
jpeg(file="CANON.taxa.abundance.jpg", quality=100, res=600, height=8, width=6, units="in")
par(mfrow=c(5,2), mar=c(2, 3, 2, 1))
for(i in 1:10) {
	color <- mypal[i]
	axis <- as.numeric(t)
	title <- names[i]
	print(title)
	plot(axis, res2[,i], col=color, type="n", ylab=NA, lwd=2, main=title, cex.main=1.5, axes=F)
	rect(xleft = c(-2, 19, 43, 67, 91, 115), xright=c(7, 31, 55, 79, 103, 127), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(100, 100, 100, 100, 100, 100), col='grey90', border=NA)
	par(new=T)
	plot(axis, res2[,i], col=color, type="l", ylab=NA, lwd=2, main=title, cex.main=1.5, axes=F)
		#if(i == 7) { axis(2, at=c(1, 2, 3, 4), labels=c('1', '2', '3', '4'), las=2, cex.axis=1.2, font=2) }
		#else { 
		axis(2, las=2, cex.axis=1.2, font=2)
		# }
	#axis(1, cex.axis=1.2)
	axis(1, side=1,  at=c(12, 36, 60, 84, 108), labels=c("9/12", "9/13", "9/14", "9/15", "9/16"), cex.axis=1.1, font.axis=2)
	#axis(1, at = c(0.1, 5, 10, 15, 20, 25, 30, 35), labels = c(0, 5, 10, 15, 20, 25, 30, 35), cex.axis=1.4)
	box(lty=1)
}
dev.off()

###########################
####### END ###############
###########################


#################### Try plotting with ggplot2
############### Set color palette using RColorBrewer
library(RColorBrewer)
mypal <- c("green4", brewer.pal(9, "Set1")); mypal[7] = "gold"

counts1 <- as.numeric(c("259558", "161837", "14930", "153188", "191217", "141794", "134940", "11495", "710590", "21206", "304290", "40380", "10199", "304664", "200826", "49166", "54780", "125220", "99397", "299937", "191111", "224567", "231658", "70186", "125668", "101283", "23796", "136668", "122119", "237315", "163874", "168570", "136397", "92250", "94468"))

list <- c("Ostreococcus", "SAR11_", "SAR116_", "SAR86", "Euryarchaota", "Roseobacter", "SAR92", "SAR406", "Flavobacteria", "ARCTIC")
list <- c("Ostreococcus", "SAR11_", "SAR116_", "SAR86", "Roseobacter", "SAR406", "Euryarchaeota", "SAR92", "Flavobacteria", "ARCTIC")
result <- c()
for(i in 1:length(list)) {counts <- apply(combined[(combined$Organism == list[i]),31:65], 2, sum) ; result <- data.frame(rbind(result, counts))}
row.names(result) <- list
res <- 100*scale(result, scale=counts1, center=F)
#dat1 <- melt(t(res)); colnames(dat1) <- c("time", "taxa", "abund")
#ggplot(data=dat1, aes(x=time, y=abund, group=taxa, colour=taxa)) + geom_line(size=1.5)+ scale_colour_manual(values=mypal)

################ Plot with regular R functions
names <- c("Ostreococcus", "SAR11", "SAR116", "SAR86", "GII Euryarchaota", "SAR92", "Flavobacteria", "Roseobacter", "ARCTIC96-BD19", "SAR406")
names <- c("Ostreococcus", "SAR11", "SAR116", "SAR86", "Roseobacter",  "SAR406", "GII Euryarchaeota", "SAR92", "Flavobacteria", "ARCTIC96-BD19")
res2 <- t(res)
other <- c(1, 2, 3, 4, 5, 6, 8, 9, 10)
jpeg(file="CANON.taxa.abundance.jpg", quality=100, res=600, height=8, width=6, units="in")
par(mfrow=c(5,2), mar=c(2, 3, 2, 1))
for(i in 1:10) {
	color <- mypal[i]
	axis <- c(1:35)
	title <- names[i]
	print(title)
	plot(axis, res2[,i], col=color, type="n", ylab=NA, lwd=2, main=title, cex.main=1.5, axes=F)
	rect(xleft = c(0, 4.5, 10.5, 15.5, 22.5, 30.5), xright=c(1.5, 7.5, 12.5, 19.5, 26.5, 33.5), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(100, 100, 100, 100, 100, 100), col='grey90', border=NA)
	par(new=T)
	plot(axis, res2[,i], col=color, type="l", ylab=NA, lwd=2, main=title, cex.main=1.5, axes=F)
		#if(i == 7) { axis(2, at=c(1, 2, 3, 4), labels=c('1', '2', '3', '4'), las=2, cex.axis=1.2, font=2) }
		#else { 
		axis(2, las=2, cex.axis=1.2, font=2)
		# }
	axis(1, cex.axis=1.2)
	#axis(1, at = c(0.1, 5, 10, 15, 20, 25, 30, 35), labels = c(0, 5, 10, 15, 20, 25, 30, 35), cex.axis=1.4)
	box(lty=1)
}
dev.off()

#####################################################

names <- c("Ostreococcus", "SAR11", "SAR116", "SAR86", "GII Euryarchaota", "SAR92", "Flavobacteria", "Roseobacter", "ARCTIC96-BD19", "SAR406")
res2 <- t(res)
jpeg(file="canon.taxa.abundance.combined.jpg", quality=100, res=600, height=5, width=7, units="in")
#par(mfrow=c(4,2), mar=c(2, 3, 2, 1))
axis <- c(1:35)
dum <- seq(from=0, to=25, length.out=35)
plot(axis, dum, type="n", ylab=NA, axes=F)
rect(xleft = c(0, 9, 18.5, 27), xright=c(4.5, 12.5, 22.5, 31), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(100, 100, 100, 100, 100, 100), col='grey90', border=NA)
for(i in 1:7) {
	color <- mypal[i]
	#axis <- c(1:30)
	title <- names[i]
	print(title)
	lines(axis, res2[,i], col=color, lwd=2)
	#plot(axis, res2[,i], col=color, type="n", ylab=NA, lwd=2, main=title, cex.main=1.5, axes=F)
	#rect(xleft = c(0, 9, 18.5, 27), xright=c(4.5, 12.5, 22.5, 31), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(100, 100, 100, 100, 100, 100), col='grey90', border=NA)
	#par(new=T)
	#plot(axis, res2[,i], col=color, type="l", ylab=NA, lwd=2, main=title, cex.main=1.5, axes=F)
		#if(i == 7) { axis(2, at=c(1, 2, 3, 4), labels=c('1', '2', '3', '4'), las=2, cex.axis=1.2, font=2) }
		#else { axis(2, las=2, cex.axis=1.2, font=2) }
	axis(2, las=2, cex.axis=1.2, font=2)
	axis(1, cex.axis=1.2)
	#axis(1, at = c(0.1, 5, 10, 15, 20, 25, 30, 35), labels = c(0, 5, 10, 15, 20, 25, 30, 35), cex.axis=1.4)
	box(lty=1)
}
dev.off()


