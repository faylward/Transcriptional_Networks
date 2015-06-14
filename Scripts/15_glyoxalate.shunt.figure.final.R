
##############################################################################################################
################################## For ICL and Malate synthease across different genomes in CANON ############
##############################################################################################################



canon <- read.table(file="R_network_pipeline/CANON/canon.final.combined.list3", sep="\t", quote="", header=TRUE, row.names=1)

list1 <- c("SAR11_", "SAR116_", "SAR86", "Roseobacter", "SAR406", "ARCTIC")
list <- c("SAR11_", "SAR116_", "SAR86", "Roseobacter", "SAR406", "SAR92", "ARCTIC")

msn <- c("SAR11_Cluster_942", "SAR116_Cluster_499", "SAR86_Cluster_988", "SAR406_Cluster_163", "Roseo_Cluster_517", "SAR92_Cluster_767", "Arctic_Cluster_949")
icl <- c("SAR11_Cluster_1261", "SAR116_Cluster_2430", "SAR86_Cluster_809", "SAR406_Cluster_162", "Roseo_Cluster_1798", "SAR92_Cluster_767", "Arctic_Cluster_952")

icl.result <- c(); msn.result <- c()

for(i in 1:length(icl)) {
	counts <- apply(canon[(canon$Organism == list[i]),31:65], 2, sum)
	spec.icl <- 100*(canon[icl[i], 31:65]/counts)
	icl.result <- data.frame(rbind(icl.result, spec.icl))
}

for(i in 1:length(msn)) {
	counts <- apply(canon[(canon$Organism == list[i]),31:65], 2, sum)
	spec.msn <- 100*(canon[msn[i], 31:65]/counts)
	msn.result <- data.frame(rbind(msn.result, spec.msn))
}

t=c(2, 6, 8, 12, 16, 20, 24, 28, 32, 36, 40, 48, 52, 56, 62, 66, 70, 72, 74, 78, 82, 86, 90, 94, 96, 98, 100, 104, 106, 110, 114, 118, 122, 126, 128)
labels <- c("SAR11", "SAR116", "SAR86", "Roseobacter", "SAR406", "SAR92", "ARCTIC96-BD19")
jpeg(file="canon.glyoxylate.enzymes.jpg", quality=100, res=600, height=10, width=10, units="in")
#colors <- c('red', 'green4', 'blue', 'gold', 'orange', 'darkblue')
par(mfrow=c(4,2), mar=c(2, 3, 2, 1))
for(j in 1:7) {
axis <- as.numeric(t)
plot(axis, as.numeric(icl.result[icl[j],]), type='n', axes=F, main=labels[j], cex.main=1.6, ylim=c(0, max(icl.result[icl[j],], msn.result[msn[j],], na.rm=TRUE)))
rect(xleft = c(-2, 19, 43, 67, 91, 115), xright=c(7, 31, 55, 79, 103, 127), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(100, 100, 100, 100, 100, 100), col='grey90', border=NA)
lines(axis, as.numeric(msn.result[msn[j],]), col="grey65", lwd=2)
lines(axis, as.numeric(icl.result[icl[j],]), col="dodgerblue2", lwd=2)
axis(2, las=2, font=2, cex.axis=1.2)
axis(1, side=1,  at=c(12, 36, 60, 84, 108), labels=c("9/12", "9/13", "9/14", "9/15", "9/16"), cex.axis=1.4, font.axis=2)
box()
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", ylim=c(0, 10), xlim=c(0, 10))
text(4.7, 7.2, labels="*", cex=5, col='dodgerblue2')
text(10.1, 9.9, labels="*", cex=5, col='dodgerblue2')
legend(x=5.5, y=2, cex=1.5, legend=c('Isocitrate lyase', 'Malate synthase'), fill=c('dodgerblue2', 'grey55'))
dev.off()

########################################################################################
######################## BLC glyoxylate cycle enzymes ##################################
########################################################################################

blc <- read.table(file="R_network_pipeline/BioLincs/blc.final.combined.list3", sep="\t", header=TRUE, row.names=1, as.is=TRUE, quote="")
#msn <- subset(blc, grepl('malate_synthase', blc$annote))

msn <- c("SAR11_Cluster_942", "SAR116_Cluster_499", "SAR86_Cluster_988", "SAR406_Cluster_163", "Roseo_Cluster_517")
icl <- c("SAR11_Cluster_1261", "SAR116_Cluster_2430", "SAR86_Cluster_809", "SAR406_Cluster_162")

times <- c(22, 25, 27, 29, 31, 33, 38, 41, 43, 45, 47, 49, 54, 57, 59, 61, 63, 65, 70, 73, 75, 77, 79, 81, 86, 89, 91, 93, 95, 97)
t <- times - 20

## PLOT BLC Enzymes
list <- c("SAR11_", "SAR116_", "SAR86", "SAR406", "Roseo")
icl.result <- c(); msn.result <- c()
for(i in 1:length(icl)) {counts <- apply(blc[(blc$Organism == list[i]),28:57], 2, sum); spec.icl <- 100*(blc[icl[i], 28:57]/counts); icl.result <- data.frame(rbind(icl.result, spec.icl))}
for(i in 1:length(msn)) {counts <- apply(blc[(blc$Organism == list[i]),28:57], 2, sum); spec.msn <- 100*(blc[msn[i], 28:57]/counts); msn.result <- data.frame(rbind(msn.result, spec.msn))}

labels <- c("SAR11", "SAR116", "SAR86", "SAR406", "Roseobacter")
jpeg(file="blc.glyoxylate.enzymes.jpg", quality=100, res=600, height=9, width=10, units="in")
par(mfrow=c(3,2), mar=c(2, 3.5, 2, 1))
for(j in 1:length(labels)) {
axis <- as.numeric(t)
plot(axis, as.numeric(icl.result[icl[j],]), type='n', axes=F, ylab=NA, xlab=NA, main=labels[j], cex.main=1.6, ylim=c(0, max(icl.result[icl[j],], msn.result[msn[j],], na.rm=TRUE)))
rect(xleft = c(0, 23, 47, 71), xright=c(10, 34, 58, 80), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(100, 100, 100, 100, 100, 100), col='grey90', border=NA)
lines(axis, as.numeric(msn.result[msn[j],]), col="grey60", lwd=2)
lines(axis, as.numeric(icl.result[icl[j],]), col="dodgerblue2", lwd=2)
axis(2, las=2, font=2, cex.axis=1.2)
axis(1, side=1,  at=c(16, 40, 64), labels=c("9/8", "9/9", "9/10"), cex.axis=1.4, font.axis=2)
#legend(1, 1, legend=c('Isocitrate lyase', 'Malate synthase'), fill=c('dodgerblue2', 'grey55'))
box()
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", ylim=c(0, 10), xlim=c(0, 10))
text(4.7, 9.8, labels="*", cex=5, col='dodgerblue2')
text(10.1, 9.8, labels="*", cex=5, col='dodgerblue2')
legend(x=5.5, y=2.7, cex=1.5, legend=c('Isocitrate lyase', 'Malate synthase'), fill=c('dodgerblue2', 'grey55'))
dev.off()



