


####### CANON Glyoxylate cycle enzymes analysis
diel <- read.table(file="R_network_pipeline/CANON/All.canon.hra.out", sep="\t", header=TRUE, row.names=1)
diel.transcript <- (diel$Perm.FDR <= 0.1 & diel$Regression.FDR <= 0.1)
diel.2 <- cbind(diel, diel.transcript)

canon <- read.table(file="R_network_pipeline/CANON/canon.combined.summary2", sep="\t", header=TRUE, row.names=1)
final.canon <- merge(diel.2, canon, by="row.names", all=T)
row.names(final.canon) <- final.canon$Row.names

#####################################################################################
################# Plot SAR11 SucD, ICL, and PR for CANON  ###########################
#####################################################################################
ICL <- c("SAR11_Cluster_1331")
PR <- c("SAR11_Cluster_790")
SUC <- c("SAR11_Cluster_12")
ACON <- c("SAR11_Cluster_640")

enz <- c("SAR11_Cluster_1331", "SAR11_Cluster_790", "SAR11_Cluster_12")
#suc <- c(); pr <- c(); icl <- c() 
#for(i in 1:length(enz)) {
# CANON
	counts <- apply(canon[(canon$Organism == 'SAR11_'),25:59], 2, sum)
	icl.c <- as.numeric(100*(canon[ICL, 25:59]/counts))
	pr.c <- as.numeric(100*(canon[PR, 25:59]/(counts*10)))
	suc.c <- as.numeric(100*(canon[SUC, 25:59]/counts))
	acon.c <- as.numeric(100*(canon[ACON, 25:59]/counts))
	final <- cbind(icl.c, pr.c, acon.c)

##############################################################
######### Plot whole time series in spark lines ##############
##############################################################
	jpeg(file="CANON.TCA.PR.jpg", quality=100, res=600, height=10, width=7.5, units="in")
	par(mfrow=c(4,1), mar=c(0.2,2.2,1,0), oma=c(3,0.2,0.2,0.2))

	plot(pr.c, col="orange", lwd=2.5, ylab=NA, axes=F, type='n')
	rect(xleft = c(0, 4.5, 10.5, 15.5, 22.5, 30.5), xright=c(1.5, 7.5, 12.5, 19.5, 26.5, 33.5), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(100, 100, 100, 100, 100, 100), col='grey90', border=NA)
	lines(pr.c, col="orange", lwd=2.5)
	title('PR', line= -1.5, cex.main=2)

	plot(icl.c, col="dodgerblue2", lwd=2.5, ylab=NA, axes=F, type='n')
	rect(xleft = c(0, 4.5, 10.5, 15.5, 22.5, 30.5), xright=c(1.5, 7.5, 12.5, 19.5, 26.5, 33.5), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(100, 100, 100, 100, 100, 100), col='grey90', border=NA)
	lines(icl.c, col="dodgerblue2", lwd=2.5)
	title('ICL', line= -1.5, cex.main=2)

	plot(acon.c, col="green4", lwd=2.5, ylab=NA, axes=F, type='n')
	rect(xleft = c(0, 4.5, 10.5, 15.5, 22.5, 30.5), xright=c(1.5, 7.5, 12.5, 19.5, 26.5, 33.5), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(100, 100, 100, 100, 100, 100), col='grey90', border=NA)
	lines(acon.c, col="green4", lwd=2.5)
	title('AcnA', line= -1.5, cex.main=2)

	plot(suc.c, col="blue", lwd=2.5, ylab=NA, axes=F, type='n')
	rect(xleft = c(0, 4.5, 10.5, 15.5, 22.5, 30.5), xright=c(1.5, 7.5, 12.5, 19.5, 26.5, 33.5), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(100, 100, 100, 100, 100, 100), col='grey90', border=NA)
	lines(suc.c, col="blue", lwd=2.5)
	title('2OD', line= -1.5, cex.main=2)
	axis(1, font=2, cex.axis=1.5)

	dev.off()
##############################################

###########################################################
######################## BLC  SAR11 #######################
###########################################################
ICL <- c("SAR11_Cluster_1331")
PR <- c("SAR11_Cluster_790")
SUC <- c("SAR11_Cluster_13")
SUC <- c("SAR11_Cluster_12")
ACON <- c("SAR11_Cluster_640")

diel.blc <- read.table(file="R_network_pipeline/BioLincs/All.blc.hra", sep="\t", header=TRUE, row.names=1)
diel.t.blc <- (diel.blc$Perm.FDR <= 0.1 & diel.blc$Regression.FDR <= 0.1)
diel.2.blc <- cbind(diel.blc, diel.t.blc)

blc <- read.table(file="R_network_pipeline/BioLincs/blc.combined.summary", sep="\t", header=TRUE, row.names=1, as.is=TRUE, quote="")
final.blc <- merge(diel.2.blc, blc, by="row.names", all=T)
row.names(final.blc) <- final.blc$Row.names

	counts <- apply(final.blc[(final.blc$Organism == 'SAR11_'),30:59], 2, sum)
	icl.c <- as.numeric(100*(final.blc[ICL, 30:59]/counts))
	pr.c <- as.numeric(100*(final.blc[PR, 30:59]/(counts*10)))
	suc.c <- as.numeric(100*(final.blc[SUC, 30:59]/counts))
	acon.c <- as.numeric(100*(final.blc[ACON, 30:59]/counts))
	final.b <- cbind(icl.c, pr.c, acon.c)

##############################################################
######### Plot whole time series in spark lines ##############
##############################################################
	jpeg(file="BLC.TCA.PR.jpg", quality=100, res=600, height=10, width=7.5, units="in")
	par(mfrow=c(4,1), mar=c(0.2,2.2,1,0), oma=c(3,0.2,0.2,0.2))

	plot(pr.c, col="orange", lwd=2.5, ylab=NA, axes=F, type='n')
	rect(xleft = c(0, 9, 18.5, 27), xright=c(4.5, 12.5, 22.5, 31), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(100, 100, 100, 100, 100, 100), col='grey90', border=NA)
	lines(pr.c, col="orange", lwd=2.5)
	title('PR', line= -1.5, cex.main=2)

	plot(icl.c, col="dodgerblue2", lwd=2.5, ylab=NA, axes=F, type='n')
	rect(xleft = c(0, 9, 18.5, 27), xright=c(4.5, 12.5, 22.5, 31), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(100, 100, 100, 100, 100, 100), col='grey90', border=NA)
	lines(icl.c, col="dodgerblue2", lwd=2.5)
	title('ICL', line= -1.5, cex.main=2)

	plot(acon.c, col="green4", lwd=2.5, ylab=NA, axes=F, type='n')
	rect(xleft = c(0, 9, 18.5, 27), xright=c(4.5, 12.5, 22.5, 31), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(100, 100, 100, 100, 100, 100), col='grey90', border=NA)
	lines(acon.c, col="green4", lwd=2.5)
	title('AcnA', line= -1.5, cex.main=2)

	plot(suc.c, col="blue", lwd=2.5, ylab=NA, axes=F, type='n')
	rect(xleft = c(0, 9, 18.5, 27), xright=c(4.5, 12.5, 22.5, 31), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(100, 100, 100, 100, 100, 100), col='grey90', border=NA)
	lines(suc.c, col="blue", lwd=2.5)
	title('2OD', line= -1.5, cex.main=2)
	axis(1, font=2, cex.axis=1.5)

	dev.off()
##############################################

###########################################################
######################## BLC  SAR116 #######################
###########################################################
ICL <- c("SAR116_Cluster_2430")
PR <- c("SAR116_Cluster_989")
SUC <- c("SAR116_Cluster_2272")
SUC <- c("SAR116_Cluster_1496")
ACON <- c("SAR116_Cluster_928")

diel.blc <- read.table(file="R_network_pipeline/BioLincs/All.blc.hra", sep="\t", header=TRUE, row.names=1)
diel.t.blc <- (diel.blc$Perm.FDR <= 0.1 & diel.blc$Regression.FDR <= 0.1)
diel.2.blc <- cbind(diel.blc, diel.t.blc)

blc <- read.table(file="R_network_pipeline/BioLincs/blc.combined.summary", sep="\t", header=TRUE, row.names=1, as.is=TRUE, quote="")
final.blc <- merge(diel.2.blc, blc, by="row.names", all=T)
row.names(final.blc) <- final.blc$Row.names

	counts <- apply(final.blc[(final.blc$Organism == 'SAR116_'),30:59], 2, sum)
	icl.c <- as.numeric(100*(final.blc[ICL, 30:59]/counts))
	pr.c <- as.numeric(100*(final.blc[PR, 30:59]/(counts*10)))
	suc.c <- as.numeric(100*(final.blc[SUC, 30:59]/counts))
	acon.c <- as.numeric(100*(final.blc[ACON, 30:59]/counts))
	final.b <- cbind(icl.c, pr.c, acon.c)

##############################################################
######### Plot whole time series in spark lines ##############
##############################################################
	jpeg(file="BLC.TCA.PR.SAR116.jpg", quality=100, res=600, height=10, width=7.5, units="in")
	par(mfrow=c(4,1), mar=c(0.2,2.2,1,0), oma=c(3,0.2,0.2,0.2))

	plot(pr.c, col="orange", lwd=2.5, ylab=NA, axes=F, type='n')
	rect(xleft = c(0, 9, 18.5, 27), xright=c(4.5, 12.5, 22.5, 31), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(100, 100, 100, 100, 100, 100), col='grey90', border=NA)
	lines(pr.c, col="orange", lwd=2.5)
	title('PR', line= -1.5, cex.main=2)

	plot(icl.c, col="dodgerblue2", lwd=2.5, ylab=NA, axes=F, type='n')
	rect(xleft = c(0, 9, 18.5, 27), xright=c(4.5, 12.5, 22.5, 31), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(100, 100, 100, 100, 100, 100), col='grey90', border=NA)
	lines(icl.c, col="dodgerblue2", lwd=2.5)
	title('ICL', line= -1.5, cex.main=2)

	plot(acon.c, col="green4", lwd=2.5, ylab=NA, axes=F, type='n')
	rect(xleft = c(0, 9, 18.5, 27), xright=c(4.5, 12.5, 22.5, 31), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(100, 100, 100, 100, 100, 100), col='grey90', border=NA)
	lines(acon.c, col="green4", lwd=2.5)
	title('AcnA', line= -1.5, cex.main=2)

	plot(suc.c, col="blue", lwd=2.5, ylab=NA, axes=F, type='n')
	rect(xleft = c(0, 9, 18.5, 27), xright=c(4.5, 12.5, 22.5, 31), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(100, 100, 100, 100, 100, 100), col='grey90', border=NA)
	lines(suc.c, col="blue", lwd=2.5)
	title('2OD', line= -1.5, cex.main=2)
	axis(1, font=2, cex.axis=1.5)

	dev.off()

################################################

#####################################################################################
################# Plot SAR11 SucD, ICL, and PR for CANON  ###########################
#####################################################################################
ICL <- c("SAR116_Cluster_2430")
PR <- c("SAR116_Cluster_989")
SUC <- c("SAR116_Cluster_1496")
ACON <- c("SAR116_Cluster_928")

enz <- c("SAR11_Cluster_1331", "SAR11_Cluster_790", "SAR11_Cluster_12")
#suc <- c(); pr <- c(); icl <- c() 
#for(i in 1:length(enz)) {
# CANON
	counts <- apply(canon[(canon$Organism == 'SAR11_'),25:59], 2, sum)
	icl.c <- as.numeric(100*(canon[ICL, 25:59]/counts))
	pr.c <- as.numeric(100*(canon[PR, 25:59]/(counts*10)))
	suc.c <- as.numeric(100*(canon[SUC, 25:59]/counts))
	acon.c <- as.numeric(100*(canon[ACON, 25:59]/counts))
	final <- cbind(icl.c, pr.c, acon.c)

##############################################################
######### Plot whole time series in spark lines ##############
##############################################################
	jpeg(file="CANON.TCA.PR.SAR116.jpg", quality=100, res=600, height=10, width=7.5, units="in")
	par(mfrow=c(4,1), mar=c(0.2,2.2,1,0), oma=c(3,0.2,0.2,0.2))

	plot(pr.c, col="orange", lwd=2.5, ylab=NA, axes=F, type='n')
	rect(xleft = c(0, 4.5, 10.5, 15.5, 22.5, 30.5), xright=c(1.5, 7.5, 12.5, 19.5, 26.5, 33.5), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(100, 100, 100, 100, 100, 100), col='grey90', border=NA)
	lines(pr.c, col="orange", lwd=2.5)
	title('PR', line= -1.5, cex.main=2)

	plot(icl.c, col="dodgerblue2", lwd=2.5, ylab=NA, axes=F, type='n')
	rect(xleft = c(0, 4.5, 10.5, 15.5, 22.5, 30.5), xright=c(1.5, 7.5, 12.5, 19.5, 26.5, 33.5), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(100, 100, 100, 100, 100, 100), col='grey90', border=NA)
	lines(icl.c, col="dodgerblue2", lwd=2.5)
	title('ICL', line= -1.5, cex.main=2)

	plot(acon.c, col="green4", lwd=2.5, ylab=NA, axes=F, type='n')
	rect(xleft = c(0, 4.5, 10.5, 15.5, 22.5, 30.5), xright=c(1.5, 7.5, 12.5, 19.5, 26.5, 33.5), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(100, 100, 100, 100, 100, 100), col='grey90', border=NA)
	lines(acon.c, col="green4", lwd=2.5)
	title('AcnA', line= -1.5, cex.main=2)

	plot(suc.c, col="blue", lwd=2.5, ylab=NA, axes=F, type='n')
	rect(xleft = c(0, 4.5, 10.5, 15.5, 22.5, 30.5), xright=c(1.5, 7.5, 12.5, 19.5, 26.5, 33.5), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(100, 100, 100, 100, 100, 100), col='grey90', border=NA)
	lines(suc.c, col="blue", lwd=2.5)
	title('2OD', line= -1.5, cex.main=2)
	axis(1, font=2, cex.axis=1.5)

	dev.off()
##############################################


#### Plot whole time series

	plot(icl.c, type='n', axes=F, main='TCA + PR', ylim=c(0, max(as.numeric(c(icl.c, pr.c, suc.c)), na.rm=TRUE)), axes=F)
	rect(xleft = c(0, 9, 18.5, 27), xright=c(4.5, 12.5, 22.5, 31), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(100, 100, 100, 100, 100, 100), col='grey90', border=NA)
	par(new=T)
	plot(icl.c, col="dodgerblue2", lwd=2, ylab=NA, axes=F, type='l')
	par(new=T)
	plot(pr.c, col="orange", lwd=2, ylab=NA, axes=F, type='l')
	par(new=T)
	plot(acon.c, col="green4", lwd=2, ylab=NA, axes=F, type='l')
	axis(2, las=2, font=2, cex.axis=1.2)
	axis(1, font=2, cex.axis=1.2)

##### Plot aggregated time seris (24 hrs)
timelist <- as.factor(as.character(c(22, '01', '03', '05', '07', '09', 12, 17, 19, 21, 23, '01', '06', '09', 11, 13, 15, 17, 22, '01', '03', '05', '07', '09', 14, 17, 19, 21, 23, '01')))
times <- c(22, 1, 3, 5, 7, 9, 12, 17, 19, 21, 23, 1, 6, 9, 11, 13, 15, 17, 22, 1, 3, 5, 7, 9, 14, 17, 19, 21, 23, 1)
ag <- aggregate(final.b, list(timelist), mean)

	plot(as.numeric(ag[,1]), ag[,2], type='n', axes=F, main='TCA + PR', ylim=c(0, max(as.numeric(c(ag[,2], ag[,3],ag[,4])), na.rm=TRUE)))
	rect(xleft = c(0, 12), xright=c(5, 16), ybottom=c(-1, -1), ytop=c(100, 100), col='grey90', border=NA)
	lines(as.numeric(ag[,2]), col="dodgerblue2", lwd=2)
	lines(as.numeric(ag[,3]), col="orange", lwd=2)
	points(pr.c, times, col='orange')
	lines(as.numeric(ag[,4]), col="green4", lwd=2)
	axis(2, las=2, font=2, cex.axis=1.2)
	axis(1, font=2, cex.axis=1.2)


################################################
################# End ##########################
################################################













########### Rest of citrate cycle

suc <- c("Arctic_Cluster_1182", "Arctic_Cluster_1184", "Roseo_Cluster_1374", "Roseo_Cluster_897", "SAR116_Cluster_1300", "SAR11_Cluster_12",    "SAR11_Cluster_871") 
list <- c("ARCTIC", "ARCTIC", "Roseobacter", "Roseobacter", "SAR116_", "SAR11_", "SAR11_")
suc.result <- c()
for(i in 1:length(suc)) {
	counts <- apply(canon[(canon$Organism == list[i]),25:59], 2, sum)
	spec.suc <- 100*(canon[suc[i], 25:59]/counts)
	suc.result <- data.frame(rbind(suc.result, spec.suc))
}

par(mfrow=c(4,2), mar=c(2, 3, 2, 1))
for(j in 1:length(suc)) {
plot(as.numeric(suc.result[suc[j],]), type='n', axes=F, main=suc[j], cex.main=1.6, ylim=c(0, max(suc.result[suc[j],], na.rm=TRUE)))
rect(xleft = c(0, 4.5, 10.5, 15.5, 22.5, 30.5), xright=c(1.5, 7.5, 12.5, 19.5, 26.5, 33.5), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(100, 100, 100, 100, 100, 100), col='grey90', border=NA)
lines(as.numeric(suc.result[suc[j],]), col="dodgerblue2", lwd=2)
axis(2, las=2, font=2, cex.axis=1.2)
axis(1, font=2, cex.axis=1.2)
box()
}


## SAR11 Specific: Comparison of Isocitrate lyase and succinate dehydrogenase 
com <- c('SAR11_Cluster_871', 'SAR11_Cluster_1331')
counts <- apply(canon[(canon$Organism == 'SAR11_'),25:59], 2, sum)
suc <- 100*(canon[com[1], 25:59]/counts)
icl <- 100*(canon[com[2], 25:59]/counts)
plot(as.numeric(suc), type='n', axes=F, main='SAR11 SucA_ICL', cex.main=1.6, ylim=c(0, max(as.numeric(c(suc, icl)), na.rm=TRUE)))
rect(xleft = c(0, 4.5, 10.5, 15.5, 22.5, 30.5), xright=c(1.5, 7.5, 12.5, 19.5, 26.5, 33.5), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(100, 100, 100, 100, 100, 100), col='grey90', border=NA)
lines(as.numeric(suc), col="dodgerblue2", lwd=2)
lines(as.numeric(icl), col="red", lwd=2)
axis(2, las=2, font=2, cex.axis=1.2)
axis(1, font=2, cex.axis=1.2)
box()

## SAR116 Specific: Comparison of Isocitrate lyase and succinate dehydrogenase 
com <- c('SAR116_Cluster_1300', 'SAR116_Cluster_2430')
counts <- apply(canon[(canon$Organism == 'SAR116_'),25:59], 2, sum)
suc <- 100*(canon[com[1], 25:59]/counts)
icl <- 100*(canon[com[2], 25:59]/counts)
plot(as.numeric(suc), type='n', axes=F, main='SAR11 SucA_ICL', cex.main=1.6, ylim=c(0, max(as.numeric(c(suc, icl)), na.rm=TRUE)))
rect(xleft = c(0, 4.5, 10.5, 15.5, 22.5, 30.5), xright=c(1.5, 7.5, 12.5, 19.5, 26.5, 33.5), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(100, 100, 100, 100, 100, 100), col='grey90', border=NA)
lines(as.numeric(suc), col="dodgerblue2", lwd=2)
lines(as.numeric(icl), col="red", lwd=2)
axis(2, las=2, font=2, cex.axis=1.2)
axis(1, font=2, cex.axis=1.2)
box()
###########################################################################


###########################################################################



row.names(canon)[which($canon$Organism == "SAR116_" & grepl('malate synthase', canon$Annotation)
msn <- subset(canon, grepl('malate synthase', canon$Annotation))


msn <- c("SAR11_Cluster_267", "SAR116_Cluster_499", "SAR86_Cluster_988", "Roseo_Cluster_517", "SAR406_Cluster_163", "SAR92_Cluster_767", "Arctic_Cluster_949")
icl <- c("SAR11_Cluster_1331", "SAR116_Cluster_2430", "SAR86_Cluster_809", "Roseo_Cluster_1798", "SAR406_Cluster_162", "SAR92_Cluster_767", "Arctic_Cluster_952")

msn <- c("SAR11_Cluster_267", "SAR116_Cluster_499", "SAR86_Cluster_988", "Roseo_Cluster_517", "SAR406_Cluster_163", "Arctic_Cluster_949")
icl <- c("SAR11_Cluster_1331", "SAR116_Cluster_2430", "SAR86_Cluster_809", "Roseo_Cluster_1798", "SAR406_Cluster_162", "SAR92_Cluster_767")

canon["Arctic_Cluster_949", ]

###############################################################################################
##### Perform Mann-Whitney U test to see if ICL expression is enriched in daytime samples #####
###############################################################################################
times <- c(2, 6, 8, 12, 16, 20, 24, 4, 8, 12, 16, 24, 4, 8, 14, 18, 22, 24, 2, 6, 10, 14, 18, 22, 24, 2, 4, 8, 10, 14, 18, 22, 2, 6, 8)
day <- times > 6 & times < 20
night <- times > 18 | times < 7
canon.perc <- 100*(canon[,25:59]/apply(canon[,25:59], 2, sum))

list <- c()
for(j in 1:length(icl)) {
	day.val <- as.numeric(canon.perc[icl[j],][day])
	night.val <- as.numeric(canon.perc[icl[j],][night])
	wc <- wilcox.test(day.val, night.val, alternative='g')
	list[j] <- wc$p.value
}
##############################################################################################

acon <- c('SAR11_Cluster_640', 'Roseo_Cluster_1403', 'SAR406_Cluster_333', 'SAR116_Cluster_680', 'SAR86_Cluster_1085' 'SAR92_Cluster_841', 'Arctic_Cluster_485')
acon <- c('SAR11_Cluster_640', 'SAR116_Cluster_680', 'SAR86_Cluster_1085', 'Roseo_Cluster_1403', 'SAR406_Cluster_333', 'Arctic_Cluster_485')

'SAR11_Cluster_497'
counts <- apply(canon[(canon$Organism == list[i]),25:59], 2, sum)
spec.icl <- 100*(canon['SAR11_Cluster_1105', 25:59]/counts)


result <- c()
list <- c("SAR11_", "SAR116_", "SAR86", "Roseobacter", "SAR406", "SAR92")
for(i in 1:length(icl)) {
	counts <- apply(canon[(canon$Organism == list[i]),25:59], 2, sum)
	spec.icl <- 100*(canon[icl[i], 25:59]/counts)
	spec.acon <- 100*(canon[acon[i], 25:59]/counts)
	pear <- as.numeric(cor(as.numeric(spec.icl), as.numeric(spec.acon), method='pearson'))
	result[i] <- pear
}


	plot(icl.c, type='n', axes=F, main='TCA + PR', ylim=c(0, max(as.numeric(c(icl.c, pr.c, suc.c)), na.rm=TRUE)))
	rect(xleft = c(0, 4.5, 10.5, 15.5, 22.5, 30.5), xright=c(1.5, 7.5, 12.5, 19.5, 26.5, 33.5), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(100, 100, 100, 100, 100, 100), col='grey90', border=NA)
	par(new=T)
	plot(icl.c, col="dodgerblue2", lwd=2.5, ylab=NA, axes=F, type='l')
	par(new=T)
	plot(pr.c, col="orange", lwd=2.5, ylab=NA, axes=F, type='l')
	par(new=T)
	plot(acon.c, col="green4", lwd=2.5, ylab=NA, axes=F, type='l')
	axis(2, las=2, font=2, cex.axis=1.2)
	axis(1, font=2, cex.axis=1.2)

##### Plot aggregated time seris (24 hrs)
timelist <- as.factor(as.character(c('02', '06', '08', 12, 16, 20, 24, '04', '08', 12, 16, 24, '04', '08', 14, 18, 22, 24, '02', '06', 10, 14, 18, 22, 24, '02', '04', '08', 10, 14, 18, 22, '02', '06', '08')))
times <- c(2, 6, 8, 12, 16, 20, 24, 4, 8, 12, 16, 24, 4, 8, 14, 18, 22, 24, 2, 6, 10, 14, 18, 22, 24, 2, 4, 8, 10, 14, 18, 22, 2, 6, 8)
ag <- aggregate(final, list(timelist), mean)

	plot(as.numeric(ag[,2]), type='n', axes=F, main='TCA + PR', ylim=c(0, max(as.numeric(c(ag[,2])), na.rm=TRUE)))
	rect(xleft = c(0, 9), xright=c(3, 12), ybottom=c(-1, -1), ytop=c(100, 100), col='grey90', border=NA)
	lines(as.numeric(ag[,2]), col="dodgerblue2", lwd=2)
	lines(as.numeric(ag[,3]), col="orange", lwd=2)
	lines(as.numeric(ag[,4]), col="green4", lwd=2)
	axis(2, las=2, font=2, cex.axis=1.2)
	axis(1, font=2, cex.axis=1.2)


##############################################################################################################
################################## For ICL and Malate synthease across different genomes in CANON ############
##############################################################################################################






list1 <- c("SAR11_", "SAR116_", "SAR86", "Roseobacter", "SAR406", "ARCTIC")
list <- c("SAR11_", "SAR116_", "SAR86", "Roseobacter", "SAR406", "SAR92")
icl.result <- c(); msn.result <- c()
for(i in 1:length(icl)) {counts <- apply(canon[(canon$Organism == list[i]),25:59], 2, sum); spec.icl <- 100*(canon[icl[i], 25:59]/counts); icl.result <- data.frame(rbind(icl.result, spec.icl))}
for(i in 1:length(msn)) {counts <- apply(canon[(canon$Organism == list1[i]),25:59], 2, sum); spec.msn <- 100*(canon[msn[i], 25:59]/counts); msn.result <- data.frame(rbind(msn.result, spec.msn))}

msn <- c("SAR11_Cluster_267", "SAR116_Cluster_499", "SAR86_Cluster_988", "SAR406_Cluster_163", "Roseo_Cluster_517", "SAR92_Cluster_767", "Arctic_Cluster_949")
icl <- c("SAR11_Cluster_1331", "SAR116_Cluster_2430", "SAR86_Cluster_809", "SAR406_Cluster_162", "Roseo_Cluster_1798", "SAR92_Cluster_767", "Arctic_Cluster_952")
labels <- c("SAR11", "SAR116", "SAR86", "Roseobacter", "SAR406", "SAR92", "ARCTIC96-BD19")
jpeg(file="canon.glyoxylate.enzymes.jpg", quality=100, res=600, height=10, width=10, units="in")
#colors <- c('red', 'green4', 'blue', 'gold', 'orange', 'darkblue')
par(mfrow=c(4,2), mar=c(2, 3, 2, 1))
for(j in 1:length(labels)) {
plot(as.numeric(icl.result[icl[j],]), type='n', axes=F, col=colors[i], main=labels[j], cex.main=1.6, ylim=c(0, max(icl.result[icl[j],], msn.result[msn[j],], na.rm=TRUE)))
rect(xleft = c(0, 4.5, 10.5, 15.5, 22.5, 30.5), xright=c(1.5, 7.5, 12.5, 19.5, 26.5, 33.5), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(100, 100, 100, 100, 100, 100), col='grey90', border=NA)
lines(as.numeric(msn.result[msn[j],]), col="grey65", lwd=2)
lines(as.numeric(icl.result[icl[j],]), col="dodgerblue2", lwd=2)
axis(2, las=2, font=2, cex.axis=1.2)
axis(1, font=2, cex.axis=1.2)
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

diel.blc <- read.table(file="R_network_pipeline/BioLincs/All.blc.hra", sep="\t", header=TRUE, row.names=1)
diel.t.blc <- (diel.blc$Perm.FDR <= 0.1 & diel.blc$Regression.FDR <= 0.1)
diel.2.blc <- cbind(diel.blc, diel.t.blc)

blc <- read.table(file="R_network_pipeline/BioLincs/blc.combined.summary", sep="\t", header=TRUE, row.names=1, as.is=TRUE, quote="")
final.blc <- merge(diel.2.blc, blc, by="row.names", all=T)
row.names(final.blc) <- final.blc$Row.names

diel.c <- row.names(final.blc)[which(final.blc$diel.t.blc == TRUE & final.blc$Final_Annote == "Citrate_Cycle")]

#msn <- subset(blc, grepl('malate_synthase', blc$annote))

msn <- c("SAR11_Cluster_267", "SAR116_Cluster_499", "SAR86_Cluster_988", "SAR406_Cluster_163", "Roseo_Cluster_517")
icl <- c("SAR11_Cluster_1331", "SAR116_Cluster_2430", "SAR86_Cluster_809", "SAR406_Cluster_162")

###############################################################################################
##### Perform Mann-Whitney U test to see if ICL expression is enriched in daytime samples #####
###############################################################################################
times <- c(22, 1, 3, 5, 7, 9, 12, 17, 19, 21, 23, 1, 6, 9, 11, 13, 15, 17, 22, 1, 3, 5, 7, 9, 14, 17, 19, 21, 23, 1)
day <- times > 6 & times < 19
night <- times > 18 | times < 7
blc.perc <- 100*(blc[,18:47]/apply(blc[,18:47], 2, sum))

list <- c()
for(j in 1:length(icl)) {
	day.val <- as.numeric(blc.perc[icl[j],][day])
	night.val <- as.numeric(blc.perc[icl[j],][night])
	wc <- wilcox.test(day.val, night.val, alternative='g')
	list[j] <- wc$p.value
}

plot(day.val, seq(from=1, to=1, length=length(day.val)), ylim=c(0, 3), xlim=c(0, max(day.val, night.val)))
par(new=T)
plot(night.val, seq(from=2, to=2, length=length(night.val)), ylim=c(0, 3), xlim=c(0, max(day.val, night.val)))
##############################################################################################

## PLOT BLC Enzymes
list <- c("SAR11_", "SAR116_", "SAR86", "SAR406", "Roseo")
icl.result <- c(); msn.result <- c()
for(i in 1:length(icl)) {counts <- apply(blc[(blc$Organism == list[i]),18:47], 2, sum); spec.icl <- 100*(blc[icl[i], 18:47]/counts); icl.result <- data.frame(rbind(icl.result, spec.icl))}
for(i in 1:length(msn)) {counts <- apply(blc[(blc$Organism == list[i]),18:47], 2, sum); spec.msn <- 100*(blc[msn[i], 18:47]/counts); msn.result <- data.frame(rbind(msn.result, spec.msn))}

labels <- c("SAR11", "SAR116", "SAR86", "SAR406", "Roseobacter")
jpeg(file="blc.glyoxylate.enzymes.jpg", quality=100, res=600, height=9, width=10, units="in")
par(mfrowc(3,2), mar=c(2, 3.5, 2, 1))
for(j in 1:length(labels)) {
plot(as.numeric(icl.result[icl[j],]), type='n', axes=F, ylab=NA, xlab=NA, col=colors[i], main=labels[j], cex.main=1.6, ylim=c(0, max(icl.result[icl[j],], msn.result[msn[j],], na.rm=TRUE)))
rect(xleft = c(0, 9, 18.5, 27), xright=c(4.5, 12.5, 22.5, 31), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(100, 100, 100, 100, 100, 100), col='grey90', border=NA)
lines(as.numeric(msn.result[msn[j],]), col="grey60", lwd=2)
lines(as.numeric(icl.result[icl[j],]), col="dodgerblue2", lwd=2)
axis(2, las=2, font=2, cex.axis=1.2)
axis(1, font=2, cex.axis=1.2)
#legend(1, 1, legend=c('Isocitrate lyase', 'Malate synthase'), fill=c('dodgerblue2', 'grey55'))
box()
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n", ylim=c(0, 10), xlim=c(0, 10))
text(4.7, 9.8, labels="*", cex=5, col='dodgerblue2')
text(10.1, 9.8, labels="*", cex=5, col='dodgerblue2')
legend(x=5.5, y=2.7, cex=1.5, legend=c('Isocitrate lyase', 'Malate synthase'), fill=c('dodgerblue2', 'grey55'))
dev.off()


# Plot individual enzymes ICL and succinate dehydrogenase
## SAR11 Specific: Comparison of Isocitrate lyase and succinate dehydrogenase 
com <- c('SAR11_Cluster_871', 'SAR11_Cluster_1331')
counts <- apply(final.blc[(final.blc$Organism == 'SAR11_'),30:59], 2, sum)
suc <- 100*(final.blc[com[1], 30:59]/counts)
icl <- 100*(final.blc[com[2], 30:59]/counts)
plot(as.numeric(suc), type='n', axes=F, main='SAR11 SucA_ICL', cex.main=1.6, ylim=c(0, max(as.numeric(c(suc, icl)), na.rm=TRUE)))
rect(xleft = c(0, 9, 18.5, 27), xright=c(4.5, 12.5, 22.5, 31), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(1, 1, 1, 1, 1, 1), col='grey90', border=NA)
lines(as.numeric(suc), col="dodgerblue2", lwd=2)
lines(as.numeric(icl), col="red", lwd=2)
axis(2, las=2, font=2, cex.axis=1.2)
axis(1, font=2, cex.axis=1.2)
box()

#########

## SAR116 Specific: Comparison of Isocitrate lyase and succinate dehydrogenase 
com <- c('SAR116_Cluster_1300', 'SAR116_Cluster_2430')
counts <- apply(final.blc[(final.blc$Organism == 'SAR11_'),30:59], 2, sum)
suc <- 100*(final.blc[com[1], 30:59]/counts)
icl <- 100*(final.blc[com[2], 30:59]/counts)
plot(as.numeric(suc), type='n', axes=F, main='SAR11 SucA_ICL', cex.main=1.6, ylim=c(0, max(as.numeric(c(suc, icl)), na.rm=TRUE)))
rect(xleft = c(0, 9, 18.5, 27), xright=c(4.5, 12.5, 22.5, 31), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(1, 1, 1, 1, 1, 1), col='grey90', border=NA)
lines(as.numeric(suc), col="dodgerblue2", lwd=2)
lines(as.numeric(icl), col="red", lwd=2)
axis(2, las=2, font=2, cex.axis=1.2)
axis(1, font=2, cex.axis=1.2)
box()


########################################### for BLC
icl <- c("SAR11_Cluster_1331", "SAR116_Cluster_2430", "SAR86_Cluster_809", "SAR406_Cluster_162")
suc <- c('SAR11_Cluster_871', 'SAR116_Cluster_1300', 'SAR86_Cluster_375', 'SAR406_Cluster_757')
 
list <- c("SAR11_", "SAR116_", "SAR86", "SAR406")
suc.result <- c(); icl.result <- c()
for(i in 1:length(suc)) {
	counts <- apply(final.blc[(final.blc$Organism == list[i]),30:59], 2, sum)
	spec.suc <- 100*(final.blc[suc[i], 30:59]/counts)
	spec.icl <- 100*(final.blc[icl[i], 30:59]/counts)
	suc.result <- data.frame(rbind(suc.result, spec.suc))
	icl.result <- data.frame(rbind(icl.result, spec.icl))
}

par(mfrow=c(4,2), mar=c(2, 3, 2, 1))
for(j in 1:length(suc)) {
maxn <- max(as.numeric(c(suc.result[suc[j],], icl.result[icl[j],])), na.rm=TRUE)
plot(as.numeric(suc.result[suc[j],]), type='n', axes=F, main=suc[j], cex.main=1.6, ylim=c(0, maxn))
rect(xleft = c(0, 9, 18.5, 27), xright=c(4.5, 12.5, 22.5, 31), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(1, 1, 1, 1, 1, 1), col='grey90', border=NA)
lines(as.numeric(suc.result[suc[j],]), col="dodgerblue2", lwd=2)
lines(as.numeric(icl.result[icl[j],]), col="red", lwd=2)
axis(2, las=2, font=2, cex.axis=1.2)
axis(1, font=2, cex.axis=1.2)
box()
}

######################## For CANON

icl <- c("SAR11_Cluster_1331", "SAR116_Cluster_2430", "SAR86_Cluster_809", "SAR406_Cluster_162")
suc <- c('SAR11_Cluster_871', 'SAR116_Cluster_1300', 'SAR86_Cluster_375', 'SAR406_Cluster_757')
 
list <- c("SAR11_", "SAR116_", "SAR86", "SAR406")
suc.result <- c(); icl.result <- c()
for(i in 1:length(suc)) {
	counts <- apply(final.canon[(final.canon$Organism == list[i]),37:71], 2, sum)
	spec.suc <- 100*(final.canon[suc[i], 37:71]/counts)
	spec.icl <- 100*(final.canon[icl[i], 37:71]/counts)
	suc.result <- data.frame(rbind(suc.result, spec.suc))
	icl.result <- data.frame(rbind(icl.result, spec.icl))
}

par(mfrow=c(4,2), mar=c(2, 3, 2, 1))
for(j in 1:length(suc)) {
maxn <- max(as.numeric(c(suc.result[suc[j],], icl.result[icl[j],])), na.rm=TRUE)
plot(as.numeric(suc.result[suc[j],]), type='n', axes=F, main=suc[j], cex.main=1.6, ylim=c(0, maxn))
rect(xleft = c(0, 4.5, 10.5, 15.5, 22.5, 30.5), xright=c(1.5, 7.5, 12.5, 19.5, 26.5, 33.5), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(100, 100, 100, 100, 100, 100), col='grey90', border=NA)
lines(as.numeric(suc.result[suc[j],]), col="dodgerblue2", lwd=2)
lines(as.numeric(icl.result[icl[j],]), col="red", lwd=2)
axis(2, las=2, font=2, cex.axis=1.2)
axis(1, font=2, cex.axis=1.2)
box()
}

