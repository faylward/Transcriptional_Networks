
#######################################################################################
################################ Calculate Eigengenes #################################
#######################################################################################

##################################################################################################################################################
############################## Load Necessary Packages
library(WGCNA)
library(igraph)

############################ Enable multithreading for WGCNA
enableWGCNAThreads()


########## Load big summary file
# CANON

jpeg(file="abundance.vs.centrality.jpg", quality=100, res=600, height=17, width=7, units="in")
#par(mfrow=c(2,1), mar=c(0.2,2.2,0.2,0), oma=c(0.2,0.2,0.2,0.2))
par(mfrow=c(2,1))
canon <- read.table(file="R_network_pipeline/CANON/canon.final.combined.list3", sep="\t", header=TRUE, row.names=1, quote="")
sum <- apply(canon[,31:65], 1, sum)
conn <- as.numeric(canon$Scaled.Connectivity.TC)
#conn <- as.numeric(canon$Scaled.Connectivity.DESeq)
reg <- lm(conn ~ log(sum))
plot(log(sum), conn, col='dodgerblue4', ylab='Degree centrality (connectedness)', xlab='log transcript abundance (natural log)')
title("California Coast", line = 0.2)
abline(reg, lwd=2, col='darkorange')
mtext("R^2 = 0.005", at=9.7, line=-2, cex=1.3)
#dev.off()

blc <- read.table(file="R_network_pipeline/BioLincs/blc.final.combined.list3", sep="\t", header=TRUE, row.names=1, quote="")
sum <- apply(blc[,28:57], 1, sum)
conn <- as.numeric(blc$Scaled.Connectivity.TC)
#conn <- as.numeric(blc$Scaled.Connectivity.DESeq)
reg <- lm(conn ~ log(sum))
plot(log(sum), conn, col='dodgerblue4', ylab='Degree centrality (connectedness)', xlab='Transcript abundance (natural log)')
title("North Pacific Subtropical Gyre", line = 0.2)
abline(reg, lwd=2, col='darkorange')
rsq <- summary(reg)$r.square
mtext("R^2 = 0.372", at=11, line=-2, cex=1.3)
dev.off()

blc <- read.table(file="R_network_pipeline/BioLincs/blc.final.combined.list3", sep="\t", header=TRUE, row.names=1, quote="")
non_pro <- row.names(blc)[which(blc$Organism != 'Prochlorococcus')]
non_pro_counts <- apply(blc[non_pro, 28:57], 1, sum)
non_pro_conn <- as.numeric(blc[non_pro, 3])
reg <- lm(non_pro_conn ~ log(non_pro_counts))
plot(log(non_pro_counts), non_pro_conn)
abline(reg)

sum <- apply(non_pro_counts, 1, sum)
conn <- as.numeric(blc$Scaled.Connectivity.TC[non_pro,])
reg <- lm(conn ~ log(sum))

sum <- apply(blc[,28:57], 1, sum)
conn <- as.numeric(blc$Scaled.Connectivity.TC)
#conn <- as.numeric(blc$Scaled.Connectivity.DESeq)
reg <- lm(conn ~ log(sum))
plot(log(sum), conn, col='dodgerblue4', ylab='Degree centrality (connectedness)', xlab='Transcript abundance (natural log)')
title("North Pacific Subtropical Gyre", line = 0.2)
abline(reg, lwd=2, col='darkorange')
rsq <- summary(reg)$r.square
mtext("R^2 = 0.372", at=11, line=-2, cex=1.3)



################################# Find common translation transcripts to map against each other


sar11t <- row.names(blc)[which(blc$Organism == 'SAR11_' & blc$Final_Annote == 'Translation')]
sar11c <- row.names(canon)[which(canon$Organism == 'SAR11_' & canon$Final_Annote == 'Translation')]
sar11 <- unique(c(sar11t, sar11c))

roseot <- row.names(blc)[which(blc$Organism == 'Roseo' & blc$Final_Annote == 'Translation')]
roseoc <- row.names(canon)[which(canon$Organism == 'Roseo' & canon$Final_Annote == 'Translation')]
roseo <- unique(c(roseot, roseoc))

sar116t <- row.names(blc)[which(blc$Organism == 'SAR116_' & blc$Final_Annote == 'Translation')]
sar116c <- row.names(canon)[which(canon$Organism == 'SAR116_' & canon$Final_Annote == 'Translation')]
sar116 <- unique(c(sar116t, sar116c))

sar406t <- row.names(blc)[which(blc$Organism == 'SAR406' & blc$Final_Annote == 'Translation')]
sar406c <- row.names(canon)[which(canon$Organism == 'SAR406' & canon$Final_Annote == 'Translation')]
sar406 <- unique(c(sar406t, sar406c))

sar86t <- row.names(blc)[which(blc$Organism == 'SAR86' & blc$Final_Annote == 'Translation')]
sar86c <- row.names(canon)[which(canon$Organism == 'SAR86' & canon$Final_Annote == 'Translation')]
sar86 <- unique(c(sar86t, sar86c))

prot <- row.names(blc)[which(blc$Organism == 'Prochlorococcus' & blc$Final_Annote == 'Translation')]
proc <- row.names(canon)[which(canon$Organism == 'Prochlorococcus' & canon$Final_Annote == 'Translation')]
pro <- unique(c(prot, proc))

names <- c(sar11, roseo, sar116, sar406, sar86, pro)

write.table(names, file='translation.orthologs.list', sep='\t', quote=F, row.names=F, col.names=F)










