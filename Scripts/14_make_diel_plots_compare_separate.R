
#######################################################################################
################################ Calculate Eigengenes #################################
#######################################################################################

##################################################################################################################################################
############################## Load Necessary Packages
library(ggplot2)

########################################################################
############## Load big summary files CANON ############################
########################################################################

diel <- read.table(file="R_network_pipeline/CANON/canon.all.hra.out", sep="\t", header=TRUE, row.names=1, quote="")
diel.transcript <- (diel$Perm.FDR <= 0.1 & diel$Regression.FDR <= 0.1)
diel.2 <- cbind(diel, diel.transcript)

canon <- read.table(file="R_network_pipeline/CANON/canon.final.combined.list3", sep="\t", header=TRUE, row.names=1, quote="")
final.canon <- merge(diel.2, canon, by="row.names", all=T)
row.names(final.canon) <- final.canon$Row.names

########################################################################
########## Load big summary files BioLincs #############################
########################################################################

diel <- read.table(file="R_network_pipeline/BioLincs/all.blc.hra_11_2014", sep="\t", header=TRUE, row.names=1)
diel.transcript <- (diel$Perm.FDR <= 0.1 & diel$Regression.FDR <= 0.1)
diel.2 <- cbind(diel, diel.transcript)

blc <- read.table(file="R_network_pipeline/BioLincs/blc.final.combined.list2", sep="\t", header=TRUE, row.names=1, quote="")
final.blc <- merge(diel.2, blc, by="row.names", all=T)
row.names(final.blc) <- final.blc$Row.names
##########################################################################


canon <- read.table(file="R_network_pipeline/CANON/canon.final.combined.list3", sep="\t", header=TRUE, row.names=1, quote="")
diel.transcript <- (canon$Perm.FDR <= 0.1 & canon$Regression.FDR <= 0.1)
final.canon <- cbind(canon, diel.transcript)

blc <- read.table(file="R_network_pipeline/BioLincs/blc.final.combined.list3", sep="\t", header=TRUE, row.names=1, quote="")
diel.transcript <- (blc$Perm.FDR <= 0.1 & blc$Regression.FDR <= 0.1)
final.blc <- cbind(blc, diel.transcript)

##########################################################################
############# Specific translation plot for CANON ########################
##########################################################################

functions <- c("Translation")
taxa <- c('SAR11_', 'Ostreococcus', 'SAR116_', 'Roseobacter', 'SAR86', 'SAR406', 'ARCTIC', 'SAR92', 'Flavobacteria', 'Euryarchaota')
# Ostreococcus
taxon <- taxa[2]
org <- row.names(final.canon)[which(final.canon$Organism == taxon & final.canon$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.canon[org,]$Peak.Time, as.character(final.canon[org,]$Final_Annote), as.character(final.canon[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '5', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
Ostreo <- Org

# SAR11
taxon <- taxa[1]
org <- row.names(final.canon)[which(final.canon$Organism == taxon & final.canon$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.canon[org,]$Peak.Time, as.character(final.canon[org,]$Final_Annote), as.character(final.canon[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '3', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
SAR11 <- Org

# SAR116
taxon <- taxa[3]
org <- row.names(final.canon)[which(final.canon$Organism == taxon & final.canon$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.canon[org,]$Peak.Time, as.character(final.canon[org,]$Final_Annote), as.character(final.canon[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '3.5', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
SAR116 <- Org

# Roseo
taxon <- taxa[4]
org <- row.names(final.canon)[which(final.canon$Organism == taxon & final.canon$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.canon[org,]$Peak.Time, as.character(final.canon[org,]$Final_Annote), as.character(final.canon[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '4.5', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
roseo <- Org

# SAR86
taxon <- taxa[5]
org <- row.names(final.canon)[which(final.canon$Organism == taxon & final.canon$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.canon[org,]$Peak.Time, as.character(final.canon[org,]$Final_Annote), as.character(final.canon[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '2.5', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
SAR86 <- Org

# SAR406
taxon <- taxa[6]
org <- row.names(final.canon)[which(final.canon$Organism == taxon & final.canon$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.canon[org,]$Peak.Time, as.character(final.canon[org,]$Final_Annote), as.character(final.canon[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '2', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
SAR406 <- Org

# arctic
taxon <- taxa[7]
org <- row.names(final.canon)[which(final.canon$Organism == taxon & final.canon$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.canon[org,]$Peak.Time, as.character(final.canon[org,]$Final_Annote), as.character(final.canon[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '4', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
arctic <- Org

# SAR92
taxon <- taxa[8]
org <- row.names(final.canon)[which(final.canon$Organism == taxon & final.canon$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.canon[org,]$Peak.Time, as.character(final.canon[org,]$Final_Annote), as.character(final.canon[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '1', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
SAR92 <- Org

# Flavo
taxon <- taxa[9]
org <- row.names(final.canon)[which(final.canon$Organism == taxon & final.canon$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.canon[org,]$Peak.Time, as.character(final.canon[org,]$Final_Annote), as.character(final.canon[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '1.5', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
flavo <- Org

# Euryarchaeota
taxon <- taxa[10]
org <- row.names(final.canon)[which(final.canon$Organism == taxon & final.canon$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.canon[org,]$Peak.Time, as.character(final.canon[org,]$Final_Annote), as.character(final.canon[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '0.5', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
eury <- Org
##################

Ostreo.diel  <- rownames(Ostreo)[which(Ostreo$Diel.Transcripts == TRUE)]
Ostreo.ndiel <- rownames(Ostreo)[which(Ostreo$Diel.Transcripts == FALSE)]
SAR11.diel   <- rownames(SAR11)[which(SAR11$Diel.Transcripts == TRUE)]
SAR11.ndiel  <- rownames(SAR11)[which(SAR11$Diel.Transcripts == FALSE)]
SAR116.diel  <- rownames(SAR116)[which(SAR116$Diel.Transcripts == TRUE)]
SAR116.ndiel <- rownames(SAR116)[which(SAR116$Diel.Transcripts == FALSE)]
roseo.diel   <- rownames(roseo)[which(roseo$Diel.Transcripts == TRUE)]
roseo.ndiel  <- rownames(roseo)[which(roseo$Diel.Transcripts == FALSE)]
SAR86.diel  <- rownames(SAR86)[which(SAR86$Diel.Transcripts == TRUE)]
SAR86.ndiel <- rownames(SAR86)[which(SAR86$Diel.Transcripts == FALSE)]
SAR406.ndiel <- rownames(SAR406)[which(SAR406$Diel.Transcripts == FALSE)]

arctic.diel  <- rownames(arctic)[which(arctic$Diel.Transcripts == TRUE)]
arctic.ndiel <- rownames(arctic)[which(arctic$Diel.Transcripts == FALSE)]

eury.diel  <- rownames(eury)[which(eury$Diel.Transcripts == TRUE)]
eury.ndiel <- rownames(eury)[which(eury$Diel.Transcripts == FALSE)]

SAR92.ndiel <- rownames(SAR92)[which(SAR92$Diel.Transcripts == FALSE)]
flavo.ndiel <- rownames(flavo)[which(flavo$Diel.Transcripts == FALSE)]


#combined <- rbind(Ostreo.diel, Ostreo.ndiel, SAR11.diel, SAR11.ndiel, SAR116.diel, SAR116.ndiel)

y <- geom_vline(aes(xintercept=12), linetype='dashed', colour='grey80')
z <- geom_vline(aes(xintercept=18), linetype='dashed', colour='grey80')

a <- geom_rect(data=Ostreo[Ostreo.ndiel,], aes(xmin=0, xmax=6.75, ymin=0, ymax=5.5), fill='grey90', colour=NA)
b <- geom_rect(data=Ostreo[Ostreo.ndiel,], aes(xmin=19.33, xmax=24.5, ymin=0, ymax=5.5), fill='grey90', colour=NA)
c <- geom_point(data=Ostreo[Ostreo.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")
d <- geom_point(data=Ostreo[Ostreo.diel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.65, colour="green4")
e <- geom_point(data=  SAR11[SAR11.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")
f <- geom_point(data=SAR11[SAR11.diel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.65, colour="red")
g <- geom_point(data=SAR116[SAR116.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")
h <- geom_point(data=SAR116[SAR116.diel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.65, colour="dodgerblue2")
i <- geom_point(data=roseo[roseo.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")
j <- geom_point(data=roseo[roseo.diel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.65, colour="brown")
k <- geom_point(data=SAR86[SAR86.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")
l <- geom_point(data=SAR86[SAR86.diel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.65, colour="gold")
m <- geom_point(data=SAR406[SAR406.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")

q <- geom_point(data=eury[eury.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")
r <- geom_point(data=eury[eury.diel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.75, colour="purple")

s <- geom_point(data=arctic[arctic.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")
t <- geom_point(data=arctic[arctic.diel,],  position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.75, colour="darkorange")

u <- geom_point(data=flavo[flavo.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")
v <- geom_point(data=SAR92[SAR92.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")

n <- scale_x_continuous(breaks=c(0, 6, 12, 18, 24))
o <- scale_y_continuous(breaks=c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), labels=c('GII Euryarchaeota', 'SAR92', 'Flavobacteria', 'SAR406', 'SAR86', 'SAR11', 'SAR116', 'ARCTIC96-BD19', 'Roseobacter', 'Ostreococcus'))
p <- theme(axis.text.y = element_text(colour = 'grey30', size=9, face='bold'), axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(colour='grey60', size=12), panel.background = element_rect(fill='white', colour='white'), legend.position="none")
jpeg("translation.canon.jpg", height=4, width=10, qual=100, res=400, units='in')
ggplot() + a + b + y + z + c + d + e + f + g + h + i + j + k + l + m + n + o + p + q + r + s + t + u + v
dev.off() 

##############################################################################################################
########################################### End ##############################################################
##############################################################################################################

##########################################################################
########### Specific translation plot for BioLincs #######################
##########################################################################

functions <- c("Translation")
taxa <- c('SAR11_', 'Prochlorococcus', 'SAR116_', 'Roseo', 'SAR86', 'SAR406', 'SAR324')
# Ostreococcus
taxon <- taxa[2]
org <- row.names(final.blc)[which(final.blc$Organism == taxon & final.blc$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.blc[org,]$Peak.Time, as.character(final.blc[org,]$Final_Annote), as.character(final.blc[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '3.5', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
Ostreo <- Org

# SAR11
taxon <- taxa[1]
org <- row.names(final.blc)[which(final.blc$Organism == taxon & final.blc$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.blc[org,]$Peak.Time, as.character(final.blc[org,]$Final_Annote), as.character(final.blc[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '2', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
SAR11 <- Org

# SAR116
taxon <- taxa[3]
org <- row.names(final.blc)[which(final.blc$Organism == taxon & final.blc$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.blc[org,]$Peak.Time, as.character(final.blc[org,]$Final_Annote), as.character(final.blc[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '2.5', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
SAR116 <- Org

# Roseo
taxon <- taxa[4]
org <- row.names(final.blc)[which(final.blc$Organism == taxon & final.blc$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.blc[org,]$Peak.Time, as.character(final.blc[org,]$Final_Annote), as.character(final.blc[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '3', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
roseo <- Org

# SAR86
taxon <- taxa[5]
org <- row.names(final.blc)[which(final.blc$Organism == taxon & final.blc$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.blc[org,]$Peak.Time, as.character(final.blc[org,]$Final_Annote), as.character(final.blc[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '1.5', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
SAR86 <- Org

# SAR406
taxon <- taxa[6]
org <- row.names(final.blc)[which(final.blc$Organism == taxon & final.blc$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.blc[org,]$Peak.Time, as.character(final.blc[org,]$Final_Annote), as.character(final.blc[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '1', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
SAR406 <- Org

# SAR324
taxon <- taxa[7]
org <- row.names(final.blc)[which(final.blc$Organism == taxon & final.blc$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.blc[org,]$Peak.Time, as.character(final.blc[org,]$Final_Annote), as.character(final.blc[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '0.5', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
SAR324 <- Org

##################

Ostreo.diel  <- rownames(Ostreo)[which(Ostreo$Diel.Transcripts == TRUE)]
Ostreo.ndiel <- rownames(Ostreo)[which(Ostreo$Diel.Transcripts == FALSE)]
SAR11.diel   <- rownames(SAR11)[which(SAR11$Diel.Transcripts == TRUE)]
SAR11.ndiel  <- rownames(SAR11)[which(SAR11$Diel.Transcripts == FALSE)]
SAR116.diel  <- rownames(SAR116)[which(SAR116$Diel.Transcripts == TRUE)]
SAR116.ndiel <- rownames(SAR116)[which(SAR116$Diel.Transcripts == FALSE)]
roseo.diel   <- rownames(roseo)[which(roseo$Diel.Transcripts == TRUE)]
roseo.ndiel  <- rownames(roseo)[which(roseo$Diel.Transcripts == FALSE)]
SAR86.diel  <- rownames(SAR86)[which(SAR86$Diel.Transcripts == TRUE)]
SAR86.ndiel <- rownames(SAR86)[which(SAR86$Diel.Transcripts == FALSE)]
SAR406.ndiel <- rownames(SAR406)[which(SAR406$Diel.Transcripts == FALSE)]
SAR406.diel <- rownames(SAR406)[which(SAR406$Diel.Transcripts == TRUE)]

SAR324.ndiel <- rownames(SAR324)[which(SAR324$Diel.Transcripts == FALSE)]

#combined <- rbind(Ostreo.diel, Ostreo.ndiel, SAR11.diel, SAR11.ndiel, SAR116.diel, SAR116.ndiel)

#### Change sunset time for BioLincs
y <- geom_vline(aes(xintercept=12), linetype='dashed', colour='grey80')
z <- geom_vline(aes(xintercept=18), linetype='dashed', colour='grey80')

a <- geom_rect(data=Ostreo[Ostreo.ndiel,1:2], aes(xmin=0, xmax=6.33, ymin=0, ymax=4), fill='grey90', colour=NA)
b <- geom_rect(data=Ostreo[Ostreo.ndiel,], aes(xmin=18.5, xmax=24.5, ymin=0, ymax=4), fill='grey90', colour=NA)
c <- geom_point(data=Ostreo[Ostreo.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")
d <- geom_point(data=Ostreo[Ostreo.diel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.65, colour="green4")
e <- geom_point(data=  SAR11[SAR11.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")
f <- geom_point(data=SAR11[SAR11.diel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.65, colour="red")
g <- geom_point(data=SAR116[SAR116.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")
h <- geom_point(data=SAR116[SAR116.diel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.65, colour="dodgerblue2")
i <- geom_point(data=roseo[roseo.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")
j <- geom_point(data=roseo[roseo.diel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.65, colour="brown")
k <- geom_point(data=SAR86[SAR86.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")
l <- geom_point(data=SAR86[SAR86.diel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.65, colour="gold")
m <- geom_point(data=SAR406[SAR406.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")
n <- geom_point(data=SAR406[SAR406.diel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.65, colour="darkblue")

s <- geom_point(data=SAR324[SAR324.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")

p <- theme(axis.text.y = element_text(colour = 'grey30', size=9, face='bold'), axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(colour='grey60', size=12), panel.background = element_rect(fill='white', colour='white'), legend.position="none")
q <- scale_x_continuous(breaks=c(0, 6, 12, 18, 24))
r <- scale_y_continuous(breaks=c(0.5, 1, 1.5, 2, 2.5, 3, 3.5), labels=c('SAR324', 'SAR406', 'SAR86', 'SAR11', 'SAR116', 'Roseobacter', 'Prochlorococ'))
jpeg("translation.blc.jpg", height=3.3, width=10, qual=100, res=400, units='in')
ggplot() + y + z +  a + b + c + d + e + f + g + h + i + j + k + l + m + n + s + p + q + r
dev.off() 

##############################################################################################################
########################################### End ##############################################################
##############################################################################################################


##########################################################################
############# Specific translation plot for CANON ########################
##########################################################################

functions <- c("Oxidative_phosphorylation")
taxa <- c('SAR11_', 'Ostreococcus', 'SAR116_', 'Roseobacter', 'SAR86', 'SAR406', 'ARCTIC', 'SAR92', 'Flavobacteria', 'Euryarchaota')

# Ostreococcus
taxon <- taxa[2]
org <- row.names(final.canon)[which(final.canon$Organism == taxon & final.canon$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.canon[org,]$Peak.Time, as.character(final.canon[org,]$Final_Annote), as.character(final.canon[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '5', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
Ostreo <- Org

# SAR11
taxon <- taxa[1]
org <- row.names(final.canon)[which(final.canon$Organism == taxon & final.canon$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.canon[org,]$Peak.Time, as.character(final.canon[org,]$Final_Annote), as.character(final.canon[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '3', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
SAR11 <- Org

# SAR116
taxon <- taxa[3]
org <- row.names(final.canon)[which(final.canon$Organism == taxon & final.canon$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.canon[org,]$Peak.Time, as.character(final.canon[org,]$Final_Annote), as.character(final.canon[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '3.5', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
SAR116 <- Org

# Roseo
taxon <- taxa[4]
org <- row.names(final.canon)[which(final.canon$Organism == taxon & final.canon$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.canon[org,]$Peak.Time, as.character(final.canon[org,]$Final_Annote), as.character(final.canon[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '4.5', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
roseo <- Org

# SAR86
taxon <- taxa[5]
org <- row.names(final.canon)[which(final.canon$Organism == taxon & final.canon$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.canon[org,]$Peak.Time, as.character(final.canon[org,]$Final_Annote), as.character(final.canon[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '2.5', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
SAR86 <- Org

# SAR406
taxon <- taxa[6]
org <- row.names(final.canon)[which(final.canon$Organism == taxon & final.canon$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.canon[org,]$Peak.Time, as.character(final.canon[org,]$Final_Annote), as.character(final.canon[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '2', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
SAR406 <- Org

# arctic
taxon <- taxa[7]
org <- row.names(final.canon)[which(final.canon$Organism == taxon & final.canon$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.canon[org,]$Peak.Time, as.character(final.canon[org,]$Final_Annote), as.character(final.canon[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '4', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
arctic <- Org

# SAR92
taxon <- taxa[8]
org <- row.names(final.canon)[which(final.canon$Organism == taxon & final.canon$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.canon[org,]$Peak.Time, as.character(final.canon[org,]$Final_Annote), as.character(final.canon[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '1', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
SAR92 <- Org

# Flavo
taxon <- taxa[9]
org <- row.names(final.canon)[which(final.canon$Organism == taxon & final.canon$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.canon[org,]$Peak.Time, as.character(final.canon[org,]$Final_Annote), as.character(final.canon[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '1.5', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
flavo <- Org

# Euryarchaeota
taxon <- taxa[10]
org <- row.names(final.canon)[which(final.canon$Organism == taxon & final.canon$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.canon[org,]$Peak.Time, as.character(final.canon[org,]$Final_Annote), as.character(final.canon[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '0.5', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
eury <- Org
##################

Ostreo.diel  <- rownames(Ostreo)[which(Ostreo$Diel.Transcripts == TRUE)]
Ostreo.ndiel <- rownames(Ostreo)[which(Ostreo$Diel.Transcripts == FALSE)]
SAR11.diel   <- rownames(SAR11)[which(SAR11$Diel.Transcripts == TRUE)]
SAR11.ndiel  <- rownames(SAR11)[which(SAR11$Diel.Transcripts == FALSE)]
SAR116.diel  <- rownames(SAR116)[which(SAR116$Diel.Transcripts == TRUE)]
SAR116.ndiel <- rownames(SAR116)[which(SAR116$Diel.Transcripts == FALSE)]
roseo.diel   <- rownames(roseo)[which(roseo$Diel.Transcripts == TRUE)]
roseo.ndiel  <- rownames(roseo)[which(roseo$Diel.Transcripts == FALSE)]
SAR86.diel  <- rownames(SAR86)[which(SAR86$Diel.Transcripts == TRUE)]
SAR86.ndiel <- rownames(SAR86)[which(SAR86$Diel.Transcripts == FALSE)]
SAR406.ndiel <- rownames(SAR406)[which(SAR406$Diel.Transcripts == FALSE)]

arctic.diel  <- rownames(arctic)[which(arctic$Diel.Transcripts == TRUE)]
arctic.ndiel <- rownames(arctic)[which(arctic$Diel.Transcripts == FALSE)]

eury.ndiel <- rownames(eury)[which(eury$Diel.Transcripts == FALSE)]

SAR92.ndiel <- rownames(SAR92)[which(SAR92$Diel.Transcripts == FALSE)]
flavo.ndiel <- rownames(flavo)[which(flavo$Diel.Transcripts == FALSE)]



#combined <- rbind(Ostreo.diel, Ostreo.ndiel, SAR11.diel, SAR11.ndiel, SAR116.diel, SAR116.ndiel)
y <- geom_vline(aes(xintercept=12), linetype='dashed', colour='grey80')
z <- geom_vline(aes(xintercept=18), linetype='dashed', colour='grey80')

a <- geom_rect(data=Ostreo[Ostreo.ndiel,1:2], aes(xmin=0, xmax=6.75, ymin=0, ymax=5.5), fill='grey90', colour=NA)
b <- geom_rect(data=Ostreo[Ostreo.ndiel,], aes(xmin=19.33, xmax=24.5, ymin=0, ymax=5.5), fill='grey90', colour=NA)
c <- geom_point(data=Ostreo[Ostreo.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")
d <- geom_point(data=Ostreo[Ostreo.diel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.65, colour="green4")
e <- geom_point(data=  SAR11[SAR11.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")
f <- geom_point(data=SAR11[SAR11.diel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.65, colour="red")
g <- geom_point(data=SAR116[SAR116.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")
h <- geom_point(data=SAR116[SAR116.diel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.65, colour="dodgerblue2")
i <- geom_point(data=roseo[roseo.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")
j <- geom_point(data=roseo[roseo.diel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.65, colour="brown")
k <- geom_point(data=SAR86[SAR86.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")
l <- geom_point(data=SAR86[SAR86.diel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.65, colour="gold")
m <- geom_point(data=SAR406[SAR406.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")

q <- geom_point(data=eury[eury.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")

s <- geom_point(data=arctic[arctic.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")
t <- geom_point(data=arctic[arctic.diel,],  position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.75, colour="darkorange")

u <- geom_point(data=flavo[flavo.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")
v <- geom_point(data=SAR92[SAR92.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")

n <- scale_x_continuous(breaks=c(0, 6, 12, 18, 24))
o <- scale_y_continuous(breaks=c(0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5), labels=c('GII Euryarchaeota', 'SAR92', 'Flavobacteria', 'SAR406', 'SAR86', 'SAR11', 'SAR116', 'ARCTIC96-BD19', 'Roseobacter', 'Ostreococcus'))
p <- theme(axis.text.y = element_text(colour = 'grey30', size=9, face='bold'), axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(colour='grey60', size=12), panel.background = element_rect(fill='white', colour='white'), legend.position="none")
jpeg("OxPhos.CANON.jpg", height=4, width=10, qual=100, res=400, units='in')
ggplot() + y + z + a + b + c + d + e + f + g + h + i + j + k + l + m + q + s + t + u + v + n + o + p
dev.off() 

##############################################################################################################
########################################### End ##############################################################
##############################################################################################################

##########################################################################
########### Specific translation plot for BioLincs #######################
##########################################################################

functions <- c("Oxidative_phosphorylation")
taxa <- c('SAR11_', 'Prochlorococcus', 'SAR116_', 'Roseo', 'SAR86', 'SAR406', 'SAR324')
# Ostreococcus
taxon <- taxa[2]
org <- row.names(final.blc)[which(final.blc$Organism == taxon & final.blc$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.blc[org,]$Peak.Time, as.character(final.blc[org,]$Final_Annote), as.character(final.blc[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '3.5', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
Ostreo <- Org

# SAR11
taxon <- taxa[1]
org <- row.names(final.blc)[which(final.blc$Organism == taxon & final.blc$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.blc[org,]$Peak.Time, as.character(final.blc[org,]$Final_Annote), as.character(final.blc[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '2', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
SAR11 <- Org

# SAR116
taxon <- taxa[3]
org <- row.names(final.blc)[which(final.blc$Organism == taxon & final.blc$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.blc[org,]$Peak.Time, as.character(final.blc[org,]$Final_Annote), as.character(final.blc[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '2.5', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
SAR116 <- Org

# Roseo
taxon <- taxa[4]
org <- row.names(final.blc)[which(final.blc$Organism == taxon & final.blc$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.blc[org,]$Peak.Time, as.character(final.blc[org,]$Final_Annote), as.character(final.blc[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '3', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
roseo <- Org

# SAR86
taxon <- taxa[5]
org <- row.names(final.blc)[which(final.blc$Organism == taxon & final.blc$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.blc[org,]$Peak.Time, as.character(final.blc[org,]$Final_Annote), as.character(final.blc[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '1.5', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
SAR86 <- Org

# SAR406
taxon <- taxa[6]
org <- row.names(final.blc)[which(final.blc$Organism == taxon & final.blc$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.blc[org,]$Peak.Time, as.character(final.blc[org,]$Final_Annote), as.character(final.blc[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '1', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
SAR406 <- Org

# SAR324
taxon <- taxa[7]
org <- row.names(final.blc)[which(final.blc$Organism == taxon & final.blc$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.blc[org,]$Peak.Time, as.character(final.blc[org,]$Final_Annote), as.character(final.blc[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '0.5', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
SAR324 <- Org

##################

Ostreo.diel  <- rownames(Ostreo)[which(Ostreo$Diel.Transcripts == TRUE)]
Ostreo.ndiel <- rownames(Ostreo)[which(Ostreo$Diel.Transcripts == FALSE)]
SAR11.diel   <- rownames(SAR11)[which(SAR11$Diel.Transcripts == TRUE)]
SAR11.ndiel  <- rownames(SAR11)[which(SAR11$Diel.Transcripts == FALSE)]
SAR116.diel  <- rownames(SAR116)[which(SAR116$Diel.Transcripts == TRUE)]
SAR116.ndiel <- rownames(SAR116)[which(SAR116$Diel.Transcripts == FALSE)]
roseo.diel   <- rownames(roseo)[which(roseo$Diel.Transcripts == TRUE)]
roseo.ndiel  <- rownames(roseo)[which(roseo$Diel.Transcripts == FALSE)]
SAR86.diel  <- rownames(SAR86)[which(SAR86$Diel.Transcripts == TRUE)]
SAR86.ndiel <- rownames(SAR86)[which(SAR86$Diel.Transcripts == FALSE)]
SAR406.ndiel <- rownames(SAR406)[which(SAR406$Diel.Transcripts == FALSE)]
SAR406.diel <- rownames(SAR406)[which(SAR406$Diel.Transcripts == TRUE)]

SAR324.ndiel <- rownames(SAR324)[which(SAR324$Diel.Transcripts == FALSE)]

#combined <- rbind(Ostreo.diel, Ostreo.ndiel, SAR11.diel, SAR11.ndiel, SAR116.diel, SAR116.ndiel)

#### Change sunset time for BioLincs
y <- geom_vline(aes(xintercept=12), linetype='dashed', colour='grey80')
z <- geom_vline(aes(xintercept=18), linetype='dashed', colour='grey80')

a <- geom_rect(data=Ostreo[Ostreo.ndiel,1:2], aes(xmin=0, xmax=6.33, ymin=0, ymax=4), fill='grey90', colour=NA)
b <- geom_rect(data=Ostreo[Ostreo.ndiel,], aes(xmin=18.5, xmax=24.5, ymin=0, ymax=4), fill='grey90', colour=NA)
c <- geom_point(data=Ostreo[Ostreo.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")
d <- geom_point(data=Ostreo[Ostreo.diel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.65, colour="green4")
e <- geom_point(data=  SAR11[SAR11.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")
f <- geom_point(data=SAR11[SAR11.diel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.65, colour="red")
g <- geom_point(data=SAR116[SAR116.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")
h <- geom_point(data=SAR116[SAR116.diel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.65, colour="dodgerblue2")
i <- geom_point(data=roseo[roseo.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")
j <- geom_point(data=roseo[roseo.diel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.65, colour="brown")
k <- geom_point(data=SAR86[SAR86.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")
l <- geom_point(data=SAR86[SAR86.diel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.65, colour="gold")
m <- geom_point(data=SAR406[SAR406.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")
n <- geom_point(data=SAR406[SAR406.diel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.65, colour="darkblue")

s <- geom_point(data=SAR324[SAR324.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")

p <- theme(axis.text.y = element_text(colour = 'grey30', size=9, face='bold'), axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(colour='grey60', size=12), panel.background = element_rect(fill='white', colour='white'), legend.position="none")
q <- scale_x_continuous(breaks=c(0, 6, 12, 18, 24))
r <- scale_y_continuous(breaks=c(0.5, 1, 1.5, 2, 2.5, 3, 3.5), labels=c('SAR324', 'SAR406', 'SAR86', 'SAR11', 'SAR116', 'Roseobacter', 'Prochlorococ'))
jpeg("OxPhos.blc.jpg", height=3.3, width=10, qual=100, res=400, units='in')
ggplot() + y + z +  a + b + c + d + e + f + g + h + i + j + k + l + m + n + s + p + q + r
dev.off() 


##############################################################################################################
########################################### End ##############################################################
##############################################################################################################


##########################################################################
############# Specific translation plot for CANON ########################
##########################################################################

functions <- c("Photosynthesis")
taxa <- c('SAR11_', 'Ostreococcus', 'SAR116_', 'Roseobacter', 'SAR86', 'SAR406')
# Ostreococcus
taxon <- taxa[2]
org <- row.names(final.canon)[which(final.canon$Organism == taxon & final.canon$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.canon[org,]$Peak.Time, as.character(final.canon[org,]$Final_Annote), as.character(final.canon[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '1', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
Ostreo <- Org

Ostreo.diel  <- rownames(Ostreo)[which(Ostreo$Diel.Transcripts == TRUE)]
Ostreo.ndiel <- rownames(Ostreo)[which(Ostreo$Diel.Transcripts == FALSE)]

#combined <- rbind(Ostreo.diel, Ostreo.ndiel, SAR11.diel, SAR11.ndiel, SAR116.diel, SAR116.ndiel)
y <- geom_vline(aes(xintercept=12), linetype='dashed', colour='grey80')
z <- geom_vline(aes(xintercept=18), linetype='dashed', colour='grey80')

a <- geom_rect(data=Ostreo[Ostreo.ndiel,], aes(xmin=0, xmax=6.75, ymin=0, ymax=2), fill='grey90', colour=NA)
b <- geom_rect(data=Ostreo[Ostreo.ndiel,], aes(xmin=19.33, xmax=24.5, ymin=0, ymax=2), fill='grey90', colour=NA)
c <- geom_point(data=Ostreo[Ostreo.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")
d <- geom_point(data=Ostreo[Ostreo.diel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.65, colour="green4")
e <- scale_x_continuous(breaks=c(0, 6, 12, 18, 24))
f <- scale_y_continuous(breaks=c(1), labels=c('Ostreococcus'))
g <- theme(axis.text.y = element_text(colour = 'grey30', size=9, face='bold'), axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(colour='grey60', size=12), panel.background = element_rect(fill='white', colour='white'), legend.position="none")
jpeg("Photosynthesis.CANON.jpg", height=1, width=10, qual=100, res=400, units='in')
ggplot() + y + z + a + b + c + d + e + f + g
dev.off() 

##############################################################################################################
########################################### End ##############################################################
##############################################################################################################

##########################################################################
########### Specific translation plot for BioLincs #######################
##########################################################################

functions <- c("Photosynthesis")
taxa <- c('SAR11_', 'Prochlorococcus', 'SAR116_', 'Roseo', 'SAR86', 'SAR406')
# Ostreococcus
taxon <- taxa[2]
org <- row.names(final.blc)[which(final.blc$Organism == taxon & final.blc$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.blc[org,]$Peak.Time, as.character(final.blc[org,]$Final_Annote), as.character(final.blc[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '1', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
Ostreo <- Org

Ostreo.diel  <- rownames(Ostreo)[which(Ostreo$Diel.Transcripts == TRUE)]
Ostreo.ndiel <- rownames(Ostreo)[which(Ostreo$Diel.Transcripts == FALSE)]

#combined <- rbind(Ostreo.diel, Ostreo.ndiel, SAR11.diel, SAR11.ndiel, SAR116.diel, SAR116.ndiel)

#### Change sunset time for BioLincs
y <- geom_vline(aes(xintercept=12), linetype='dashed', colour='grey80')
z <- geom_vline(aes(xintercept=18), linetype='dashed', colour='grey80')

a <- geom_rect(data=Ostreo[Ostreo.ndiel,1:2], aes(xmin=0, xmax=6.33, ymin=0, ymax=2), fill='grey90', colour=NA)
b <- geom_rect(data=Ostreo[Ostreo.ndiel,], aes(xmin=18.5, xmax=24.5, ymin=0, ymax=2), fill='grey90', colour=NA)
c <- geom_point(data=Ostreo[Ostreo.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")
d <- geom_point(data=Ostreo[Ostreo.diel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.65, colour="green4")
e <- theme(axis.text.y = element_text(colour = 'grey30', size=9, face='bold'), axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(colour='grey60', size=12), panel.background = element_rect(fill='white', colour='white'), legend.position="none")
f <- scale_x_continuous(breaks=c(0, 6, 12, 18, 24))
g <- scale_y_continuous(breaks=c(1), labels=c('Prochlorococ'))
jpeg("Photosynthesis.blc.jpg", height=1, width=10, qual=100, res=400, units='in')
ggplot() + y + z + a + b + c + d + e + f + g
dev.off() 

##############################################################################################################
########################################### End ##############################################################
##############################################################################################################

##########################################################################
############# Specific translation plot for CANON ########################
##########################################################################

functions <- c("Carbon_fixation")
taxa <- c('SAR11_', 'Ostreococcus', 'SAR116_', 'Roseobacter', 'SAR86', 'SAR406')
# Ostreococcus
taxon <- taxa[2]
org <- row.names(final.canon)[which(final.canon$Organism == taxon & final.canon$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.canon[org,]$Peak.Time, as.character(final.canon[org,]$Final_Annote), as.character(final.canon[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '1', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
Ostreo <- Org

Ostreo.diel  <- rownames(Ostreo)[which(Ostreo$Diel.Transcripts == TRUE)]
Ostreo.ndiel <- rownames(Ostreo)[which(Ostreo$Diel.Transcripts == FALSE)]

#combined <- rbind(Ostreo.diel, Ostreo.ndiel, SAR11.diel, SAR11.ndiel, SAR116.diel, SAR116.ndiel)
y <- geom_vline(aes(xintercept=12), linetype='dashed', colour='grey80')
z <- geom_vline(aes(xintercept=18), linetype='dashed', colour='grey80')

a <- geom_rect(data=Ostreo[Ostreo.ndiel,], aes(xmin=0, xmax=6.75, ymin=0, ymax=2), fill='grey90', colour=NA)
b <- geom_rect(data=Ostreo[Ostreo.ndiel,], aes(xmin=19.33, xmax=24.5, ymin=0, ymax=2), fill='grey90', colour=NA)
c <- geom_point(data=Ostreo[Ostreo.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")
d <- geom_point(data=Ostreo[Ostreo.diel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.65, colour="green4")
e <- scale_x_continuous(breaks=c(0, 6, 12, 18, 24))
f <- scale_y_continuous(breaks=c(1), labels=c('Ostreococcus'))
g <- theme(axis.text.y = element_text(colour = 'grey30', size=9, face='bold'), axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(colour='grey60', size=12), panel.background = element_rect(fill='white', colour='white'), legend.position="none")
jpeg("Carbon_fixation.CANON.jpg", height=1, width=10, qual=100, res=400, units='in')
ggplot() +  y + z + a + b + c + d + e + f + g
dev.off() 

##############################################################################################################
########################################### End ##############################################################
##############################################################################################################

##########################################################################
########### Specific translation plot for BioLincs #######################
##########################################################################

functions <- c("Carbon_fixation")
taxa <- c('SAR11_', 'Prochlorococcus', 'SAR116_', 'Roseo', 'SAR86', 'SAR406')
# Ostreococcus
taxon <- taxa[2]
org <- row.names(final.blc)[which(final.blc$Organism == taxon & final.blc$Final_Annote %in% functions)]
Org <- data.frame(cbind(final.blc[org,]$Peak.Time, as.character(final.blc[org,]$Final_Annote), as.character(final.blc[org,]$diel.transcript)))
row.names(Org) <- org; colnames(Org) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Org$Annotation <- as.numeric(gsub(functions, '1', Org$Annotation))
Org$Peak.Time <- as.numeric(as.character(Org$Peak.Time))
Ostreo <- Org

Ostreo.diel  <- rownames(Ostreo)[which(Ostreo$Diel.Transcripts == TRUE)]
Ostreo.ndiel <- rownames(Ostreo)[which(Ostreo$Diel.Transcripts == FALSE)]

#combined <- rbind(Ostreo.diel, Ostreo.ndiel, SAR11.diel, SAR11.ndiel, SAR116.diel, SAR116.ndiel)

#### Change sunset time for BioLincs
y <- geom_vline(aes(xintercept=12), linetype='dashed', colour='grey80')
z <- geom_vline(aes(xintercept=18), linetype='dashed', colour='grey80')

a <- geom_rect(data=Ostreo[Ostreo.ndiel,1:2], aes(xmin=0, xmax=6.33, ymin=0, ymax=2), fill='grey90', colour=NA)
b <- geom_rect(data=Ostreo[Ostreo.ndiel,], aes(xmin=18.5, xmax=24.5, ymin=0, ymax=2), fill='grey90', colour=NA)
c <- geom_point(data=Ostreo[Ostreo.ndiel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.45, colour="grey60")
d <- geom_point(data=Ostreo[Ostreo.diel,], position = position_jitter(w = 0, h = 0.1), aes(Peak.Time, Annotation), shape=16, size=3, alpha=0.65, colour="green4")
e <- theme(axis.text.y = element_text(colour = 'grey30', size=9, face='bold'), axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(colour='grey60', size=12), panel.background = element_rect(fill='white', colour='white'), legend.position="none")
f <- scale_x_continuous(breaks=c(0, 6, 12, 18, 24))
g <- scale_y_continuous(breaks=c(1), labels=c('Prochlorococ'))
jpeg("Carbon_fixation.blc.jpg", height=1, width=10, qual=100, res=400, units='in')
ggplot() + y + z + a + b + c + d + e + f + g
dev.off() 

##############################################################################################################
########################################### End ##############################################################
##############################################################################################################



