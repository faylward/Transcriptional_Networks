
#######################################################################################
################################ Calculate Eigengenes #################################
#######################################################################################

##################################################################################################################################################
############################## Load Necessary Packages
library(ggplot2)

########################################################################
############## Load big summary files CANON ############################
########################################################################

diel <- read.table(file="R_network_pipeline/CANON/All.canon.hra.out", sep="\t", header=TRUE, row.names=1)
diel.transcript <- (diel$Perm.FDR <= 0.1 & diel$Regression.FDR <= 0.1)
diel.2 <- cbind(diel, diel.transcript)

canon <- read.table(file="R_network_pipeline/CANON/canon.combined.summary2", sep="\t", header=TRUE, row.names=1)
final.canon <- merge(diel.2, canon, by="row.names", all=T)
row.names(final.canon) <- final.canon$Row.names

########################################################################
########## Load big summary files BioLincs #############################
########################################################################

diel <- read.table(file="R_network_pipeline/BioLincs/All.blc.hra", sep="\t", header=TRUE, row.names=1)
diel.transcript <- (diel$Perm.FDR <= 0.1 & diel$Regression.FDR <= 0.1)
diel.2 <- cbind(diel, diel.transcript)

blc <- read.table(file="R_network_pipeline/BioLincs/blc.combined.summary", sep="\t", header=TRUE, row.names=1, quote="")
final.blc <- merge(diel.2, blc, by="row.names", all=T)
row.names(final.blc) <- final.blc$Row.names

#################################################################################
################### Ostreococcus/Prochlorococcus ################################
#################################################################################
functions <- c("Amino_Acid_metabolism", "Carbohydrate_Metabolism", "Carbon_Fixation", "Cell_Cycle", "Cofactor_vitamin_metabolism", "DNA/Chromosome_Processing", "Fatty_acid_metabolism", "Nucleotide_metabolism", "Oxidative_phosphorylation", "Photosynthesis", "Protein_RNA_Processing", "Signal_Transduction", "Transcription", "Translation", "Transporters")

# Ostreococcus
diel.canon <- row.names(final.canon)[which(final.canon$Organism == "Ostreococcus" & final.canon$Final_Annote %in% functions)]
Ostreo <- data.frame(cbind(final.canon[diel.canon,]$Peak.Time, as.character(final.canon[diel.canon,]$Final_Annote), as.character(final.canon[diel.canon,]$diel.transcript)))
row.names(Ostreo) <- diel.canon; colnames(Ostreo) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Ostreo$Peak.Time <- as.numeric(as.character(Ostreo$Peak.Time))

Ostreo$Annotation <- gsub("Protein_RNA_Processing", "29", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Fatty_acid_metabolism", "27", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Translation", "25", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Oxidative_phosphorylation", "23", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Carbon_Fixation", "21", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Cofactor_vitamin_metabolism", "19", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Amino_Acid_metabolism", "17", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Photosynthesis", "15", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Nucleotide_metabolism", "13", Ostreo$Annotation)
Ostreo$Annotation <- gsub("DNA/Chromosome_Processing", "11", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Cell_Cycle", "09", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Carbohydrate_Metabolism", "07", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Transcription", "05", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Signal_Transduction", "03", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Transporters", "01", Ostreo$Annotation)

Ostreo$Annotation <- as.numeric(Ostreo$Annotation)

# Prochlorococcus
diel.blc <- row.names(final.blc)[which(final.blc$Organism == "Pro" & final.blc$Final_Annote %in% functions)]
pro <- data.frame(cbind(final.blc[diel.blc,]$Peak.Time, as.character(final.blc[diel.blc,]$Final_Annote), as.character(final.blc[diel.blc,]$diel.transcript)))
row.names(pro) <- diel.blc; colnames(pro) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
pro$Peak.Time <- as.numeric(as.character(pro$Peak.Time))
pro$Annotation <- as.character(pro$Annotation)

pro$Annotation <- gsub("Protein_RNA_Processing", "30", pro$Annotation)
pro$Annotation <- gsub("Fatty_acid_metabolism", "28", pro$Annotation)
pro$Annotation <- gsub("Translation", "26", pro$Annotation)
pro$Annotation <- gsub("Oxidative_phosphorylation", "24", pro$Annotation)
pro$Annotation <- gsub("Carbon_Fixation", "22", pro$Annotation)
pro$Annotation <- gsub("Cofactor_vitamin_metabolism", "20", pro$Annotation)
pro$Annotation <- gsub("Amino_Acid_metabolism", "18", pro$Annotation)
pro$Annotation <- gsub("Photosynthesis", "16", pro$Annotation)
pro$Annotation <- gsub("Nucleotide_metabolism", "14", pro$Annotation)
pro$Annotation <- gsub("DNA/Chromosome_Processing", "12", pro$Annotation)
pro$Annotation <- gsub("Cell_Cycle", "10", pro$Annotation)
pro$Annotation <- gsub("Carbohydrate_Metabolism", "08", pro$Annotation)
pro$Annotation <- gsub("Transcription", "06", pro$Annotation)
pro$Annotation <- gsub("Signal_Transduction", "04", pro$Annotation)
pro$Annotation <- gsub("Transporters", "02", pro$Annotation)

pro$Annotation <- as.numeric(pro$Annotation)

##################
diel.b <- row.names(pro)[which(pro$Diel.Transcripts == TRUE)]
not.diel.b <- row.names(pro)[which(pro$Diel.Transcripts == FALSE)]

diel.c <- row.names(Ostreo)[which(Ostreo$Diel.Transcripts == TRUE)]
not.diel.c <- row.names(Ostreo)[which(Ostreo$Diel.Transcripts == FALSE)]

######## Labels
functions <- c("Protein/RNA Processing", "Fatty Acid Metabolism", "Translation", "Oxidative Phosphorylation", "Carbon Fixation", "Cofactor Vitamin Metabolism", "Amino Acid Metabolism", "Photosynthesis", "Nucleotide Metabolism", "DNA/Chromosome Processing", "Cell Cycle", "Carbohydrate Metabolism", "Transcription", "Signal Transduction", "Transporters")
placement=c(2.5, 4.5, 6.5, 8.5, 10.5, 12.5, 14.5, 16.5, 18.5, 20.5, 22.5, 24.5, 26.5, 28.5, 30.5)
labels <- placement - 1
#Full figure for SI
#postscript("Ostreo.dotplot2.eps", height=6.5, width=13, horizontal = FALSE, onefile = FALSE)

a <- geom_rect(data=Ostreo[not.diel.c,1:2], aes(xmin=-0.5, xmax=6.5, ymin=0, ymax=30.5), fill='grey87', colour=NA)
b <- geom_rect(data=Ostreo[not.diel.c,1:2], aes(xmin=19.5, xmax=24.5, ymin=0, ymax=30.5), fill='grey87', colour=NA)
c <- geom_rect(data=Ostreo[not.diel.c,1:2], aes(xmin=18.5, xmax=19.5, ymin=0, ymax=30.5), fill='grey92', colour=NA)
d <- geom_point(data=Ostreo[not.diel.c,1:2], position = position_jitter(w = 0, h = 0.2), aes(Peak.Time, Annotation), alpha=0.45, shape=16, size=4, colour="grey60")
e <- geom_point(data=pro[not.diel.b,1:2], position = position_jitter(w = 0, h = 0.2), aes(Peak.Time, Annotation), alpha=0.45, shape=16, size=4, colour="grey60")
f <- geom_point(data=Ostreo[diel.c,1:2], position = position_jitter(w = 0, h = 0.2), aes(Peak.Time, Annotation), alpha=0.45, shape=16, size=4, colour="dodgerblue2")
g <- geom_point(data=pro[diel.b,1:2], position = position_jitter(w = 0, h = 0.2), aes(Peak.Time, Annotation), alpha=0.45, shape=16, size=4, colour="darkorange")
h <- theme(axis.text.y = element_text(colour = 'grey30', size=10, face='bold'), axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(colour='grey60', size=12), panel.background = element_rect(fill='white', colour='white'), legend.position="none")
i <- scale_x_continuous(breaks=c(0, 6, 12, 18, 24))
j <- scale_y_continuous(breaks=labels, labels=rev(functions))
k <- geom_vline(aes(xintercept=12), linetype='dashed', colour='grey80')
l <- geom_vline(aes(xintercept=18), linetype='dashed', colour='grey80')
m <- geom_hline(aes(yintercept=placement), linetype='dashed', colour='grey80')

jpeg("Ostreo.dotplot2.jpg", height=8.5, width=10, qual=100, res=400, units='in')
ggplot() + a + b + c + d + e + f + g + h + i + j + k + l + m
dev.off()

#############################################################################################
######################################### End ###############################################
#############################################################################################

##########################################

#SAR116
#CANON
functions <- c("Amino_Acid_metabolism", "Carbohydrate_Metabolism", "Citrate_Cycle", "Cofactor_vitamin_metabolism", "DNA/Chromosome_Processing", "Nucleotide_metabolism", "Oxidative_phosphorylation", "Protein_RNA_Processing", "Signal_Transduction", "Transcription", "Translation", "Transporters", "Chemotaxis_and_Flagellar_Assembly")

diel.t <- row.names(final.canon)[which(final.canon$Organism == "SAR116_" & final.canon$Final_Annote %in% functions)]
Ostreo <- data.frame(cbind(final.canon[diel.t,]$Peak.Time, as.character(final.canon[diel.t,]$Final_Annote), as.character(final.canon[diel.t,]$diel.transcript)))
row.names(Ostreo) <- diel.t; colnames(Ostreo) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Ostreo$Peak.Time <- as.numeric(as.character(Ostreo$Peak.Time))
Ostreo$Annotation <- as.character(Ostreo$Annotation)
Ostreo$Annotation <- gsub("Transporters", "25", Ostreo$Annotation)
Ostreo$Annotation <- gsub("DNA/Chromosome_Processing", "23", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Signal_Transduction", "21", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Chemotaxis_and_Flagellar_Assembly", "19", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Amino_Acid_metabolism", "17", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Translation", "15", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Transcription", "13", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Oxidative_phosphorylation", "11", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Protein_RNA_Processing", "09", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Citrate_Cycle", "07", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Cofactor_vitamin_metabolism", "05", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Nucleotide_metabolism", "03", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Carbohydrate_Metabolism", "01", Ostreo$Annotation)

Ostreo$Annotation <- as.numeric(Ostreo$Annotation)

# BioLincs
diel.blc <- row.names(final.blc)[which(final.blc$Organism == "SAR116_" & final.blc$Final_Annote %in% functions)]
pro <- data.frame(cbind(final.blc[diel.blc,]$Peak.Time, as.character(final.blc[diel.blc,]$Final_Annote), as.character(final.blc[diel.blc,]$diel.transcript)))
row.names(pro) <- diel.blc; colnames(pro) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
pro$Peak.Time <- as.numeric(as.character(pro$Peak.Time))
pro$Annotation <- as.character(pro$Annotation)
pro$Annotation <- gsub("Transporters", "26", pro$Annotation)
pro$Annotation <- gsub("DNA/Chromosome_Processing", "24", pro$Annotation)
pro$Annotation <- gsub("Signal_Transduction", "22", pro$Annotation)
pro$Annotation <- gsub("Chemotaxis_and_Flagellar_Assembly", "20", pro$Annotation)
pro$Annotation <- gsub("Amino_Acid_metabolism", "18", pro$Annotation)
pro$Annotation <- gsub("Translation", "16", pro$Annotation)
pro$Annotation <- gsub("Transcription", "14", pro$Annotation)
pro$Annotation <- gsub("Oxidative_phosphorylation", "12", pro$Annotation)
pro$Annotation <- gsub("Protein_RNA_Processing", "10", pro$Annotation)
pro$Annotation <- gsub("Citrate_Cycle", "08", pro$Annotation)
pro$Annotation <- gsub("Cofactor_vitamin_metabolism", "06", pro$Annotation)
pro$Annotation <- gsub("Nucleotide_metabolism", "04", pro$Annotation)
pro$Annotation <- gsub("Carbohydrate_Metabolism", "02", pro$Annotation)

pro$Annotation <- as.numeric(pro$Annotation)

diel.c <- row.names(Ostreo)[which(Ostreo$Diel.Transcripts == TRUE)]
not.diel.c <- row.names(Ostreo)[which(Ostreo$Diel.Transcripts == FALSE)]
diel.b <- row.names(pro)[which(pro$Diel.Transcripts == TRUE)]
not.diel.b <- row.names(pro)[which(pro$Diel.Transcripts == FALSE)]

functions <- c( "Carbohydrate Metabolism",  "Nucleotide Metabolism", "Cofactor Vitamin Metabolism", "Citrate Cycle",  "Protein/RNA Processing", "Oxidative Phosphorylation", "Transcription", "Translation",  "Amino Acid Metabolism", "Chemotaxis and Flagellar Assembly", "Signal Transduction", "DNA/Chromosome Processing",  "Transporters")

placement=c(2.5, 4.5, 6.5, 8.5, 10.5, 12.5, 14.5, 16.5, 18.5, 20.5, 22.5, 24.5, 26.5)
labels <- placement - 1
#Full figure for SI
#postscript("Ostreo.dotplot2.eps", height=6.5, width=13, horizontal = FALSE, onefile = FALSE)

a <- geom_rect(data=Ostreo[not.diel.c,1:2], aes(xmin=-0.5, xmax=6.5, ymin=0, ymax=26.5), fill='grey87', colour=NA)
b <- geom_rect(data=Ostreo[not.diel.c,1:2], aes(xmin=19.5, xmax=24.5, ymin=0, ymax=26.5), fill='grey87', colour=NA)
c <- geom_rect(data=Ostreo[not.diel.c,1:2], aes(xmin=18.5, xmax=19.5, ymin=0, ymax=26.5), fill='grey92', colour=NA)
d <- geom_point(data=Ostreo[not.diel.c,1:2], position = position_jitter(w = 0, h = 0.2), aes(Peak.Time, Annotation), alpha=0.45, shape=16, size=4, colour="grey60")
e <- geom_point(data=pro[not.diel.b,1:2], position = position_jitter(w = 0, h = 0.2), aes(Peak.Time, Annotation), alpha=0.45, shape=16, size=4, colour="grey60")
f <- geom_point(data=Ostreo[diel.c,1:2], position = position_jitter(w = 0, h = 0.2), aes(Peak.Time, Annotation), alpha=0.45, shape=16, size=4, colour="dodgerblue2")
g <- geom_point(data=pro[diel.b,1:2], position = position_jitter(w = 0, h = 0.2), aes(Peak.Time, Annotation), alpha=0.45, shape=16, size=4, colour="darkorange")
h <- theme(axis.text.y = element_text(colour = 'grey30', size=10, face='bold'), axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(colour='grey60', size=12), panel.background = element_rect(fill='white', colour='white'), legend.position="none")
i <- scale_x_continuous(breaks=c(0, 6, 12, 18, 24))
j <- scale_y_continuous(breaks=labels, labels=functions)
k <- geom_vline(aes(xintercept=12), linetype='dashed', colour='grey80')
l <- geom_vline(aes(xintercept=18), linetype='dashed', colour='grey80')
m <- geom_hline(aes(yintercept=placement), linetype='dashed', colour='grey80')

jpeg("SAR116.dotplot2.jpg", height=8.5, width=10, qual=100, res=400, units='in')
ggplot() + a + b + c + d + e + f + g + h + i + j + k + l + m
dev.off()


#############################################################################################
######################################### End ###############################################
#############################################################################################

#SAR11
functions <- c("Amino_Acid_metabolism", "Carbohydrate_Metabolism", "Citrate_Cycle", "Cofactor_vitamin_metabolism", "DNA/Chromosome_Processing", "Fatty_acid_metabolism", "Nucleotide_metabolism", "Oxidative_phosphorylation", "Protein_RNA_Processing", "Signal_Transduction", "Transcription", "Translation", "Transporters")
diel.t <- row.names(final.canon)[which(final.canon$Organism == "SAR11_" & final.canon$Final_Annote %in% functions)]
Ostreo <- data.frame(cbind(final.canon[diel.t,]$Peak.Time, as.character(final.canon[diel.t,]$Final_Annote), as.character(final.canon[diel.t,]$diel.transcript)))
row.names(Ostreo) <- diel.t; colnames(Ostreo) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Ostreo$Peak.Time <- as.numeric(as.character(Ostreo$Peak.Time))
Ostreo$Annotation <- gsub("Fatty_acid_metabolism", "25", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Transporters", "23", Ostreo$Annotation)
Ostreo$Annotation <- gsub("DNA/Chromosome_Processing", "21", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Signal_Transduction", "19", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Amino_Acid_metabolism", "17", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Translation", "15", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Transcription", "13", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Oxidative_phosphorylation", "11", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Protein_RNA_Processing", "09", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Citrate_Cycle", "07", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Cofactor_vitamin_metabolism", "05", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Nucleotide_metabolism", "03", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Carbohydrate_Metabolism", "01", Ostreo$Annotation)

Ostreo$Annotation <- as.numeric(Ostreo$Annotation)
# BioLincs
diel.blc <- row.names(final.blc)[which(final.blc$Organism == "SAR11_" & final.blc$Final_Annote %in% functions)]
pro <- data.frame(cbind(final.blc[diel.blc,]$Peak.Time, as.character(final.blc[diel.blc,]$Final_Annote), as.character(final.blc[diel.blc,]$diel.transcript)))
row.names(pro) <- diel.blc; colnames(pro) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
pro$Peak.Time <- as.numeric(as.character(pro$Peak.Time))
pro$Annotation <- gsub("Fatty_acid_metabolism", "26", pro$Annotation)
pro$Annotation <- gsub("Transporters", "24", pro$Annotation)
pro$Annotation <- gsub("DNA/Chromosome_Processing", "22", pro$Annotation)
pro$Annotation <- gsub("Signal_Transduction", "20", pro$Annotation)
pro$Annotation <- gsub("Amino_Acid_metabolism", "18", pro$Annotation)
pro$Annotation <- gsub("Translation", "16", pro$Annotation)
pro$Annotation <- gsub("Transcription", "14", pro$Annotation)
pro$Annotation <- gsub("Oxidative_phosphorylation", "12", pro$Annotation)
pro$Annotation <- gsub("Protein_RNA_Processing", "10", pro$Annotation)
pro$Annotation <- gsub("Citrate_Cycle", "08", pro$Annotation)
pro$Annotation <- gsub("Cofactor_vitamin_metabolism", "06", pro$Annotation)
pro$Annotation <- gsub("Nucleotide_metabolism", "04", pro$Annotation)
pro$Annotation <- gsub("Carbohydrate_Metabolism", "02", pro$Annotation)

pro$Annotation <- as.numeric(pro$Annotation)

diel.c <- row.names(Ostreo)[which(Ostreo$Diel.Transcripts == TRUE)]
not.diel.c <- row.names(Ostreo)[which(Ostreo$Diel.Transcripts == FALSE)]

diel.b <- row.names(pro)[which(pro$Diel.Transcripts == TRUE)]
not.diel.b <- row.names(pro)[which(pro$Diel.Transcripts == FALSE)]
######################################

functions <- c( "Carbohydrate Metabolism",  "Nucleotide Metabolism", "Cofactor Vitamin Metabolism", "Citrate_Cycle",  "Protein/RNA Processing", "Oxidative Phosphorylation", "Transcription", "Translation",  "Amino Acid Metabolism", "Signal Transduction", "DNA/Chromosome Processing",  "Transporters", "Fatty_acid_metabolism")

placement=c(2.5, 4.5, 6.5, 8.5, 10.5, 12.5, 14.5, 16.5, 18.5, 20.5, 22.5, 24.5, 26.5)
labels <- placement - 1
#Full figure for SI
#postscript("Ostreo.dotplot2.eps", height=6.5, width=13, horizontal = FALSE, onefile = FALSE)

a <- geom_rect(data=Ostreo[not.diel.c,1:2], aes(xmin=-0.5, xmax=6.5, ymin=0, ymax=26.5), fill='grey87', colour=NA)
b <- geom_rect(data=Ostreo[not.diel.c,1:2], aes(xmin=19.5, xmax=24.5, ymin=0, ymax=26.5), fill='grey87', colour=NA)
c <- geom_rect(data=Ostreo[not.diel.c,1:2], aes(xmin=18.5, xmax=19.5, ymin=0, ymax=26.5), fill='grey92', colour=NA)
d <- geom_point(data=Ostreo[not.diel.c,1:2], position = position_jitter(w = 0, h = 0.2), aes(Peak.Time, Annotation), alpha=0.45, shape=16, size=4, colour="grey60")
e <- geom_point(data=pro[not.diel.b,1:2], position = position_jitter(w = 0, h = 0.2), aes(Peak.Time, Annotation), alpha=0.45, shape=16, size=4, colour="grey60")
f <- geom_point(data=Ostreo[diel.c,1:2], position = position_jitter(w = 0, h = 0.2), aes(Peak.Time, Annotation), alpha=0.45, shape=16, size=4, colour="dodgerblue2")
g <- geom_point(data=pro[diel.b,1:2], position = position_jitter(w = 0, h = 0.2), aes(Peak.Time, Annotation), alpha=0.45, shape=16, size=4, colour="darkorange")
h <- theme(axis.text.y = element_text(colour = 'grey30', size=10, face='bold'), axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(colour='grey60', size=12), panel.background = element_rect(fill='white', colour='white'), legend.position="none")
i <- scale_x_continuous(breaks=c(0, 6, 12, 18, 24))
j <- scale_y_continuous(breaks=labels, labels=functions)
k <- geom_vline(aes(xintercept=12), linetype='dashed', colour='grey80')
l <- geom_vline(aes(xintercept=18), linetype='dashed', colour='grey80')
m <- geom_hline(aes(yintercept=placement), linetype='dashed', colour='grey80')

jpeg("SAR11.dotplot2.jpg", height=8.5, width=10, qual=100, res=400, units='in')
ggplot() + a + b + c + d + e + f + g + h + i + j + k + l + m
dev.off()

######################################
#############################################################################################
######################################### End ###############################################
#############################################################################################

#Roseo
functions <- c("Amino_Acid_metabolism", "Carbohydrate_Metabolism", "Chemotaxis_and_Flagellar_Assembly", "Citrate_Cycle", "Cofactor_vitamin_metabolism", "DNA/Chromosome_Processing", "Nucleotide_metabolism", "Oxidative_phosphorylation", "Protein_RNA_Processing", "Signal_Transduction", "Sulfur_metabolism", "Transcription", "Translation", "Transporters")
diel.t <- row.names(final.canon)[which(final.canon$Organism == "Roseobacter" & final.canon$Final_Annote %in% functions)]
Ostreo <- data.frame(cbind(final.canon[diel.t,]$Peak.Time, as.character(final.canon[diel.t,]$Final_Annote), as.character(final.canon[diel.t,]$diel.transcript)))
row.names(Ostreo) <- diel.t; colnames(Ostreo) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Ostreo$Peak.Time <- as.numeric(as.character(Ostreo$Peak.Time))
Ostreo$Annotation <- gsub("Transporters", "25", Ostreo$Annotation)
Ostreo$Annotation <- gsub("DNA/Chromosome_Processing", "23", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Signal_Transduction", "21", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Chemotaxis_and_Flagellar_Assembly", "19", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Amino_Acid_metabolism", "17", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Translation", "15", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Transcription", "13", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Oxidative_phosphorylation", "11", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Protein_RNA_Processing", "09", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Citrate_Cycle", "07", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Cofactor_vitamin_metabolism", "05", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Nucleotide_metabolism", "03", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Carbohydrate_Metabolism", "01", Ostreo$Annotation)

Ostreo$Annotation <- as.numeric(Ostreo$Annotation)

diel.blc <- row.names(final.blc)[which(final.blc$Organism == "Roseo" & final.blc$Final_Annote %in% functions)]
pro <- data.frame(cbind(final.blc[diel.blc,]$Peak.Time, as.character(final.blc[diel.blc,]$Final_Annote), as.character(final.blc[diel.blc,]$diel.transcript)))
row.names(pro) <- diel.blc; colnames(pro) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
pro$Peak.Time <- as.numeric(as.character(pro$Peak.Time))
pro$Annotation <- gsub("Transporters", "24", pro$Annotation)
pro$Annotation <- gsub("DNA/Chromosome_Processing", "24", pro$Annotation)
pro$Annotation <- gsub("Signal_Transduction", "22", pro$Annotation)
pro$Annotation <- gsub("Chemotaxis_and_Flagellar_Assembly", "20", pro$Annotation)
pro$Annotation <- gsub("Amino_Acid_metabolism", "18", pro$Annotation)
pro$Annotation <- gsub("Translation", "16", pro$Annotation)
pro$Annotation <- gsub("Transcription", "14", pro$Annotation)
pro$Annotation <- gsub("Oxidative_phosphorylation", "12", pro$Annotation)
pro$Annotation <- gsub("Protein_RNA_Processing", "10", pro$Annotation)
pro$Annotation <- gsub("Citrate_Cycle", "08", pro$Annotation)
pro$Annotation <- gsub("Cofactor_vitamin_metabolism", "06", pro$Annotation)
pro$Annotation <- gsub("Nucleotide_metabolism", "04", pro$Annotation)
pro$Annotation <- gsub("Carbohydrate_Metabolism", "02", pro$Annotation)

pro$Annotation <- as.numeric(pro$Annotation)

diel.c <- row.names(Ostreo)[which(Ostreo$Diel.Transcripts == TRUE)]
not.diel.c <- row.names(Ostreo)[which(Ostreo$Diel.Transcripts == FALSE)]
diel.b <- row.names(pro)[which(pro$Diel.Transcripts == TRUE)]
not.diel.b <- row.names(pro)[which(pro$Diel.Transcripts == FALSE)]

functions <- c( "Carbohydrate Metabolism",  "Nucleotide Metabolism", "Cofactor Vitamin Metabolism", "Citrate Cycle",  "Protein/RNA Processing", "Oxidative Phosphorylation", "Transcription", "Translation",  "Amino Acid Metabolism", "Chemotaxis and Flagellar Assembly", "Signal Transduction", "DNA/Chromosome Processing",  "Transporters")

placement=c(2.5, 4.5, 6.5, 8.5, 10.5, 12.5, 14.5, 16.5, 18.5, 20.5, 22.5, 24.5, 26.5)
labels <- placement - 1
#Full figure for SI
#postscript("Ostreo.dotplot2.eps", height=6.5, width=13, horizontal = FALSE, onefile = FALSE)

a <- geom_rect(data=Ostreo[not.diel.c,1:2], aes(xmin=-0.5, xmax=6.5, ymin=0, ymax=26.5), fill='grey87', colour=NA)
b <- geom_rect(data=Ostreo[not.diel.c,1:2], aes(xmin=19.5, xmax=24.5, ymin=0, ymax=26.5), fill='grey87', colour=NA)
c <- geom_rect(data=Ostreo[not.diel.c,1:2], aes(xmin=18.5, xmax=19.5, ymin=0, ymax=26.5), fill='grey92', colour=NA)
d <- geom_point(data=Ostreo[not.diel.c,1:2], position = position_jitter(w = 0, h = 0.2), aes(Peak.Time, Annotation), alpha=0.45, shape=16, size=4, colour="grey60")
e <- geom_point(data=pro[not.diel.b,1:2], position = position_jitter(w = 0, h = 0.2), aes(Peak.Time, Annotation), alpha=0.45, shape=16, size=4, colour="grey60")
f <- geom_point(data=Ostreo[diel.c,1:2], position = position_jitter(w = 0, h = 0.2), aes(Peak.Time, Annotation), alpha=0.45, shape=16, size=4, colour="dodgerblue2")
g <- geom_point(data=pro[diel.b,1:2], position = position_jitter(w = 0, h = 0.2), aes(Peak.Time, Annotation), alpha=0.45, shape=16, size=4, colour="darkorange")
h <- theme(axis.text.y = element_text(colour = 'grey30', size=10, face='bold'), axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(colour='grey60', size=12), panel.background = element_rect(fill='white', colour='white'), legend.position="none")
i <- scale_x_continuous(breaks=c(0, 6, 12, 18, 24))
j <- scale_y_continuous(breaks=labels, labels=functions)
k <- geom_vline(aes(xintercept=12), linetype='dashed', colour='grey80')
l <- geom_vline(aes(xintercept=18), linetype='dashed', colour='grey80')
m <- geom_hline(aes(yintercept=placement), linetype='dashed', colour='grey80')

jpeg("Roseobacter.dotplot2.jpg", height=8.5, width=10, qual=100, res=400, units='in')
ggplot() + a + b + c + d + e + f + g + h + i + j + k + l + m
dev.off()

#############################################################################################
######################################### End ###############################################
#############################################################################################


#SAR86
functions <- c("Amino_Acid_metabolism", "Carbohydrate_Metabolism", "Citrate_Cycle", "Cofactor_vitamin_metabolism", "DNA/Chromosome_Processing", "Nucleotide_metabolism", "Oxidative_phosphorylation", "Protein_RNA_Processing", "Transcription", "Translation", "Transporters", "Signal_Transduction")
diel.t <- row.names(final.canon)[which(final.canon$Organism == "SAR86" & final.canon$Final_Annote %in% functions)]
Ostreo <- data.frame(cbind(final.canon[diel.t,]$Peak.Time, as.character(final.canon[diel.t,]$Final_Annote), as.character(final.canon[diel.t,]$diel.transcript)))
row.names(Ostreo) <- diel.t; colnames(Ostreo) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Ostreo$Peak.Time <- as.numeric(as.character(Ostreo$Peak.Time))
Ostreo$Annotation <- gsub("Transporters", "23", Ostreo$Annotation)
Ostreo$Annotation <- gsub("DNA/Chromosome_Processing", "21", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Signal_Transduction", "19", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Amino_Acid_metabolism", "17", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Translation", "15", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Transcription", "13", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Oxidative_phosphorylation", "11", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Protein_RNA_Processing", "09", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Citrate_Cycle", "07", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Cofactor_vitamin_metabolism", "05", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Nucleotide_metabolism", "03", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Carbohydrate_Metabolism", "01", Ostreo$Annotation)

Ostreo$Annotation <- as.numeric(Ostreo$Annotation)

diel.blc <- row.names(final.blc)[which(final.blc$Organism == "SAR86" & final.blc$Final_Annote %in% functions)]
pro <- data.frame(cbind(final.blc[diel.blc,]$Peak.Time, as.character(final.blc[diel.blc,]$Final_Annote), as.character(final.blc[diel.blc,]$diel.transcript)))
row.names(pro) <- diel.blc; colnames(pro) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
pro$Peak.Time <- as.numeric(as.character(pro$Peak.Time))
pro$Annotation <- gsub("Transporters", "24", pro$Annotation)
pro$Annotation <- gsub("DNA/Chromosome_Processing", "22", pro$Annotation)
pro$Annotation <- gsub("Signal_Transduction", "20", pro$Annotation)
pro$Annotation <- gsub("Amino_Acid_metabolism", "18", pro$Annotation)
pro$Annotation <- gsub("Translation", "16", pro$Annotation)
pro$Annotation <- gsub("Transcription", "14", pro$Annotation)
pro$Annotation <- gsub("Oxidative_phosphorylation", "12", pro$Annotation)
pro$Annotation <- gsub("Protein_RNA_Processing", "10", pro$Annotation)
pro$Annotation <- gsub("Citrate_Cycle", "08", pro$Annotation)
pro$Annotation <- gsub("Cofactor_vitamin_metabolism", "06", pro$Annotation)
pro$Annotation <- gsub("Nucleotide_metabolism", "04", pro$Annotation)
pro$Annotation <- gsub("Carbohydrate_Metabolism", "02", pro$Annotation)

pro$Annotation <- as.numeric(pro$Annotation)

diel.c <- row.names(Ostreo)[which(Ostreo$Diel.Transcripts == TRUE)]
not.diel.c <- row.names(Ostreo)[which(Ostreo$Diel.Transcripts == FALSE)]

diel.b <- row.names(pro)[which(pro$Diel.Transcripts == TRUE)]
not.diel.b <- row.names(pro)[which(pro$Diel.Transcripts == FALSE)]

functions <- c( "Carbohydrate Metabolism",  "Nucleotide Metabolism", "Cofactor Vitamin Metabolism", "Citrate Cycle", "Protein/RNA Processing", "Oxidative Phosphorylation", "Transcription", "Translation",  "Amino Acid Metabolism", "Signal Transduction", "DNA/Chromosome Processing",  "Transporters")

placement=c(2.5, 4.5, 6.5, 8.5, 10.5, 12.5, 14.5, 16.5, 18.5, 20.5, 22.5, 24.5)
labels <- placement - 1
#Full figure for SI
#postscript("Ostreo.dotplot2.eps", height=6.5, width=13, horizontal = FALSE, onefile = FALSE)

a <- geom_rect(data=Ostreo[not.diel.c,1:2], aes(xmin=-0.5, xmax=6.5, ymin=0, ymax=24.5), fill='grey87', colour=NA)
b <- geom_rect(data=Ostreo[not.diel.c,1:2], aes(xmin=19.5, xmax=24.5, ymin=0, ymax=24.5), fill='grey87', colour=NA)
c <- geom_rect(data=Ostreo[not.diel.c,1:2], aes(xmin=18.5, xmax=19.5, ymin=0, ymax=24.5), fill='grey92', colour=NA)
d <- geom_point(data=Ostreo[not.diel.c,1:2], position = position_jitter(w = 0, h = 0.2), aes(Peak.Time, Annotation), alpha=0.45, shape=16, size=4, colour="grey60")
e <- geom_point(data=pro[not.diel.b,1:2], position = position_jitter(w = 0, h = 0.2), aes(Peak.Time, Annotation), alpha=0.45, shape=16, size=4, colour="grey60")
f <- geom_point(data=Ostreo[diel.c,1:2], position = position_jitter(w = 0, h = 0.2), aes(Peak.Time, Annotation), alpha=0.45, shape=16, size=4, colour="dodgerblue2")
g <- geom_point(data=pro[diel.b,1:2], position = position_jitter(w = 0, h = 0.2), aes(Peak.Time, Annotation), alpha=0.45, shape=16, size=4, colour="darkorange")
h <- theme(axis.text.y = element_text(colour = 'grey30', size=10, face='bold'), axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(colour='grey60', size=12), panel.background = element_rect(fill='white', colour='white'), legend.position="none")
i <- scale_x_continuous(breaks=c(0, 6, 12, 18, 24))
j <- scale_y_continuous(breaks=labels, labels=functions)
k <- geom_vline(aes(xintercept=12), linetype='dashed', colour='grey80')
l <- geom_vline(aes(xintercept=18), linetype='dashed', colour='grey80')
m <- geom_hline(aes(yintercept=placement), linetype='dashed', colour='grey80')

jpeg("SAR86.dotplot2.jpg", height=8.5, width=10, qual=100, res=400, units='in')
ggplot() + a + b + c + d + e + f + g + h + i + j + k + l + m
dev.off()


########################################################################
########## Load big summary files BioLincs #############################
########################################################################

diel <- read.table(file="R_network_pipeline/BioLincs/All.blc.hra", sep="\t", header=TRUE, row.names=1)
diel.transcript <- (diel$Perm.FDR <= 0.1 & diel$Regression.FDR <= 0.1)
diel.2 <- cbind(diel, diel.transcript)

blc <- read.table(file="R_network_pipeline/BioLincs/blc.combined.summary", sep="\t", header=TRUE, row.names=1, quote="")
final.blc <- merge(diel.2, blc, by="row.names", all=T)
row.names(final.blc) <- final.blc$Row.names

#SAR11

functions <- c("Amino_Acid_metabolism", "Carbohydrate_Metabolism", "Citrate_Cycle", "Cofactor_vitamin_metabolism", "DNA/Chromosome_Processing", "Fatty_acid_metabolism", "Nucleotide_metabolism", "Oxidative_phosphorylation", "Protein_RNA_Processing", "Signal_Transduction", "Transcription", "Translation", "Transporters")

diel.t <- row.names(final.blc)[which(final.blc$Organism == "SAR11_" & final.blc$Final_Annote %in% functions)]
Ostreo <- data.frame(cbind(final.blc[diel.t,]$Peak.Time, as.character(final.blc[diel.t,]$Final_Annote), as.character(final.blc[diel.t,]$diel.transcript)))
row.names(Ostreo) <- diel.t; colnames(Ostreo) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Ostreo$Peak.Time <- as.numeric(as.character(Ostreo$Peak.Time))
Ostreo$Annotation <- gsub("Transporters", "15_Transporters", Ostreo$Annotation)
Ostreo$Annotation <- gsub("DNA/Chromosome_Processing", "14_DNA/Chromosome_Processing", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Signal_Transduction", "13_Signal_Transduction", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Chemotaxis_and_Flagellar_Assembly", "12_Chemotaxis_and_Flagellar_Assembly", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Amino_Acid_metabolism", "11_Amino_Acid_metabolism", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Translation", "10_Translation", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Transcription", "09_Transcription", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Oxidative_phosphorylation", "08_Oxidative_phosphorylation", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Protein_RNA_Processing", "07_Protein_RNA_Processing", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Citrate_Cycle", "06_Citrate_Cycle", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Cofactor_vitamin_metabolism", "05_Cofactor_vitamin_metabolism", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Nucleotide_metabolism", "04_Nucleotide_metabolism", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Carbohydrate_Metabolism", "03_Carbohydrate_Metabolism", Ostreo$Annotation)

diel <- row.names(Ostreo)[which(Ostreo$Diel.Transcripts == TRUE)]
not.diel <- row.names(Ostreo)[which(Ostreo$Diel.Transcripts == FALSE)]

postscript(file="R_network_pipeline/Figures/SAR11.blc.dotplot.eps", height=4.5, width=9, horizontal = FALSE, onefile = FALSE)
ggplot() + geom_point(data=Ostreo[not.diel,], aes(Peak.Time, Annotation), shape=21, size=5, colour="grey60") + geom_rect(data=Ostreo, aes(xmin=-0.5, xmax=6.5, ymin=-Inf, ymax=+Inf), fill='grey90', colour=NA) + geom_rect(data=Ostreo, aes(xmin=18.5, xmax=24.5, ymin=-Inf, ymax=+Inf), fill='grey90', colour=NA) + geom_point(data=Ostreo[not.diel,], aes(Peak.Time, Annotation), shape=21, size=5, colour="grey60") + geom_point(shape=21, size=5.3) + geom_point(shape=21, size=5) + geom_point(data=Ostreo[diel,], aes(Peak.Time, Annotation), shape=21, size=5, colour="blue") + theme(axis.text.y = element_text(colour = 'grey40', size=10), axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(colour='grey60', size=12), panel.background = element_rect(fill='white', colour='white'), legend.position="none") + scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) + geom_vline(aes(xintercept=12), linetype='dashed', colour='grey70') + geom_vline(aes(xintercept=18), linetype='dashed', colour='grey70')
dev.off()

########################SAR116 

functions <- c("Amino_Acid_metabolism", "Carbohydrate_Metabolism", "Chemotaxis_and_Flagellar_Assembly", "Citrate_Cycle", "Cofactor_vitamin_metabolism", "DNA/Chromosome_Processing", "Nucleotide_metabolism", "Oxidative_phosphorylation", "Protein_RNA_Processing", "Signal_Transduction", "Sulfur_metabolism", "Transcription", "Translation", "Transporters")
diel.t <- row.names(final.blc)[which(final.blc$Organism == "SAR116_" & final.blc$Final_Annote %in% functions)]
Ostreo <- data.frame(cbind(final.blc[diel.t,]$Peak.Time, as.character(final.blc[diel.t,]$Final_Annote), as.character(final.blc[diel.t,]$diel.transcript)))
row.names(Ostreo) <- diel.t; colnames(Ostreo) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Ostreo$Peak.Time <- as.numeric(as.character(Ostreo$Peak.Time))
Ostreo$Annotation <- gsub("Transporters", "15_Transporters", Ostreo$Annotation)
Ostreo$Annotation <- gsub("DNA/Chromosome_Processing", "14_DNA/Chromosome_Processing", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Signal_Transduction", "13_Signal_Transduction", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Chemotaxis_and_Flagellar_Assembly", "12_Chemotaxis_and_Flagellar_Assembly", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Amino_Acid_metabolism", "11_Amino_Acid_metabolism", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Translation", "10_Translation", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Transcription", "09_Transcription", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Oxidative_phosphorylation", "08_Oxidative_phosphorylation", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Protein_RNA_Processing", "07_Protein_RNA_Processing", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Citrate_Cycle", "06_Citrate_Cycle", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Cofactor_vitamin_metabolism", "05_Cofactor_vitamin_metabolism", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Nucleotide_metabolism", "04_Nucleotide_metabolism", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Carbohydrate_Metabolism", "03_Carbohydrate_Metabolism", Ostreo$Annotation)

diel <- row.names(Ostreo)[which(Ostreo$Diel.Transcripts == TRUE)]
not.diel <- row.names(Ostreo)[which(Ostreo$Diel.Transcripts == FALSE)]

postscript(file="R_network_pipeline/Figures/SAR116.blc.dotplot.eps", height=4.5, width=9, horizontal = FALSE, onefile = FALSE)
ggplot() + geom_point(data=Ostreo[not.diel,], aes(Peak.Time, Annotation), shape=21, size=5, colour="grey60") + geom_rect(data=Ostreo, aes(xmin=-0.5, xmax=6.5, ymin=-Inf, ymax=+Inf), fill='grey90', colour=NA) + geom_rect(data=Ostreo, aes(xmin=18.5, xmax=24.5, ymin=-Inf, ymax=+Inf), fill='grey90', colour=NA) + geom_point(data=Ostreo[not.diel,], aes(Peak.Time, Annotation), shape=21, size=5, colour="grey60") + geom_point(shape=21, size=5.3) + geom_point(shape=21, size=5) + geom_point(data=Ostreo[diel,], aes(Peak.Time, Annotation), shape=21, size=5, colour="blue") + theme(axis.text.y = element_text(colour = 'grey40', size=10), axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(colour='grey60', size=12), panel.background = element_rect(fill='white', colour='white'), legend.position="none") + scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) + geom_vline(aes(xintercept=12), linetype='dashed', colour='grey70') + geom_vline(aes(xintercept=18), linetype='dashed', colour='grey70')
dev.off()

###################################### Prochlorococcus

functions <- c("Amino_Acid_metabolism", "Carbohydrate_Metabolism", "Carbon_Fixation", "Cell_Cycle", "Cofactor_vitamin_metabolism", "DNA/Chromosome_Processing", "Fatty_acid_metabolism", "Nucleotide_metabolism", "Oxidative_phosphorylation", "Photosynthesis", "Protein_RNA_Processing", "Signal_Transduction", "Transcription", "Translation", "Transporters")

diel.t <- row.names(final.blc)[which(final.blc$Organism == "Pro" & final.blc$Final_Annote %in% functions)]
Ostreo <- data.frame(cbind(final.blc[diel.t,]$Peak.Time, as.character(final.blc[diel.t,]$Final_Annote), as.character(final.blc[diel.t,]$diel.transcript)))
row.names(Ostreo) <- diel.t; colnames(Ostreo) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Ostreo$Peak.Time <- as.numeric(as.character(Ostreo$Peak.Time))
Ostreo$Annotation <- gsub("Protein_RNA_Processing", "15_Protein_RNA_Processing", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Fatty_acid_metabolism", "14_Fatty_acid_metabolism", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Translation", "13_Translation", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Oxidative_phosphorylation", "12_Oxidative_phosphorylation", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Carbon_Fixation", "11_Carbon_Fixation", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Cofactor_vitamin_metabolism", "07_Cofactor_vitamin_metabolism", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Amino_Acid_metabolism", "06_Amino_Acid_metabolism", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Photosynthesis", "10_Photosynthesis", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Nucleotide_metabolism", "08_Nucleotide_metabolism", Ostreo$Annotation)
Ostreo$Annotation <- gsub("DNA/Chromosome_Processing", "09_DNA/Chromosome_Processing", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Cell_Cycle", "05_Cell_Cycle", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Carbohydrate_Metabolism", "04_Carbohydrate_Metabolism", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Transcription", "03_Transcription", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Signal_Transduction", "02_Signal_Transduction", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Transporters", "01_Transporters", Ostreo$Annotation)

diel <- row.names(Ostreo)[which(Ostreo$Diel.Transcripts == TRUE)]
not.diel <- row.names(Ostreo)[which(Ostreo$Diel.Transcripts == FALSE)]

postscript(file="R_network_pipeline/Figures/Pro.dotplot.eps", height=5, width=9, horizontal = FALSE, onefile = FALSE)
ggplot() + geom_point(data=Ostreo[not.diel,], aes(Peak.Time, Annotation), shape=21, size=5, colour="grey60") + geom_rect(data=Ostreo, aes(xmin=-0.5, xmax=6.5, ymin=-Inf, ymax=+Inf), fill='grey90', colour=NA) + geom_rect(data=Ostreo, aes(xmin=18.5, xmax=24.5, ymin=-Inf, ymax=+Inf), fill='grey90', colour=NA) + geom_point(data=Ostreo[not.diel,], aes(Peak.Time, Annotation), shape=21, size=5, colour="grey60") + geom_point(shape=21, size=5.3) + geom_point(shape=21, size=5) + geom_point(data=Ostreo[diel,], aes(Peak.Time, Annotation), shape=21, size=5, colour="blue") + theme(axis.text.y = element_text(colour = 'grey40', size=10), axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(colour='grey60', size=12), panel.background = element_rect(fill='white', colour='white'), legend.position="none") + scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) + geom_vline(aes(xintercept=12), linetype='dashed', colour='grey70') + geom_vline(aes(xintercept=18), linetype='dashed', colour='grey70')
dev.off()

#Roseo
functions <- c("Amino_Acid_metabolism", "Carbohydrate_Metabolism", "Chemotaxis_and_Flagellar_Assembly", "Citrate_Cycle", "Cofactor_vitamin_metabolism", "DNA/Chromosome_Processing", "Nucleotide_metabolism", "Oxidative_phosphorylation", "Protein_RNA_Processing", "Signal_Transduction", "Sulfur_metabolism", "Transcription", "Translation", "Transporters")
diel.t <- row.names(final.blc)[which(final.blc$Organism == "Roseo" & final.blc$Final_Annote %in% functions)]
Ostreo <- data.frame(cbind(final.blc[diel.t,]$Peak.Time, as.character(final.blc[diel.t,]$Final_Annote), as.character(final.blc[diel.t,]$diel.transcript)))
row.names(Ostreo) <- diel.t; colnames(Ostreo) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Ostreo$Peak.Time <- as.numeric(as.character(Ostreo$Peak.Time))
Ostreo$Annotation <- gsub("Transporters", "15_Transporters", Ostreo$Annotation)
Ostreo$Annotation <- gsub("DNA/Chromosome_Processing", "14_DNA/Chromosome_Processing", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Signal_Transduction", "13_Signal_Transduction", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Chemotaxis_and_Flagellar_Assembly", "12_Chemotaxis_and_Flagellar_Assembly", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Amino_Acid_metabolism", "11_Amino_Acid_metabolism", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Translation", "10_Translation", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Transcription", "09_Transcription", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Oxidative_phosphorylation", "08_Oxidative_phosphorylation", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Protein_RNA_Processing", "07_Protein_RNA_Processing", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Citrate_Cycle", "06_Citrate_Cycle", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Cofactor_vitamin_metabolism", "05_Cofactor_vitamin_metabolism", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Nucleotide_metabolism", "04_Nucleotide_metabolism", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Carbohydrate_Metabolism", "03_Carbohydrate_Metabolism", Ostreo$Annotation)
diel <- row.names(Ostreo)[which(Ostreo$Diel.Transcripts == TRUE)]
not.diel <- row.names(Ostreo)[which(Ostreo$Diel.Transcripts == FALSE)]

postscript(file="R_network_pipeline/Figures/Roseo.blc.dotplot.eps", height=4, width=9, horizontal = FALSE, onefile = FALSE)
ggplot() + geom_point(data=Ostreo[not.diel,], aes(Peak.Time, Annotation), shape=21, size=5, colour="grey60") + geom_rect(data=Ostreo, aes(xmin=-0.5, xmax=6.5, ymin=-Inf, ymax=+Inf), fill='grey90', colour=NA) + geom_rect(data=Ostreo, aes(xmin=18.5, xmax=24.5, ymin=-Inf, ymax=+Inf), fill='grey90', colour=NA) + geom_point(data=Ostreo[not.diel,], aes(Peak.Time, Annotation), shape=21, size=5, colour="grey60") + geom_point(shape=21, size=5.3) + geom_point(shape=21, size=5) + geom_point(data=Ostreo[diel,], aes(Peak.Time, Annotation), shape=21, size=5, colour="blue") + theme(axis.text.y = element_text(colour = 'grey40', size=10), axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(colour='grey60', size=12), panel.background = element_rect(fill='white', colour='white'), legend.position="none") + scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) + geom_vline(aes(xintercept=12), linetype='dashed', colour='grey70') + geom_vline(aes(xintercept=18), linetype='dashed', colour='grey70')
dev.off()

###############################################################

#SAR86
functions <- c("Amino_Acid_metabolism", "Carbohydrate_Metabolism", "Citrate_Cycle", "DNA/Chromosome_Processing", "Nucleotide_metabolism", "Oxidative_phosphorylation", "Protein_RNA_Processing", "Sulfur_metabolism", "Transcription", "Translation", "Transporters")
diel.t <- row.names(final.blc)[which(final.blc$Organism == "SAR86" & final.blc$Final_Annote %in% functions)]
Ostreo <- data.frame(cbind(final.blc[diel.t,]$Peak.Time, as.character(final.blc[diel.t,]$Final_Annote), as.character(final.blc[diel.t,]$diel.transcript)))
row.names(Ostreo) <- diel.t; colnames(Ostreo) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Ostreo$Peak.Time <- as.numeric(as.character(Ostreo$Peak.Time))

Ostreo$Annotation <- gsub("Transporters", "15_Transporters", Ostreo$Annotation)
Ostreo$Annotation <- gsub("DNA/Chromosome_Processing", "14_DNA/Chromosome_Processing", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Signal_Transduction", "13_Signal_Transduction", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Chemotaxis_and_Flagellar_Assembly", "12_Chemotaxis_and_Flagellar_Assembly", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Amino_Acid_metabolism", "11_Amino_Acid_metabolism", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Translation", "10_Translation", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Transcription", "09_Transcription", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Oxidative_phosphorylation", "08_Oxidative_phosphorylation", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Protein_RNA_Processing", "07_Protein_RNA_Processing", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Citrate_Cycle", "06_Citrate_Cycle", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Cofactor_vitamin_metabolism", "05_Cofactor_vitamin_metabolism", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Nucleotide_metabolism", "04_Nucleotide_metabolism", Ostreo$Annotation)
Ostreo$Annotation <- gsub("Carbohydrate_Metabolism", "03_Carbohydrate_Metabolism", Ostreo$Annotation)

diel <- row.names(Ostreo)[which(Ostreo$Diel.Transcripts == TRUE)]
not.diel <- row.names(Ostreo)[which(Ostreo$Diel.Transcripts == FALSE)]

postscript(file="R_network_pipeline/Figures/SAR86.blc.dotplot.eps", height=3.5, width=9, horizontal = FALSE, onefile = FALSE)
ggplot() + geom_point(data=Ostreo[not.diel,], aes(Peak.Time, Annotation), shape=21, size=5, colour="grey60") + geom_rect(data=Ostreo, aes(xmin=-0.5, xmax=6.5, ymin=-Inf, ymax=+Inf), fill='grey90', colour=NA) + geom_rect(data=Ostreo, aes(xmin=18.5, xmax=24.5, ymin=-Inf, ymax=+Inf), fill='grey90', colour=NA) + geom_point(data=Ostreo[not.diel,], aes(Peak.Time, Annotation), shape=21, size=5, colour="grey60") + geom_point(shape=21, size=5.3) + geom_point(shape=21, size=5) + geom_point(data=Ostreo[diel,], aes(Peak.Time, Annotation), shape=21, size=5, colour="blue") + theme(axis.text.y = element_text(colour = 'grey40', size=12), axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(colour='grey60', size=12), panel.background = element_rect(fill='white', colour='white'), legend.position="none") + scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) + geom_vline(aes(xintercept=12), linetype='dashed', colour='grey70') + geom_vline(aes(xintercept=18), linetype='dashed', colour='grey70')
dev.off()

########################################











            Amino_Acid_Metabolism       Aminoacyl_tRNA_Biosynthesis 
                  Carbon_Fixation                        Cell_Cycle 
        Central_Carbon_Metabolism Chemotaxis_and_Flagellar_Assembly 
                    Citrate_Cycle       Cofactor_vitamin_metabolism 
        DNA/Chromosome_Processing             Fatty_acid_metabolism 
              Glycan_biosynthesis                        Glycolysis 
                        Hydrolase                         Isomerase  
                           Ligase                  Lipid_Metabolism 
                            Lyase               Nitrogen_metabolism 
            Nucleotide_metabolism         Oxidative_phosphorylation 
                   Oxidoreductase                         Peptidase  
                   Photosynthesis            Protein_RNA_Processing 
                  Proteorhodopsin                    RNA_Polymerase 
             Secondary_metabolism                  Secretion_System 
              Signal_Transduction                 Sulfur_Metabolism 
                      Transferase                       Translation 
                     Transporters            Xenobiotic_degradation 











mod <- sort(table(meta$Module_Color), decreasing=T)

mypal = labels2colors(1:37)
ggplot(meta, aes(meta$Module, fill=meta$Final_Annote)) + geom_bar(stat="bin", position="fill") + coord_flip()+ scale_fill_manual(values=mypal)


num <- 37
result <- c()
pathways <- c("Amino_acid_metabolism",  "Carbohydrate_metabolism", "Chromosome_associated", "Aminoacyl_tRNA_Biosynthesis", "Translation", "Nucleotide_metabolism", "Carbon_fixation", "Photosynthesis", "Cell_Cycle", "Replication_Repair", "Central_Carbon_Metabolism", "Citrate_cycle", "Cofactor_vitamin_metabolism", "Flagellar_Assembly_and_Chemotaxis", "Glycan_biosynthesis", "Glycolysis", "Hydrolase", "Isomerase", "Ligase", "Lipid_metabolism", "Lyase", "Nitrogen_metabolism",  "Oxidative_phosphorylation", "Oxidoreductase", "Proteorhodopsin", "Peptidase",  "Protein_RNA_Processing", "Secondary_metabolism", "Secretion_system", "Signal_transduction", "Sulfur_metabolism", "Transcription", "Transferase",  "Transport_Catabolism", "Transporters", "Unassigned", "Xenobiotic_degradation")


