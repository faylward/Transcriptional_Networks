
#######################################################################################
################################ Calculate Eigengenes #################################
#######################################################################################

##################################################################################################################################################
############################## Load Necessary Packages
library(ggplot2)

########## Load big summary files CANON

diel <- read.table(file="R_network_pipeline/CANON/All.canon.hra.out", sep="\t", header=TRUE, row.names=1)
diel.transcript <- (diel$Perm.FDR <= 0.1 & diel$Regression.FDR <= 0.1)
diel.2 <- cbind(diel, diel.transcript)

canon <- read.table(file="R_network_pipeline/CANON/canon.combined.summary2", sep="\t", header=TRUE, row.names=1)
final.canon <- merge(diel.2, canon, by="row.names", all=T)
row.names(final.canon) <- final.canon$Row.names

length(row.names(final.canon)[which(final.canon$diel.transcript == TRUE & final.canon$Organism == 'SAR116_')])

functions <- c("Amino_Acid_metabolism", "Aminoacyl_tRNA_Biosynthesis", "Carbohydrate_Metabolism", "Carbon_Fixation", "Cell_Cycle", "Central_Carbon_Metabolism", "Chemotaxis_and_Flagellar_Assembly", "Citrate_Cycle", "Cofactor_vitamin_metabolism", "DNA/Chromosome_Processing", "Fatty_acid_metabolism", "Glycan_biosynthesis", "Glycolysis", "Hydrolase", "Isomerase", "Ligase", "Lipid_Metabolism", "Lyase", "Nitrogen_metabolism", "Nucleotide_metabolism", "Oxidative_phosphorylation", "Oxidoreductase", "Peptidase", "Photosynthesis", "Protein_RNA_Processing", "Proteorhodopsin", "Secondary_metabolism", "Secretion_System", "Signal_Transduction", "Sulfur_metabolism", "Transcription", "Transferase", "Translation", "Transport_Catabolism", "Transporters", "Unassigned", "Xenobiotic_degradation")

#Ostreococcus

functions <- c("Amino_Acid_metabolism", "Carbohydrate_Metabolism", "Carbon_Fixation", "Cell_Cycle", "Cofactor_vitamin_metabolism", "DNA/Chromosome_Processing", "Fatty_acid_metabolism", "Nucleotide_metabolism", "Oxidative_phosphorylation", "Photosynthesis", "Protein_RNA_Processing", "Signal_Transduction", "Transcription", "Translation", "Transporters")
diel.t <- row.names(final.canon)[which(final.canon$Organism == "Ostreococcus" & final.canon$Final_Annote %in% functions)]
Ostreo <- data.frame(cbind(final.canon[diel.t,]$Peak.Time, as.character(final.canon[diel.t,]$Final_Annote), as.character(final.canon[diel.t,]$diel.transcript)))
row.names(Ostreo) <- diel.t; colnames(Ostreo) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Ostreo$Peak.Time <- as.numeric(as.character(Ostreo$Peak.Time))
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

#postscript(file="R_network_pipeline/Figures/Ostreo.dotplot.eps", height=4.5, width=9, horizontal = FALSE, onefile = FALSE)
postscript("R_network_pipeline/Figures/Ostreo.dotplot.eps", height=4.5, width=9, horizontal = FALSE, onefile = FALSE)
ggplot() + geom_point(data=Ostreo[not.diel,], aes(Peak.Time, Annotation), shape=21, size=5, colour="grey60") + geom_rect(data=Ostreo, aes(xmin=-0.5, xmax=6.5, ymin=-Inf, ymax=+Inf), fill='grey90', colour=NA) + geom_rect(data=Ostreo, aes(xmin=19.5, xmax=24.5, ymin=-Inf, ymax=+Inf), fill='grey90', colour=NA) + geom_point(data=Ostreo[not.diel,], aes(Peak.Time, Annotation), shape=21, size=5, colour="grey60") + geom_point(shape=21, size=5.3) + geom_point(shape=21, size=5) + geom_point(data=Ostreo[diel,], aes(Peak.Time, Annotation), shape=21, size=5, colour="blue") + theme(axis.text.y = element_text(colour = 'grey40', size=10, face='bold'), axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(colour='grey60', size=12), panel.background = element_rect(fill='white', colour='white'), legend.position="none") + scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) + geom_vline(aes(xintercept=12), linetype='dashed', colour='grey70') + geom_vline(aes(xintercept=18), linetype='dashed', colour='grey70')
dev.off()

#SAR116
functions <- c("Amino_Acid_metabolism", "Carbohydrate_Metabolism", "Chemotaxis_and_Flagellar_Assembly", "Citrate_Cycle", "Cofactor_vitamin_metabolism", "DNA/Chromosome_Processing", "Nucleotide_metabolism", "Oxidative_phosphorylation", "Protein_RNA_Processing", "Signal_Transduction", "Sulfur_metabolism", "Transcription", "Translation", "Transporters")
diel.t <- row.names(final.canon)[which(final.canon$Organism == "SAR116_" & final.canon$Final_Annote %in% functions)]
Ostreo <- data.frame(cbind(final.canon[diel.t,]$Peak.Time, as.character(final.canon[diel.t,]$Final_Annote), as.character(final.canon[diel.t,]$diel.transcript)))
row.names(Ostreo) <- diel.t; colnames(Ostreo) <- c("Peak.Time", "Annotation", "Diel.Transcripts")
Ostreo$Peak.Time <- as.numeric(as.character(Ostreo$Peak.Time))
Ostreo$Peak.Time <- as.numeric(as.character(Ostreo$Peak.Time))
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

postscript(file="R_network_pipeline/Figures/SAR116.canon.dotplot.eps", height=4.5, width=9, horizontal = FALSE, onefile = FALSE)
ggplot() + geom_point(data=Ostreo[not.diel,], aes(Peak.Time, Annotation), shape=21, size=5, colour="grey60") + geom_rect(data=Ostreo, aes(xmin=-0.5, xmax=6.5, ymin=-Inf, ymax=+Inf), fill='grey90', colour=NA) + geom_rect(data=Ostreo, aes(xmin=19.5, xmax=24.5, ymin=-Inf, ymax=+Inf), fill='grey90', colour=NA) + geom_point(data=Ostreo[not.diel,], aes(Peak.Time, Annotation), shape=21, size=5, colour="grey60") + geom_point(shape=21, size=5.3) + geom_point(shape=21, size=5) + geom_point(data=Ostreo[diel,], aes(Peak.Time, Annotation), shape=21, size=5, colour="blue") + theme(axis.text.y = element_text(colour = 'grey40', size=10), axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(colour='grey60', size=12), panel.background = element_rect(fill='white', colour='white'), legend.position="none") + scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) + geom_vline(aes(xintercept=12), linetype='dashed', colour='grey70') + geom_vline(aes(xintercept=18), linetype='dashed', colour='grey70')
dev.off()

#SAR11
functions <- c("Amino_Acid_metabolism", "Carbohydrate_Metabolism", "Citrate_Cycle", "Cofactor_vitamin_metabolism", "DNA/Chromosome_Processing", "Fatty_acid_metabolism", "Nucleotide_metabolism", "Oxidative_phosphorylation", "Protein_RNA_Processing", "Signal_Transduction", "Transcription", "Translation", "Transporters")
diel.t <- row.names(final.canon)[which(final.canon$Organism == "SAR11_" & final.canon$Final_Annote %in% functions)]
Ostreo <- data.frame(cbind(final.canon[diel.t,]$Peak.Time, as.character(final.canon[diel.t,]$Final_Annote), as.character(final.canon[diel.t,]$diel.transcript)))
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

postscript(file="R_network_pipeline/Figures/SAR11.canon.dotplot.eps", height=4, width=9, horizontal = FALSE, onefile = FALSE)
ggplot() + geom_point(data=Ostreo[not.diel,], aes(Peak.Time, Annotation), shape=21, size=5, colour="grey60") + geom_rect(data=Ostreo, aes(xmin=-0.5, xmax=6.5, ymin=-Inf, ymax=+Inf), fill='grey90', colour=NA) + geom_rect(data=Ostreo, aes(xmin=19.5, xmax=24.5, ymin=-Inf, ymax=+Inf), fill='grey90', colour=NA) + geom_point(data=Ostreo[not.diel,], aes(Peak.Time, Annotation), shape=21, size=5, colour="grey60") + geom_point(shape=21, size=5.3) + geom_point(shape=21, size=5) + geom_point(data=Ostreo[diel,], aes(Peak.Time, Annotation), shape=21, size=5, colour="blue") + theme(axis.text.y = element_text(colour = 'grey40', size=10), axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(colour='grey60', size=12), panel.background = element_rect(fill='white', colour='white'), legend.position="none") + scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) + geom_vline(aes(xintercept=12), linetype='dashed', colour='grey70') + geom_vline(aes(xintercept=18), linetype='dashed', colour='grey70')
dev.off()

#Roseo
functions <- c("Amino_Acid_metabolism", "Carbohydrate_Metabolism", "Chemotaxis_and_Flagellar_Assembly", "Citrate_Cycle", "Cofactor_vitamin_metabolism", "DNA/Chromosome_Processing", "Nucleotide_metabolism", "Oxidative_phosphorylation", "Protein_RNA_Processing", "Signal_Transduction", "Sulfur_metabolism", "Transcription", "Translation", "Transporters")
diel.t <- row.names(final.canon)[which(final.canon$Organism == "Roseobacter" & final.canon$Final_Annote %in% functions)]
Ostreo <- data.frame(cbind(final.canon[diel.t,]$Peak.Time, as.character(final.canon[diel.t,]$Final_Annote), as.character(final.canon[diel.t,]$diel.transcript)))
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

postscript(file="R_network_pipeline/Figures/Roseo.canon.dotplot.eps", height=4, width=9, horizontal = FALSE, onefile = FALSE)
ggplot() + geom_point(data=Ostreo[not.diel,], aes(Peak.Time, Annotation), shape=21, size=5, colour="grey60") + geom_rect(data=Ostreo, aes(xmin=-0.5, xmax=6.5, ymin=-Inf, ymax=+Inf), fill='grey90', colour=NA) + geom_rect(data=Ostreo, aes(xmin=19.5, xmax=24.5, ymin=-Inf, ymax=+Inf), fill='grey90', colour=NA) + geom_point(data=Ostreo[not.diel,], aes(Peak.Time, Annotation), shape=21, size=5, colour="grey60") + geom_point(shape=21, size=5.3) + geom_point(shape=21, size=5) + geom_point(data=Ostreo[diel,], aes(Peak.Time, Annotation), shape=21, size=5, colour="blue") + theme(axis.text.y = element_text(colour = 'grey40', size=10), axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(colour='grey60', size=12), panel.background = element_rect(fill='white', colour='white'), legend.position="none") + scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) + geom_vline(aes(xintercept=12), linetype='dashed', colour='grey70') + geom_vline(aes(xintercept=18), linetype='dashed', colour='grey70')
dev.off()

#SAR86
functions <- c("Amino_Acid_metabolism", "Carbohydrate_Metabolism", "Citrate_Cycle", "DNA/Chromosome_Processing", "Nucleotide_metabolism", "Oxidative_phosphorylation", "Protein_RNA_Processing", "Sulfur_metabolism", "Transcription", "Translation", "Transporters")
diel.t <- row.names(final.canon)[which(final.canon$Organism == "SAR86" & final.canon$Final_Annote %in% functions)]
Ostreo <- data.frame(cbind(final.canon[diel.t,]$Peak.Time, as.character(final.canon[diel.t,]$Final_Annote), as.character(final.canon[diel.t,]$diel.transcript)))
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

postscript(file="R_network_pipeline/Figures/SAR86.canon.dotplot.eps", height=3.5, width=9, horizontal = FALSE, onefile = FALSE)
ggplot() + geom_point(data=Ostreo[not.diel,], aes(Peak.Time, Annotation), shape=21, size=5, colour="grey60") + geom_rect(data=Ostreo, aes(xmin=-0.5, xmax=6.5, ymin=-Inf, ymax=+Inf), fill='grey90', colour=NA) + geom_rect(data=Ostreo, aes(xmin=19.5, xmax=24.5, ymin=-Inf, ymax=+Inf), fill='grey90', colour=NA) + geom_point(data=Ostreo[not.diel,], aes(Peak.Time, Annotation), shape=21, size=5, colour="grey60") + geom_point(shape=21, size=5.3) + geom_point(shape=21, size=5) + geom_point(data=Ostreo[diel,], aes(Peak.Time, Annotation), shape=21, size=5, colour="blue") + theme(axis.text.y = element_text(colour = 'grey40', size=10), axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x = element_text(colour='grey60', size=12), panel.background = element_rect(fill='white', colour='white'), legend.position="none") + scale_x_continuous(breaks=c(0, 6, 12, 18, 24)) + geom_vline(aes(xintercept=12), linetype='dashed', colour='grey70') + geom_vline(aes(xintercept=18), linetype='dashed', colour='grey70')
dev.off()

########################################################################
########## Load big summary files BioLincs #############################
########################################################################

diel <- read.table(file="R_network_pipeline/BioLincs/All.blc.hra", sep="\t", header=TRUE, row.names=1)
diel.transcript <- (diel$Perm.FDR <= 0.1 & diel$Regression.FDR <= 0.1)
diel.2 <- cbind(diel, diel.transcript)

blc <- read.table(file="R_network_pipeline/BioLincs/blc.combined.summary2", sep="\t", header=TRUE, row.names=1, quote="")
final.blc <- merge(diel.2, blc, by="row.names", all=T)
row.names(final.blc) <- final.blc$Row.names

length(row.names(final.blc)[which(final.blc$diel.transcript == TRUE & final.blc$Organism == 'SAR116_')])

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


