
##########################################################################################################
############################### Identify KEGG Annotations of highly-connected nodes ######################
##########################################################################################################

############################## Load Necessary Packages
library(WGCNA)
library(igraph)

############################ Enable multithreading for WGCNA
enableWGCNAThreads()

####################################################################################################### CANON:
############ Load raw count data and normalize by total counts in a given timepoint [column]
numsamples <- 35
canon.counts <- read.table(file="R_network_pipeline/CANON/canon.combined.counts.10_2014.txt", sep="\t", header=TRUE, row.names=1)
numgenes <- dim(canon.counts)[1]
Data.Counts <- canon.counts[1:numgenes,1:numsamples]
Total.Counts <- as.numeric(apply(canon.counts[1:numgenes,1:numsamples], 2, sum))
Norm=scale(Data.Counts, scale=Total.Counts, center=FALSE)
Norm_Final <- t(Norm)

#TOM = TOMsimilarityFromExpr(Norm_Final, power= pow)
cor <- (cor(Norm_Final))^5

############################################################################################
############## Load data and calculate cross taxa/pathway correlations for CANON ###########
############################################################################################

# CANON
canon <- read.table(file="R_network_pipeline/CANON/canon.final.combined.list2", sep="\t", header=TRUE, row.names=1, quote="")

functions <- c("Amino_Acid_metabolism", "DNA/Chromosome_Processing", "Translation", "Transporters", "Cofactor_vitamin_metabolism", "Carbohydrate_metabolism",   "Oxidative_phosphorylation", "Fatty_acid_metabolism", "Carbon_Fixation", "Transcription", "Photosynthesis")

functions <- c("Amino_Acid_metabolism", "Aminoacyl_tRNA_Biosynthesis", "Carbohydrate_metabolism", "Carbon_Fixation", "Cell_Cycle", "Central_Carbon_Metabolism", "Chemotaxis_and_Flagellar_Assembly", "Citrate_cycle", "Cofactor_vitamin_metabolism", "DNA/Chromosome_Processing", "Fatty_acid_metabolism", "Glycan_biosynthesis", "Glycolysis", "Hydrolase", "Isomerase", "Ligase", "Lipid_Metabolism", "Lyase", "Nitrogen_metabolism", "Nucleotide_metabolism", "Oxidative_phosphorylation", "Oxidoreductase", "Peptidase", "Photosynthesis", "Protein_RNA_Processing", "Proteorhodopsin", "Secondary_metabolism", "Secretion_System", "Signal_transduction", "Sulfur_metabolism", "Transcription", "Transferase", "Translation", "Transport_Catabolism", "Transporters", "Unassigned", "Xenobiotic_degradation")

taxa <- c("Ostreococcus", "SAR11_", "SAR116_", "SAR86", "Euryarchaota", "SAR92", "Flavobacteria", "Roseobacter", "ARCTIC", "SAR406")

labels=c("Ostreococcus.Amino_acid_metabolism", "Ostreococcus.Cofactor_vitamin_metabolism", "Ostreococcus.DNA/Chromosome_Processing", "Ostreococcus.Oxidative_phosphorylation", "Ostreococcus.Protein_RNA_Processing", "Ostreococcus.Transcription", "Ostreococcus.Translation", "Ostreococcus.Transporters", "SAR11_.Amino_acid_metabolism", "SAR11_.Central_Carbon_Metabolism", "SAR11_.Citrate_cycle", "SAR11_.Cofactor_vitamin_metabolism", "SAR11_.DNA/Chromosome_Processing", "SAR11_.Glycolysis", "SAR11_.Oxidative_phosphorylation", "SAR11_.Protein_RNA_Processing", "SAR11_.Transcription", "SAR11_.Translation", "SAR11_.Transporters", "SAR116_.Amino_acid_metabolism", "SAR116_.Central_Carbon_Metabolism", "SAR116_.Chemotaxis_and_Flagellar_Assembly", "SAR116_.Citrate_cycle", "SAR116_.Cofactor_vitamin_metabolism", "SAR116_.DNA/Chromosome_Processing", "SAR116_.Glycolysis", "SAR116_.Oxidative_phosphorylation", "SAR116_.Protein_RNA_Processing", "SAR116_.Translation", "SAR116_.Transporters", "SAR86.Amino_acid_metabolism", "SAR86.Central_Carbon_Metabolism", "SAR86.Citrate_cycle", "SAR86.Cofactor_vitamin_metabolism", "SAR86.DNA/Chromosome_Processing", "SAR86.Glycolysis", "SAR86.Oxidative_phosphorylation", "SAR86.Transcription", "SAR86.Translation", "SAR86.Transporters", "Roseobacter.Amino_acid_metabolism", "Roseobacter.Central_Carbon_Metabolism", "Roseobacter.Chemotaxis_and_Flagellar_Assembly", "Roseobacter.Citrate_cycle", "Roseobacter.Cofactor_vitamin_metabolism", "Roseobacter.Oxidative_phosphorylation", "Roseobacter.Protein_RNA_Processing", "Roseobacter.Transcription", "Roseobacter.Translation", "Roseobacter.Transporters", "SAR406.Amino_acid_metabolism", "SAR406.Citrate_cycle", "SAR406.Cofactor_vitamin_metabolism", "SAR406.Oxidative_phosphorylation", "SAR406.Protein_RNA_Processing", "SAR406.Transcription", "SAR406.Translation", "SAR406.Transporters")


list2 <- strsplit(labels, ".", fixed=TRUE)
canon.matrix <- matrix(NA, ncol=length(list2), nrow=length(list2))
colnames(canon.matrix) <- labels; row.names(canon.matrix) <- labels

for(i in 1:length(list2)) {
	c1 <- row.names(canon)[which(canon$Organism == list2[[i]][1] & canon$Final_Annote == list2[[i]][2])]
	for(j in 1:length(list2)) {
		c2 <- row.names(canon)[which(canon$Organism == list2[[j]][1] & canon$Final_Annote == list2[[j]][2])]
			if( list2[[j]][1] == list2[[i]][1]) { sum <- 0 }
			else { mat <- cor[c1, c2]; sum <- sum(mat) }
		correct <- 1000*(sum/(length(c1) * length(c2)))
		canon.matrix[i, j] <- correct
		#matrix[j, i] <- sum
	}
}

### Make heatmap
library(reshape2)
library(ggplot2)
library(grid)
dat <- melt(canon.matrix)
colnames(dat) <- c("cat1", "cat2", "val")

base_size=9
jpeg(file="CANON.crosscompare.TC.jpg", quality=100, res=600, height=11, width=20, units="in")
ggplot(dat, aes(cat1, cat2)) + geom_tile(aes(fill = val), colour = "white") + scale_fill_gradient2(low = "red", mid='white', high = "blue") + theme(legend.position = "none", axis.ticks = element_blank(), axis.text.x = element_text(size = base_size * 0.9, angle = 30, hjust = 1, vjust = 1, face = 'bold', colour = "black"), axis.text.y = element_text(face = 'bold', colour = "black")) + geom_abline(size=2, colour='grey80')
dev.off()

#write.table(matrix, file="canon.pairwise.comapre", sep="\t")

# BioLincs

############################################################################################
########## Load data and calculate cross taxa/pathway correlations for BioLincs ###########
###########################################################################################
#############################################################################################################################
############ Load raw count data and normalize by total counts in a given timepoint [column] [this is identical to script 1].


numsamples <- 30
blc.counts <- read.table(file="R_network_pipeline/BioLincs/all.biolincs.counts.10_2014.blc", sep="\t", header=TRUE, row.names=1)
numgenes <- dim(blc.counts)[1]
Data.Counts <- blc.counts[1:numgenes,1:numsamples]
Total.Counts <- as.numeric(apply(blc.counts[1:numgenes,1:numsamples], 2, sum))
Norm=scale(Data.Counts, scale=Total.Counts, center=FALSE)
Norm_Final <- t(Norm)

cor <- cor(Norm_Final)
#sign <- sign(cor)
#corexp <- cor^4
#cor <- corexp * sign
cor <- cor(Norm_Final)^5
#cor <- (cor(Norm_Final))^4

blc <- read.table(file="R_network_pipeline/BioLincs/blc.final.combined.list2", sep="\t", header=TRUE, row.names=1, quote="")
taxa <- c('Prochlorococcus', 'Roseo', 'SAR11_', 'SAR116_', 'SAR324', 'SAR406', 'SAR86')

labels=c("Prochlorococcus.Carbon_fixation", "Prochlorococcus.Cell_Cycle", "Prochlorococcus.Central_Carbon_Metabolism", "Prochlorococcus.Citrate_cycle", "Prochlorococcus.Cofactor_vitamin_metabolism", "Prochlorococcus.DNA/Chromosome_Processing", "Prochlorococcus.Oxidative_phosphorylation", "Prochlorococcus.Photosynthesis", "Prochlorococcus.Protein_RNA_Processing",  "Prochlorococcus.Transcription", "Prochlorococcus.Translation", "Prochlorococcus.Transporters", "SAR11_.Carbohydrate_metabolism", "SAR11_.Central_Carbon_Metabolism", "SAR11_.Citrate_cycle", "SAR11_.Cofactor_vitamin_metabolism", "SAR11_.DNA/Chromosome_Processing", "SAR11_.Glycolysis", "SAR11_.Oxidative_phosphorylation", "SAR11_.Protein_RNA_Processing",  "SAR11_.Transcription", "SAR11_.Translation", "SAR11_.Transporters", "SAR116_.Central_Carbon_Metabolism", "SAR116_.Chemotaxis_and_Flagellar_Assembly", "SAR116_.Citrate_cycle", "SAR116_.Cofactor_vitamin_metabolism", "SAR116_.DNA/Chromosome_Processing", "SAR116_.Glycolysis", "SAR116_.Oxidative_phosphorylation", "SAR116_.Protein_RNA_Processing", "SAR116_.Translation", "SAR116_.Transporters", "SAR86.Central_Carbon_Metabolism", "SAR86.Citrate_cycle", "SAR86.Cofactor_vitamin_metabolism", "SAR86.DNA/Chromosome_Processing", "SAR86.Glycolysis", "SAR86.Oxidative_phosphorylation", "SAR86.Transcription", "SAR86.Translation", "SAR86.Transporters", "Roseo.Central_Carbon_Metabolism", "Roseo.Chemotaxis_and_Flagellar_Assembly", "Roseo.Citrate_cycle", "Roseo.Cofactor_vitamin_metabolism", "Roseo.Oxidative_phosphorylation", "Roseo.Protein_RNA_Processing", "Roseo.Transcription", "Roseo.Translation", "Roseo.Transporters", "SAR406.Citrate_cycle", "SAR406.Cofactor_vitamin_metabolism", "SAR406.Oxidative_phosphorylation", "SAR406.Protein_RNA_Processing", "SAR406.Transcription", "SAR406.Translation", "SAR406.Transporters")

labels=c("Prochlorococcus.Amino_acid_metabolism", "Prochlorococcus.Cofactor_vitamin_metabolism", "Prochlorococcus.DNA/Chromosome_Processing", "Prochlorococcus.Oxidative_phosphorylation", "Prochlorococcus.Protein_RNA_Processing", "Prochlorococcus.Transcription", "Prochlorococcus.Translation", "Prochlorococcus.Transporters", "SAR11_.Amino_acid_metabolism", "SAR11_.Central_Carbon_Metabolism", "SAR11_.Citrate_cycle", "SAR11_.Cofactor_vitamin_metabolism", "SAR11_.DNA/Chromosome_Processing", "SAR11_.Glycolysis", "SAR11_.Oxidative_phosphorylation", "SAR11_.Protein_RNA_Processing", "SAR11_.Transcription", "SAR11_.Translation", "SAR11_.Transporters", "SAR116_.Amino_acid_metabolism", "SAR116_.Central_Carbon_Metabolism", "SAR116_.Chemotaxis_and_Flagellar_Assembly", "SAR116_.Citrate_cycle", "SAR116_.Cofactor_vitamin_metabolism", "SAR116_.DNA/Chromosome_Processing", "SAR116_.Glycolysis", "SAR116_.Oxidative_phosphorylation", "SAR116_.Protein_RNA_Processing", "SAR116_.Translation", "SAR116_.Transporters", "SAR86.Amino_acid_metabolism", "SAR86.Central_Carbon_Metabolism", "SAR86.Citrate_cycle", "SAR86.Cofactor_vitamin_metabolism", "SAR86.DNA/Chromosome_Processing", "SAR86.Glycolysis", "SAR86.Oxidative_phosphorylation", "SAR86.Transcription", "SAR86.Translation", "SAR86.Transporters", "Roseo.Amino_acid_metabolism", "Roseo.Central_Carbon_Metabolism", "Roseo.Chemotaxis_and_Flagellar_Assembly", "Roseo.Citrate_cycle", "Roseo.Cofactor_vitamin_metabolism", "Roseo.Oxidative_phosphorylation", "Roseo.Protein_RNA_Processing", "Roseo.Transcription", "Roseo.Translation", "Roseo.Transporters", "SAR406.Amino_acid_metabolism", "SAR406.Citrate_cycle", "SAR406.Cofactor_vitamin_metabolism", "SAR406.Oxidative_phosphorylation", "SAR406.Protein_RNA_Processing", "SAR406.Transcription", "SAR406.Translation", "SAR406.Transporters")
list2 <- strsplit(labels, ".", fixed=TRUE)

blc.matrix <- matrix(NA, ncol=length(labels), nrow=length(labels))
colnames(blc.matrix) <- labels; row.names(blc.matrix) <- labels

for(i in 1:length(list2)) {
	c1 <- row.names(blc)[which(blc$Organism == list2[[i]][1] & blc$Final_Annote == list2[[i]][2])]
	for(j in 1:length(list2)) {
		c2 <- row.names(blc)[which(blc$Organism == list2[[j]][1] & blc$Final_Annote == list2[[j]][2])]
			if( list2[[j]][1] == list2[[i]][1]) { sum <- 0 }
			else { mat <- cor[c1, c2]; sum <- sum(mat) }
		correct <- 1000*(sum/(length(c1) * length(c2)))
		blc.matrix[i, j] <- correct
		#matrix[j, i] <- sum
	}
}

### Make heatmap
library(reshape2)
library(ggplot2)
library(grid)
dat <- melt(blc.matrix)
colnames(dat) <- c("cat1", "cat2", "val")

base_size=9
jpeg(file="blc.crosscompare.jpg", quality=100, res=600, height=11, width=20, units="in")
ggplot(dat, aes(cat1, cat2)) + geom_tile(aes(fill = val), colour = "white") + scale_fill_gradient2(low = "red", mid='white', high = "blue") + theme(legend.position = "none", axis.ticks = element_blank(), axis.text.x = element_text(size = base_size * 0.9, angle = 30, hjust = 1, vjust = 1, face = 'bold', colour = "black"), axis.text.y = element_text(face = 'bold', colour = "black"))
dev.off()


#write.table(blc.matrix, file="blc.pairwise.comapre", sep="\t")


######################## For combined heatmap
matrix3 <- matrix(NA, ncol=length(labels), nrow=length(labels))
for(i in 1:dim(blc.matrix)[1]) {
	for(j in 1:dim(blc.matrix)[1]) {
		matrix3[i, j] <- blc.matrix[i, j]
		matrix3[j, i] <- canon.matrix[i, j]
	}
}
library(reshape2)
library(ggplot2)
dat <- melt(matrix3)
colnames(dat) <- c("cat1", "cat2", "val")

base_size=9
jpeg(file="combined.crosscompare3.jpg", quality=100, res=600, height=17, width=20, units="in")
ggplot(dat, aes(cat1, cat2)) + geom_tile(aes(fill = val), colour = "white") + scale_fill_gradient2(low = "red", mid='white', high = "blue") + theme(legend.position = "none", axis.ticks = element_blank(), axis.text.x = element_text(size = base_size * 0.9, angle = 30, hjust = 1, vjust = 1, face = 'bold', colour = "black"), axis.text.y = element_text(face = 'bold', colour = "black")) + geom_abline(size=0.5, colour='black')
dev.off()

########################################################


########################################################################################################################
################ This is a re-do using taxa-specific normalizations prior to correlation calculations ##################
########################################################################################################################

####################################################################################################### CANON:
############ Load raw count data and normalize by total counts in a given timepoint [column]

########
canon <- read.table(file="R_network_pipeline/CANON/canon.combined.summary2", sep="\t", quote="", header=TRUE, row.names=1)
list <- c("Ostreococcus", "SAR11_", "SAR116_", "SAR86", "Roseobacter", "SAR406", "Euryarchaota", "SAR92", "Flavobacteria", "ARCTIC")
result <- c()
taxa.mat <- c()
for(i in 1:length(list)) {
taxa.spec <- canon[(canon$Organism == list[i]),25:59]
counts <- apply(taxa.spec, 2, sum)
mat <- scale(taxa.spec, scale=counts, center=F)
taxa.mat <- data.frame(cbind(taxa.mat, t(mat)))
}
cor <- (cor(taxa.mat))^5

labels=c("Ostreococcus.Carbon_Fixation", "Ostreococcus.Cell_Cycle", "Ostreococcus.Central_Carbon_Metabolism", "Ostreococcus.Citrate_cycle", "Ostreococcus.Cofactor_vitamin_metabolism", "Ostreococcus.DNA/Chromosome_Processing", "Ostreococcus.Glycolysis", "Ostreococcus.Oxidative_phosphorylation", "Ostreococcus.Photosynthesis", "Ostreococcus.Protein_RNA_Processing", "Ostreococcus.Signal_transduction", "Ostreococcus.Transcription", "Ostreococcus.Translation", "Ostreococcus.Transporters", "SAR11_.Central_Carbon_Metabolism", "SAR11_.Citrate_cycle", "SAR11_.Cofactor_vitamin_metabolism", "SAR11_.DNA/Chromosome_Processing", "SAR11_.Glycolysis", "SAR11_.Oxidative_phosphorylation", "SAR11_.Protein_RNA_Processing", "SAR11_.Signal_transduction", "SAR11_.Transcription", "SAR11_.Translation", "SAR11_.Transporters", "SAR116_.Central_Carbon_Metabolism", "SAR116_.Chemotaxis_and_Flagellar_Assembly", "SAR116_.Citrate_cycle", "SAR116_.Cofactor_vitamin_metabolism", "SAR116_.DNA/Chromosome_Processing", "SAR116_.Glycolysis", "SAR116_.Oxidative_phosphorylation", "SAR116_.Protein_RNA_Processing", "SAR116_.Signal_transduction", "SAR116_.Translation", "SAR116_.Transporters", "SAR86.Central_Carbon_Metabolism", "SAR86.Citrate_cycle", "SAR86.Cofactor_vitamin_metabolism", "SAR86.DNA/Chromosome_Processing", "SAR86.Glycolysis", "SAR86.Oxidative_phosphorylation", "SAR86.Signal_transduction", "SAR86.Transcription", "SAR86.Translation", "SAR86.Transporters", "Roseobacter.Central_Carbon_Metabolism", "Roseobacter.Chemotaxis_and_Flagellar_Assembly", "Roseobacter.Citrate_cycle", "Roseobacter.Cofactor_vitamin_metabolism", "Roseobacter.Oxidative_phosphorylation", "Roseobacter.Protein_RNA_Processing", "Roseobacter.Transcription", "Roseobacter.Translation", "Roseobacter.Transporters", "SAR406.Citrate_cycle", "SAR406.Cofactor_vitamin_metabolism", "SAR406.Oxidative_phosphorylation", "SAR406.Protein_RNA_Processing", "SAR406.Transcription", "SAR406.Translation", "SAR406.Transporters")

list2 <- strsplit(labels, ".", fixed=TRUE)

matrix2 <- matrix(NA, ncol=length(list2), nrow=length(list2))
colnames(matrix2) <- labels; row.names(matrix2) <- labels

for(i in 1:length(list2)) {
	c1 <- row.names(canon)[which(canon$Organism == list2[[i]][1] & canon$Final_Annote == list2[[i]][2])]
	for(j in 1:length(list2)) {
		c2 <- row.names(canon)[which(canon$Organism == list2[[j]][1] & canon$Final_Annote == list2[[j]][2])]
			if( list2[[j]][1] == list2[[i]][1]) { sum <- 0 }
			else { mat <- cor[c1, c2]; sum <- sum(mat) }
		correct <- 1000*(sum/(length(c1) * length(c2)))
		matrix2[i, j] <- correct
		#matrix[j, i] <- sum
	}
}


library(reshape2)
library(ggplot2)
library(grid)
dat <- melt(matrix2)
colnames(dat) <- c("cat1", "cat2", "val")

base_size=9
jpeg(file="CANON.crosscompare.new.no.scale.jpg", quality=100, res=600, height=11, width=20, units="in")
ggplot(dat, aes(cat1, cat2)) + geom_tile(aes(fill = val), colour = "white") + scale_fill_gradient2(low = "red", mid='white', high = "blue") + theme(legend.position = "none", axis.ticks = element_blank(), axis.text.x = element_text(size = base_size * 0.9, angle = 30, hjust = 1, vjust = 1, face = 'bold', colour = "black"), axis.text.y = element_text(face = 'bold', colour = "black"))
dev.off()

base_size=9
jpeg(file="CANON.crosscompare.new.no.scale.jpg", quality=100, res=600, height=11, width=20, units="in")
ggplot(dat, aes(cat1, cat2)) + geom_tile(aes(fill = val), colour = "white") + scale_fill_gradient2(low = "red", mid='white', high = "blue") + theme(legend.position = "none", axis.ticks = element_blank(), axis.text.x = element_text(size = base_size * 0.9, angle = 30, hjust = 1, vjust = 1, face = 'bold', colour = "black"), axis.text.y = element_text(face = 'bold', colour = "black"))
dev.off()

########################## New Normalizaton for BioLincs


blc <- read.table(file="R_network_pipeline/BioLincs/blc.combined.summary", sep="\t", header=TRUE, row.names=1, quote="")
list <- c('Pro', 'Roseo', 'SAR11_', 'SAR116_', 'SAR406', 'SAR86')
result <- c()
taxa.mat <- c()
for(i in 1:length(list)) {
taxa.spec <- blc[(blc$Organism == list[i]),18:47]
counts <- apply(taxa.spec, 2, sum)
mat <- scale(taxa.spec, scale=counts, center=F)
taxa.mat <- data.frame(cbind(taxa.mat, t(mat)))
}
cor <- (cor(taxa.mat))^5

labels=c("Pro.Carbon_Fixation", "Pro.Citrate_cycle", "Pro.Cofactor_vitamin_metabolism", "Pro.DNA/Chromosome_Processing", "Pro.Glycolysis", "Pro.Oxidative_phosphorylation", "Pro.Photosynthesis", "Pro.Protein_RNA_Processing", "Pro.Signal_transduction", "Pro.Transcription", "Pro.Translation", "Pro.Transporters", "SAR11_.Central_Carbon_Metabolism", "SAR11_.Citrate_cycle", "SAR11_.Cofactor_vitamin_metabolism", "SAR11_.DNA/Chromosome_Processing", "SAR11_.Glycolysis", "SAR11_.Oxidative_phosphorylation", "SAR11_.Protein_RNA_Processing", "SAR11_.Signal_transduction", "SAR11_.Transcription", "SAR11_.Translation", "SAR11_.Transporters", "SAR116_.Central_Carbon_Metabolism", "SAR116_.Chemotaxis_and_Flagellar_Assembly", "SAR116_.Citrate_cycle", "SAR116_.Cofactor_vitamin_metabolism", "SAR116_.DNA/Chromosome_Processing", "SAR116_.Glycolysis", "SAR116_.Oxidative_phosphorylation", "SAR116_.Protein_RNA_Processing", "SAR116_.Signal_transduction", "SAR116_.Translation", "SAR116_.Transporters", "SAR86.Central_Carbon_Metabolism", "SAR86.Citrate_cycle", "SAR86.Cofactor_vitamin_metabolism", "SAR86.DNA/Chromosome_Processing", "SAR86.Glycolysis", "SAR86.Oxidative_phosphorylation", "SAR86.Signal_transduction", "SAR86.Transcription", "SAR86.Translation", "SAR86.Transporters", "Roseo.Central_Carbon_Metabolism", "Roseo.Chemotaxis_and_Flagellar_Assembly", "Roseo.Citrate_cycle", "Roseo.Cofactor_vitamin_metabolism", "Roseo.Oxidative_phosphorylation", "Roseo.Protein_RNA_Processing", "Roseo.Transcription", "Roseo.Translation", "Roseo.Transporters", "SAR406.Citrate_cycle", "SAR406.Cofactor_vitamin_metabolism", "SAR406.Oxidative_phosphorylation", "SAR406.Protein_RNA_Processing", "SAR406.Transcription", "SAR406.Translation", "SAR406.Transporters")

list2 <- strsplit(labels, ".", fixed=TRUE)

matrix <- matrix(NA, ncol=length(labels), nrow=length(labels))
colnames(matrix) <- labels; row.names(matrix) <- labels

for(i in 1:length(list2)) {
	c1 <- row.names(blc)[which(blc$Organism == list2[[i]][1] & blc$Final_Annote == list2[[i]][2])]
	for(j in 1:length(list2)) {
		c2 <- row.names(blc)[which(blc$Organism == list2[[j]][1] & blc$Final_Annote == list2[[j]][2])]
			if( list2[[j]][1] == list2[[i]][1]) { sum <- 0 }
			else { mat <- cor[c1, c2]; sum <- sum(mat) }
		correct <- 1000*(sum/(length(c1) * length(c2)))
		matrix[i, j] <- correct
		#matrix[j, i] <- sum
	}
}

### Make heatmap
library(reshape2)
library(ggplot2)
library(grid)
dat2 <- melt(matrix)
colnames(dat2) <- c("cat1", "cat2", "val")

base_size=9
jpeg(file="blc.crosscompare3.jpg", quality=100, res=600, height=11, width=20, units="in")
ggplot(dat2, aes(cat1, cat2)) + geom_tile(aes(fill = val), colour = "white") + scale_fill_gradient2(low = "red", mid='white', high = "blue") + theme(legend.position = "none", axis.ticks = element_blank(), axis.text.x = element_text(size = base_size * 0.9, angle = 30, hjust = 1, vjust = 1, face = 'bold', colour = "black"), axis.text.y = element_text(face = 'bold', colour = "black"))
dev.off()











