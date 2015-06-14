
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
canon.counts <- read.table(file="R_network_pipeline/CANON/Combined_Canon", sep="\t", header=TRUE, row.names=1)
numgenes <- dim(canon.counts)[1]
Data.Counts <- canon.counts[1:numgenes,1:numsamples]
Total.Counts <- as.numeric(apply(canon.counts[1:numgenes,1:numsamples], 2, sum))
Norm=scale(Data.Counts, scale=Total.Counts, center=FALSE)
Norm_Final <- t(Norm)
cor <- (cor(Norm_Final))^5

############################################################################################
############## Load data and calculate cross taxa/pathway correlations for CANON ###########
###########################################################################################

# CANON
canon <- read.table(file="R_network_pipeline/CANON/canon.combined.summary2", sep="\t", header=TRUE, row.names=1)
functions <- c("Amino_Acid_metabolism", "DNA/Chromosome_Processing", "Translation", "Transporters", "Cofactor_vitamin_metabolism", "Carbohydrate_Metabolism",   "Oxidative_phosphorylation", "Fatty_acid_metabolism", "Carbon_Fixation", "Transcription", "Photosynthesis")

functions <- c("Amino_Acid_metabolism", "Aminoacyl_tRNA_Biosynthesis", "Carbohydrate_Metabolism", "Carbon_Fixation", "Cell_Cycle", "Central_Carbon_Metabolism", "Chemotaxis_and_Flagellar_Assembly", "Citrate_Cycle", "Cofactor_vitamin_metabolism", "DNA/Chromosome_Processing", "Fatty_acid_metabolism", "Glycan_biosynthesis", "Glycolysis", "Hydrolase", "Isomerase", "Ligase", "Lipid_Metabolism", "Lyase", "Nitrogen_metabolism", "Nucleotide_metabolism", "Oxidative_phosphorylation", "Oxidoreductase", "Peptidase", "Photosynthesis", "Protein_RNA_Processing", "Proteorhodopsin", "Secondary_metabolism", "Secretion_System", "Signal_Transduction", "Sulfur_metabolism", "Transcription", "Transferase", "Translation", "Transport_Catabolism", "Transporters", "Unassigned", "Xenobiotic_degradation")


taxa <- c("Ostreococcus", "SAR11_", "SAR116_", "SAR86", "Euryarchaota", "SAR92", "Flavobacteria", "Roseobacter", "ARCTIC", "SAR406")

list <- c()
h <- 0
for(i in 1:length(taxa)) {
	for(j in 1:length(functions)) {
		name <- paste(taxa[i], functions[j], sep=".")
		h <- h + 1
		list[h] <- name
	}
}
list2 <- strsplit(list, ".", fixed=TRUE)

labels=c("Ostreococcus.Amino_Acid_metabolism", "Ostreococcus.Carbohydrate_Metabolism", "Ostreococcus.Carbon_Fixation", "Ostreococcus.Cell_Cycle", "Ostreococcus.Central_Carbon_Metabolism", "Ostreococcus.Cofactor_vitamin_metabolism", "Ostreococcus.DNA/Chromosome_Processing", "Ostreococcus.Oxidative_phosphorylation", "Ostreococcus.Oxidoreductase", "Ostreococcus.Photosynthesis", "Ostreococcus.Translation", "Ostreococcus.Transporters", "SAR11_.Amino_Acid_metabolism", "SAR11_.Carbohydrate_Metabolism", "SAR11_.Cofactor_vitamin_metabolism", "SAR11_.DNA/Chromosome_Processing", "SAR11_.Oxidative_phosphorylation", "SAR11_.Oxidoreductase", "SAR11_.Translation", "SAR11_.Transporters", "SAR116_.Amino_Acid_metabolism", "SAR116_.Carbohydrate_Metabolism", "SAR116_.Central_Carbon_Metabolism", "SAR116_.Chemotaxis_and_Flagellar_Assembly", "SAR116_.Cofactor_vitamin_metabolism", "SAR116_.DNA/Chromosome_Processing", "SAR116_.Oxidative_phosphorylation", "SAR116_.Oxidoreductase", "SAR116_.Translation", "SAR116_.Transporters", "SAR86.Amino_Acid_metabolism", "SAR86.Carbohydrate_Metabolism", "SAR86.Central_Carbon_Metabolism", "SAR86.Cofactor_vitamin_metabolism", "SAR86.DNA/Chromosome_Processing", "SAR86.Oxidative_phosphorylation", "SAR86.Oxidoreductase", "SAR86.Translation", "SAR86.Transporters", "Roseobacter.Amino_Acid_metabolism", "Roseobacter.Carbohydrate_Metabolism", "Roseobacter.Central_Carbon_Metabolism", "Roseobacter.Chemotaxis_and_Flagellar_Assembly", "Roseobacter.Cofactor_vitamin_metabolism", "Roseobacter.DNA/Chromosome_Processing", "Roseobacter.Oxidative_phosphorylation", "Roseobacter.Oxidoreductase", "Roseobacter.Translation", "Roseobacter.Transporters")

labels=c("Ostreococcus.Carbon_Fixation", "Ostreococcus.Cell_Cycle", "Ostreococcus.Central_Carbon_Metabolism", "Ostreococcus.Citrate_Cycle", "Ostreococcus.Cofactor_vitamin_metabolism", "Ostreococcus.DNA/Chromosome_Processing", "Ostreococcus.Glycolysis", "Ostreococcus.Oxidative_phosphorylation", "Ostreococcus.Photosynthesis", "Ostreococcus.Protein_RNA_Processing", "Ostreococcus.Signal_Transduction", "Ostreococcus.Transcription", "Ostreococcus.Translation", "Ostreococcus.Transporters", "SAR11_.Carbohydrate_Metabolism", "SAR11_.Central_Carbon_Metabolism", "SAR11_.Citrate_Cycle", "SAR11_.Cofactor_vitamin_metabolism", "SAR11_.DNA/Chromosome_Processing", "SAR11_.Glycolysis", "SAR11_.Oxidative_phosphorylation", "SAR11_.Protein_RNA_Processing", "SAR11_.Signal_Transduction", "SAR11_.Transcription", "SAR11_.Translation", "SAR11_.Transporters", "SAR116_.Central_Carbon_Metabolism", "SAR116_.Chemotaxis_and_Flagellar_Assembly", "SAR116_.Citrate_Cycle", "SAR116_.Cofactor_vitamin_metabolism", "SAR116_.DNA/Chromosome_Processing", "SAR116_.Glycolysis", "SAR116_.Oxidative_phosphorylation", "SAR116_.Protein_RNA_Processing", "SAR116_.Signal_Transduction", "SAR116_.Translation", "SAR116_.Transporters", "SAR86.Central_Carbon_Metabolism", "SAR86.Citrate_Cycle", "SAR86.Cofactor_vitamin_metabolism", "SAR86.DNA/Chromosome_Processing", "SAR86.Glycolysis", "SAR86.Oxidative_phosphorylation", "SAR86.Signal_Transduction", "SAR86.Transcription", "SAR86.Translation", "SAR86.Transporters", "Roseobacter.Central_Carbon_Metabolism", "Roseobacter.Chemotaxis_and_Flagellar_Assembly", "Roseobacter.Citrate_Cycle", "Roseobacter.Cofactor_vitamin_metabolism", "Roseobacter.Oxidative_phosphorylation", "Roseobacter.Protein_RNA_Processing", "Roseobacter.Transcription", "Roseobacter.Translation", "Roseobacter.Transporters", "SAR406.Citrate_Cycle", "SAR406.Cofactor_vitamin_metabolism", "SAR406.Oxidative_phosphorylation", "SAR406.Protein_RNA_Processing", "SAR406.Transcription", "SAR406.Translation", "SAR406.Transporters")

list2 <- strsplit(labels, ".", fixed=TRUE)

matrix2 <- matrix(NA, ncol=length(list2), nrow=length(list2))
colnames(matrix2) <- labels; row.names(matrix) <- labels

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

### Make heatmap
library(reshape2)
library(ggplot2)
library(grid)
dat <- melt(matrix2)
colnames(dat) <- c("cat1", "cat2", "val")

base_size=9
jpeg(file="CANON.crosscompare.jpg", quality=100, res=600, height=11, width=20, units="in")
ggplot(dat, aes(cat1, cat2)) + geom_tile(aes(fill = val), colour = "white") + scale_fill_gradient2(low = "red", mid='white', high = "blue") + theme(legend.position = "none", axis.ticks = element_blank(), axis.text.x = element_text(size = base_size * 0.9, angle = 30, hjust = 1, vjust = 1, face = 'bold', colour = "black"), axis.text.y = element_text(face = 'bold', colour = "black")) + geom_abline(size=2, colour='grey80')
dev.off()


jpeg(file="canon.crosscompare4.jpg", quality=100, res=600, height=11, width=20, units="in")
 heatmap(matrix, Rowv=NA, Colv=NA)
dev.off()


write.table(matrix, file="canon.pairwise.comapre", sep="\t")

# BioLincs

############################################################################################
########## Load data and calculate cross taxa/pathway correlations for BioLincs ###########
###########################################################################################
#############################################################################################################################
############ Load raw count data and normalize by total counts in a given timepoint [column] [this is identical to script 1].
numsamples <- 30
blc.counts <- read.table(file="R_network_pipeline/BioLincs/Biolincs.counts.combined", sep="\t", header=TRUE, row.names=1)
numgenes <- dim(blc.counts)[1]
Data.Counts <- blc.counts[1:numgenes,1:numsamples]
Total.Counts <- as.numeric(apply(blc.counts[1:numgenes,1:numsamples], 2, sum))
Norm=scale(Data.Counts, scale=Total.Counts, center=FALSE)
Norm_Final <- t(Norm)
cor <- (cor(Norm_Final))^5

blc <- read.table(file="R_network_pipeline/BioLincs/blc.combined.summary", sep="\t", header=TRUE, row.names=1, quote="")
taxa <- c('Pro', 'Roseo', 'SAR11_', 'SAR116_', 'SAR324', 'SAR406', 'SAR86')

list <- c()
h <- 0
for(i in 1:length(taxa)) {
	for(j in 1:length(functions)) {
		name <- paste(taxa[i], functions[j], sep=".")
		h <- h + 1
		list[h] <- name
	}
}
list2 <- strsplit(list, ".", fixed=TRUE)

labels=c("Pro.Carbon_Fixation", "Pro.Cell_Cycle", "Pro.Central_Carbon_Metabolism", "Pro.Citrate_Cycle", "Pro.Cofactor_vitamin_metabolism", "Pro.DNA/Chromosome_Processing", "Pro.Glycolysis", "Pro.Oxidative_phosphorylation", "Pro.Photosynthesis", "Pro.Protein_RNA_Processing", "Pro.Signal_Transduction", "Pro.Transcription", "Pro.Translation", "Pro.Transporters", "SAR11_.Carbohydrate_Metabolism", "SAR11_.Central_Carbon_Metabolism", "SAR11_.Citrate_Cycle", "SAR11_.Cofactor_vitamin_metabolism", "SAR11_.DNA/Chromosome_Processing", "SAR11_.Glycolysis", "SAR11_.Oxidative_phosphorylation", "SAR11_.Protein_RNA_Processing", "SAR11_.Signal_Transduction", "SAR11_.Transcription", "SAR11_.Translation", "SAR11_.Transporters", "SAR116_.Central_Carbon_Metabolism", "SAR116_.Chemotaxis_and_Flagellar_Assembly", "SAR116_.Citrate_Cycle", "SAR116_.Cofactor_vitamin_metabolism", "SAR116_.DNA/Chromosome_Processing", "SAR116_.Glycolysis", "SAR116_.Oxidative_phosphorylation", "SAR116_.Protein_RNA_Processing", "SAR116_.Signal_Transduction", "SAR116_.Translation", "SAR116_.Transporters", "SAR86.Central_Carbon_Metabolism", "SAR86.Citrate_Cycle", "SAR86.Cofactor_vitamin_metabolism", "SAR86.DNA/Chromosome_Processing", "SAR86.Glycolysis", "SAR86.Oxidative_phosphorylation", "SAR86.Signal_Transduction", "SAR86.Transcription", "SAR86.Translation", "SAR86.Transporters", "Roseo.Central_Carbon_Metabolism", "Roseo.Chemotaxis_and_Flagellar_Assembly", "Roseo.Citrate_Cycle", "Roseo.Cofactor_vitamin_metabolism", "Roseo.Oxidative_phosphorylation", "Roseo.Protein_RNA_Processing", "Roseo.Transcription", "Roseo.Translation", "Roseo.Transporters", "SAR406.Citrate_Cycle", "SAR406.Cofactor_vitamin_metabolism", "SAR406.Oxidative_phosphorylation", "SAR406.Protein_RNA_Processing", "SAR406.Transcription", "SAR406.Translation", "SAR406.Transporters")

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
dat <- melt(matrix)
colnames(dat) <- c("cat1", "cat2", "val")

base_size=9
jpeg(file="blc.crosscompare3.jpg", quality=100, res=600, height=11, width=20, units="in")
ggplot(dat, aes(cat1, cat2)) + geom_tile(aes(fill = val), colour = "white") + scale_fill_gradient2(low = "red", mid='white', high = "blue") + theme(legend.position = "none", axis.ticks = element_blank(), axis.text.x = element_text(size = base_size * 0.9, angle = 30, hjust = 1, vjust = 1, face = 'bold', colour = "black"), axis.text.y = element_text(face = 'bold', colour = "black"))
dev.off()

base_size=9
jpeg(file="blc.crosscompare3.jpg", quality=100, res=600, height=11, width=20, units="in")
 heatmap(matrix, Rowv=NA, Colv=NA)
dev.off()


write.table(matrix, file="blc.pairwise.comapre", sep="\t")


###############################################



for(i in 1:dim(matrix)[1]) {
	for(j in 1:dim(matrix)[1]) {
		matrix3[i, j] <- matrix[i, j]
		matrix3[j, i] <- matrix2[i, j]
	}
}
library(reshape2)
library(ggplot2)
dat <- melt(matrix3)
colnames(dat) <- c("cat1", "cat2", "val")

base_size=9
jpeg(file="blc.crosscompare2.jpg", quality=100, res=600, height=11, width=20, units="in")
ggplot(dat, aes(cat1, cat2)) + geom_tile(aes(fill = val), colour = "white") + scale_fill_gradient2(low = "red", mid='white', high = "blue") + theme(legend.position = "none", axis.ticks = element_blank(), axis.text.x = element_text(size = base_size * 0.9, angle = 30, hjust = 1, vjust = 1, face = 'bold', colour = "black"), axis.text.y = element_text(face = 'bold', colour = "black")) + geom_abline(size=2, colour='grey80')
dev.off()



########################################################


jpeg(file="CANON.KEGG.connectivity3.jpg", quality=100, res=600, height=7, width=10, units="in")
for(i in 1:12) {
	clusters <- subset(canon, grepl(functions[i], canon$Final_Annote))
	clusters <- clusters[order(clusters$Connectivity, decreasing=T),]
	plot(clusters$Scaled.Connectivity[1:1000], col=mypal[i], ylim = c(0, 1), xlim=c(1, 300), ylab="Scaled Connectivity", xlab="Connectivity-Ordered Transcripts", pch=16)
	par(new=T)
}
labels <- c("Amino Acid Metabolism", "DNA and Chromosome Processing", "Translation", "Transporters", "Cofactor and Vitamin Metabolism", "Carbohydrate Metabolism",   "Oxidative Phosphorylation", "Fatty Acid Metabolism", "Carbon Fixation", "Transcription", "Photosynthesis")
legend("topright", legend=labels[1:11], col=mypal[1:11], pch=15, cex=1.3, bty="n")
dev.off()


# BioLincs
blc <- read.table(file="R_network_pipeline/BioLincs/blc.combined.summary", sep="\t", header=TRUE, row.names=1, quote="")

#functions <- c("Amino_Acid_metabolism", "Aminoacyl_tRNA_Biosynthesis", "Carbohydrate_Metabolism", "Carbon_Fixation", "Cell_Cycle", "Central_Carbon_Metabolism", "Chemotaxis_and_Flagellar_Assembly", "Citrate_Cycle", "Cofactor_vitamin_metabolism", "DNA/Chromosome_Processing", "Fatty_acid_metabolism", "Glycan_biosynthesis", "Glycolysis", "Hydrolase", "Isomerase", "Ligase", "Lipid_Metabolism", "Lyase", "Nitrogen_metabolism", "Nucleotide_metabolism", "Oxidative_phosphorylation", "Oxidoreductase", "Peptidase", "Photosynthesis", "Protein_RNA_Processing", "Proteorhodopsin", "Secondary_metabolism", "Secretion_System", "Signal_Transduction", "Sulfur_metabolism", "Transcription", "Transferase", "Translation", "Transport_Catabolism", "Transporters", "Unassigned", "Xenobiotic_degradation")
#functions <- c("Amino_Acid_metabolism", "Carbohydrate_Metabolism", "Carbon_Fixation", "Cell_Cycle", "Central_Carbon_Metabolism", "Chemotaxis_and_Flagellar_Assembly", "Citrate_Cycle", "Cofactor_vitamin_metabolism", "DNA/Chromosome_Processing", "Fatty_acid_metabolism", "Glycolysis", "Nitrogen_metabolism", "Nucleotide_metabolism", "Oxidative_phosphorylation", "Photosynthesis", "Protein_RNA_Processing", "Secretion_System", "Signal_Transduction", "Sulfur_metabolism", "Transcription", "Translation", "Transporters")

functions <- c("Amino_Acid_metabolism", "DNA/Chromosome_Processing", "Translation", "Transporters", "Cofactor_vitamin_metabolism", "Carbohydrate_Metabolism",   "Oxidative_phosphorylation", "Fatty_acid_metabolism", "Carbon_Fixation", "Transcription", "Photosynthesis")

jpeg(file="blc.KEGG.connectivity.jpg", quality=100, res=600, height=7, width=10, units="in")
for(i in 1:12) {
	clusters <- subset(blc, grepl(functions[i], blc$Final_Annote))
	clusters <- clusters[order(clusters$Connectivity, decreasing=T),]
	plot(clusters$Scaled.Connectivity[1:1000], col=mypal[i], ylim = c(0, 1), xlim=c(1, 300), ylab="Scaled Connectivity", xlab="Connectivity-Ordered Transcripts", pch=16)
	par(new=T)
}
labels <- c("Amino Acid Metabolism", "DNA and Chromosome Processing", "Translation", "Transporters", "Cofactor and Vitamin Metabolism", "Carbohydrate Metabolism",   "Oxidative Phosphorylation", "Fatty Acid Metabolism", "Carbon Fixation", "Transcription", "Photosynthesis")
legend("topright", legend=labels[1:11], col=mypal[1:11], pch=15, cex=1.3, bty="n")
dev.off()



for(i in 1:dim(matrix)[1]) {
	for(j in 1:dim(matrix)[1]) {
		matrix3[i, j] <- matrix[i, j]
		matrix3[j, i] <- matrix2[i, j]
	}
}














#k <- read.table(file="Network_attr_Files/Canon_combined_Analysis_outputs/Alltaxa_Network_att", sep="\t", header=TRUE, row.names=1)

############# How to make the combined annotation files- don't do again
mod <- data.frame(read.table(file="R_network_pipeline/CANON/canon.module", sep="\t", header=TRUE, row.names=1))
kegg <- data.frame(read.table(file="R_network_pipeline/Cluster.kegg.annotations.list", sep="\t", header=TRUE, row.names=1, nrows=544498, as.is=TRUE, quote=""))

kegg <- data.frame(read.table(file="R_network_pipeline/CANON/canon.KEGG", sep="\t", header=TRUE, row.names=1, nrows=544498, as.is=TRUE, quote=""))
counts <- data.frame(read.table(file="R_network_pipeline/CANON/Combined_Canon", sep="\t", header=TRUE, row.names=1))
k <- data.frame(read.table(file="R_network_pipeline/CANON/canon_Network_att", sep="\t", header=TRUE, row.names=1))

con <- intersect(row.names(k), row.names(kegg))

kegg2 <- kegg[con,]
combined <- cbind(counts, mod, k, kegg2)

write.table(combined, file="R_network_pipeline/CANON/canon.combined.summary2", sep="\t", quote=FALSE)

##############################################################
############### Plot abundance of microbial groups ##########
##############################################################
combined <- read.table(file="R_network_pipeline/CANON/canon.combined.summary2", sep="\t", quote="", header=TRUE, row.names=1)
#string <- strsplit(row.names(combined), "Cluster")
#taxa <- c(1:length(row.names(combined)))
#for(i in 1:length(taxa)) {taxa[i] <- string[[i]][1]}
#combined2 <- cbind(combined, taxa)
#combined2[,1:35][which(combined2$taxa == "Ostreo_")]

row.names(combined)[which(combined$Organism == "Ostreococcus")]

#################### Try plotting with ggplot2
############### Set color palette using RColorBrewer
library(RColorBrewer)
mypal <- c("green4", brewer.pal(9, "Set1")); mypal[7] = "gold"

list <- c("Ostreococcus", "SAR11_", "SAR116_", "SAR86", "Euryarchaota", "Roseobacter", "SAR92", "SAR406", "Flavobacteria", "ARCTIC")
result <- c()
for(i in 1:length(list)) {counts <- apply(combined[(combined$Organism == list[i]),25:59], 2, sum) ; result <- data.frame(rbind(result, counts))}
row.names(result) <- list
res <- 100*scale(result, scale=apply(result, 2, sum), center=F)
#dat1 <- melt(t(res)); colnames(dat1) <- c("time", "taxa", "abund")
#ggplot(data=dat1, aes(x=time, y=abund, group=taxa, colour=taxa)) + geom_line(size=1.5)+ scale_colour_manual(values=mypal)

################ Plot with regular R functions
res2 <- t(res)
jpeg(file="CANON.taxa.abundance.tiff", quality=100, res=600, height=10, width=7, units="in")
par(mfrow=c(5,2), mar=c(2, 3, 2, 1))
for(i in 1:10) {
	color <- mypal[i]
	axis <- c(1:35)
	title <- list[i]
	print(title)
	plot(axis, res2[,i], col=color, type="l", ylab=NA, lwd=2, main=title, cex.main=1.5, axes=F)
	rect(xleft = c(0, 4.5, 10.5, 15.5, 22.5, 30.5), xright=c(1.5, 7.5, 12.5, 19.5, 26.5, 33.5), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(100, 100, 100, 100, 100, 100), col='grey90', border=NA)
	par(new=T)
	plot(axis, res2[,i], col=color, type="l", ylab=NA, lwd=2, main=title, cex.main=1.5)
	#axis(1, at = c(0.1, 5, 10, 15, 20, 25, 30, 35), labels = c(0, 5, 10, 15, 20, 25, 30, 35), cex.axis=1.4)
	box(lty=1)
}


dev.off()



######################################################

# Make boxplot- ambiguously helpful
jpeg(file="CANON.KEGG.boxplot.jpg", quality=100, res=300, height=7, width=6, units="in")
ggplot(kegg, aes(x=Final_Annote, y=Connectivity)) + geom_boxplot() + coord_flip()
dev.off()

pathways <- c("Amino_acid_metabolism",  "Carbohydrate_metabolism", "Translation", "Nucleotide_metabolism", "Carbon_fixation", "Photosynthesis", "Cell_Cycle", "Replication_Repair", "Central_Carbon_Metabolism",  "Chromosome_associated", "Citrate_cycle", "Cofactor_vitamin_metabolism", "Flagellar_assembly", "Glycan_biosynthesis", "Glycolysis", "Hydrolase", "Isomerase", "Ligase", "Lipid_metabolism", "Lyase", "Methane_metabolism", "Nitrogen_metabolism",  "Oxidative_phosphorylation", "Oxidoreductase", "Pentose_phosphate_pathway", "Peptidase",  "Protein_RNA_Processing", "Secondary_metabolism", "Secretion_system", "Signal_transduction", "Sulfur_metabolism", "Transcription", "Transferase",  "Transport_Catabolism", "Transporters", "Chemotaxis", "Unassigned", "Bacteriorhodopsin", "Xenobiotic_degradation")

topten <- c("Nucleotide_metabolism", "Amino_acid_metabolism",  "Carbohydrate_metabolism", "Aminoacyl_tRNA_Biosynthesis", "Translation", "Photosynthesis", "Cofactor_vitamin_metabolism", "DNA/Chromosome_Processing", "Oxidative_phosphorylation", "Transcription")

############### Set color palette using RColorBrewer
library(RColorBrewer)
mypal <- c("green4", brewer.pal(9, "Set1"))
#subset <- sort(kegg$Connectivity[which(kegg$Final_Annote == "Transporters")], decreasing=T)
#plot(subset, ylim=c(0, 400), xlim=c(0, 200))
#par(new=T)
#set <- subset(kegg, grepl("Carbohydrate_metabolism", kegg$Final_Annote))

jpeg(file="CANON.KEGG.connectivity.jpg", quality=100, res=300, height=7, width=10, units="in")
for(i in 1:10) {
	
	subset <- sort(kegg$Connectivity[which(kegg$Final_Annote == topten[i])], decreasing=T)
	plot(subset, col=mypal[i], ylim = c(0, 400), xlim=c(1, 500), ylab="Connectivity", xlab="Connectivity-Ordered Transcripts")
	par(new=T)
}
#labels1 <- c("Amino_acid_metabolism",  "Carbohydrate_metabolism", "Translation", "Nucleotide_metabolism", "Carbon_fixation", "Photosynthesis", "Cell_Cycle", "Replication_Repair", "Central_Carbon_Metabolism",  "Chromosome_associated")
#labels2 <- c( "Citrate_cycle", "Cofactor_vitamin_metabolism", "Flagellar_assembly", "Glycan_biosynthesis", "Glycolysis", "Hydrolase", "Isomerase", "Ligase", "Lipid_metabolism", "Lyase")
#labels3 <- c("Methane_metabolism", "Nitrogen_metabolism",  "Oxidative_phosphorylation", "Oxidoreductase", "Pentose_phosphate_pathway", "Peptidase",  "Protein_RNA_Processing", "Secondary_metabolism", "Secretion_system", "Signal_transduction")
#labels4 <- c("Sulfur_metabolism", "Transcription", "Transferase",  "Transport_Catabolism", "Transporters", "Chemotaxis", "Unassigned", "Bacteriorhodopsin", "Xenobiotic_degradation")
legend("topright", legend=topten, col=mypal[1:10], pch=15, cex=1, bty="n")
dev.off()



Result <- c()
Final <- c()
for(i in 1:length(pathways)) {
	subset <- sort(kegg$Connectivity[which(kegg$Final_Annote == name[i])], decreasing=T)
	mean <- mean(subset)
	stdev <- sd(subset)
	Result <- c(mean(subset), sd(subset), length(subset))
	Final <- rbind(Final, Result)
}
cbind(Final, pathways)
colnames(Final) <- c("Mean", "STDEV", "Members", "Pathway")
dat1 <- melt(t(Final))


bymedian <- with(InsectSprays, reorder(spray, count, median))

keggplot <- with(kegg, reorder(kegg$Final_Annote, kegg$Connectivity, median))
ggplot(keggplot, aes(x=Final_Annote, y=Connectivity)) + geom_boxplot() + coord_flip()

ggplot(kegg, aes(x=Final_Annote, y=Connectivity)) + geom_boxplot() + coord_flip()



Result <- data.frame()
#kclust <- data.frame(cbind(colnames(Norm_Final), k, colclust))
sequ <- seq(0, 1, by=0.05)
#sequ <- seq(0.9, 1, by=0.005)
#sequ <- seq(0.9, 1, by=0.001)
	for(i in 2:length(sequ)) {
	quant_right <- quantile(kegg$Connectivity, sequ[i])
	quant_left <- quantile(kegg$Connectivity, sequ[i-1])
	names <- as.character(kegg$Final_Annote[which(kegg$Connectivity > quant_left & kegg$Connectivity < quant_right )])
	#string <- strsplit(names, "_")
	#print(quant)
	#list <- c("SAR11" = 0, "SAR116" = 0, "SAR86"=0, "eury"=0, "SAR92"=0, "Flavo"=0, "Roseo"=0, "Arctic"=0, "SAR406"=0, "Ostreo" = 0)
	annote <- c("Amino_acid_metabolism"=0, "Bacteriorhodopsin"=0, "Carbohydrate_metabolism"=0, "Carbon_fixation"=0, "Cell_Cycle"=0, "Cental_Carbon_Metabolism"=0, "Chemotaxis"=0, "Chromosome_associated"=0, "Citrate_cycle"=0, "Cofactor_vitamin_metabolism"=0, "Flagellar_assembly"=0, "Glycan_biosynthesis"=0, "Glycolysis"=0, "Hydrolase"=0, "Isomerase"=0, "Ligase"=0, "Lipid_metabolism"=0, "Lyase"=0, "Methane_metabolism"=0, "Nitrogen_metabolism"=0, "Nucleotide_metabolism"=0, "Oxidative_phosphorylation"=0, "Oxidoreductase"=0, "Pentose_phosphate_pathway"=0, "Peptidase"=0, "Photosynthesis"=0, "Protein_RNA_Processing"=0, "Replication_Repair"=0, "Secondary_metabolism"=0, "Secretion_system"=0, "Signal_transduction"=0, "Sulfur_metabolism"=0, "Transcription"=0, "Transferase"=0, "Translation"=0, "Transport_Catabolism"=0, "Transporters"=0, "Unassigned"=0, "Xenobiotic_degradation"=0)

		for(j in 1:length(names)) {
			#path <- kegg$Final_Annote[which(kegg$Cluster_ID == names[j])]
			annote[[names[j]]] = annote[[names[j]]] + 1
			#list[string[[j]][1]] = list[string[[j]][1]] + 1
		}
	final <- as.numeric(annote)
	Result <- data.frame(rbind(Result, final))
	#Result <- data.frame(SAR11 = list[1], SAR116 = list[2], SAR86=list[3], eury=list[4], SAR92=list[5], Flavo=list[6], Roseo=list[7], Arctic=list[8], SAR406=list[9], Ostreo=list[10])
	#Result[i,] <- t(data.frame(list))
	}
colnames(Result) <- name
row.names(Result) <- sequ[1:(length(sequ)-1)]
stdev <- apply(Result, 1, sd)
fin <- rbind(Result, stdev)
row.names(fin)[21] <- "stdev"
write.table(fin, sep="\t", file="output")
#######################################################################################


