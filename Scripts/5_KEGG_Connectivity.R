
##########################################################################################################
############################### Identify KEGG Annotations of highly-connected nodes ######################
##########################################################################################################

#######################################################################
############## Load data and plot KEGG connectivity dotplot ###########
#######################################################################

library(RColorBrewer)
mypal <- c(brewer.pal(9, "Set1"), "darkblue", "violetred4", 'green4'); mypal[6] = "gold"

# CANON
canon <- read.table(file="R_network_pipeline/CANON/canon.final.combined.list2", sep="\t", header=TRUE, row.names=1, quote="")

#functions <- c("Amino_Acid_metabolism", "Aminoacyl_tRNA_Biosynthesis", "Carbohydrate_Metabolism", "Carbon_Fixation", "Cell_Cycle", "Central_Carbon_Metabolism", "Chemotaxis_and_Flagellar_Assembly", "Citrate_Cycle", "Cofactor_vitamin_metabolism", "DNA/Chromosome_Processing", "Fatty_acid_metabolism", "Glycan_biosynthesis", "Glycolysis", "Hydrolase", "Isomerase", "Ligase", "Lipid_Metabolism", "Lyase", "Nitrogen_metabolism", "Nucleotide_metabolism", "Oxidative_phosphorylation", "Oxidoreductase", "Peptidase", "Photosynthesis", "Protein_RNA_Processing", "Proteorhodopsin", "Secondary_metabolism", "Secretion_System", "Signal_Transduction", "Sulfur_metabolism", "Transcription", "Transferase", "Translation", "Transport_Catabolism", "Transporters", "Unassigned", "Xenobiotic_degradation")

#functions <- c("Amino_Acid_metabolism", "Carbohydrate_Metabolism", "Carbon_Fixation", "Cell_Cycle", "Central_Carbon_Metabolism", "Chemotaxis_and_Flagellar_Assembly", "Citrate_Cycle", "Cofactor_vitamin_metabolism", "DNA/Chromosome_Processing", "Fatty_acid_metabolism", "Glycolysis", "Nitrogen_metabolism", "Nucleotide_metabolism", "Oxidative_phosphorylation", "Photosynthesis", "Protein_RNA_Processing", "Secretion_System", "Signal_Transduction", "Sulfur_metabolism", "Transcription", "Translation", "Transporters")


############### First Plot
functions <- c("Amino_acid_metabolism", "DNA/Chromosome_Processing", "Translation", "Transporters", "Cofactor_vitamin_metabolism", "Carbohydrate_metabolism",   "Oxidative_phosphorylation", "Fatty_acid_metabolism", "Carbon_fixation", "Transcription", "Photosynthesis")
functions <- c("Amino_acid_metabolism", "DNA/Chromosome_Processing", "Translation", "Transporters", "Cofactor_vitamin_metabolism", "Carbohydrate_metabolism",   "Oxidative_phosphorylation", "Fatty_acid_metabolism", "Carbon_fixation", "Transcription", "Photosynthesis")


jpeg(file="CANON.KEGG.connectivity3.jpg", quality=100, res=600, height=7, width=10, units="in")
for(i in 1:11) {
	clusters <- subset(canon, grepl(functions[i], canon$Final_Annote))
	clusters <- clusters[order(clusters$Connectivity.TC, decreasing=T),]
	plot(clusters$Scaled.Connectivity.TC[1:1000], col=mypal[i], ylim = c(0, 1), xlim=c(1, 300), ylab="Scaled Connectivity", xlab="Connectivity-Ordered Transcripts", pch=16, cex.lab=1.2, font.lab=2, axes=F)
	par(new=T)
}
axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1), labels=c('0', '0.2', '0.4', '0.6', '0.8', '1'), cex.axis=1.3, las=2, font.axis='2', cex.axis=1.2)
axis(1, at=c(0, 50, 100, 150, 200, 250, 300), labels=c('0', '50', '100', '150', '200', '250', '300'), cex.axis=1.2, font.axis='2')
labels <- c("Amino Acid Metabolism", "DNA and Chromosome Processing", "Translation", "Transporters", "Cofactor and Vitamin Metabolism", "Carbohydrate Metabolism",   "Oxidative Phosphorylation", "Fatty Acid Metabolism", "Carbon Fixation", "Transcription", "Photosynthesis")
legend("topright", legend=labels[1:11], col=mypal[1:11], pch=15, cex=1.3, bty="n")
box(lwd=2)
dev.off()
################################


################### FOR DESEQ Data

############### First Plot
functions <- c("Amino_acid_metabolism", "DNA/Chromosome_Processing", "Translation", "Transporters", "Cofactor_vitamin_metabolism", "Carbohydrate_metabolism",   "Oxidative_phosphorylation", "Fatty_acid_metabolism", "Carbon_fixation", "Transcription", "Photosynthesis")

k <- read.table(file="R_network_pipeline/CANON/canon_Network_att.deseq", header=TRUE, row.names=1, sep="\t")

jpeg(file="CANON.KEGG.connectivity.jpg", quality=100, res=600, height=7, width=10, units="in")
for(i in 1:11) {
	clusters <- row.names(subset(canon, grepl(functions[i], canon$Final_Annote)))
	#clusters <- clusters[order(clusters$Connectivity, decreasing=T),]
	clusters <- sort(k[clusters,]$Scaled.Connectivity, decreasing=T)
		if(length(clusters) >= 1000) {n <- 1000}
		if(length(clusters) < 1000) {n <- length(clusters)}
	#plot(clusters[1:n], col=mypal[i], ylim = c(0, 1), xlim=c(1, 300), ylab="Scaled Connectivity", xlab="Connectivity-Ordered Transcripts", pch=16, cex.lab=1.2, font.lab=2, axes=F)
	plot(clusters[1:n], col=mypal[i], ylim = c(0, 1), xlim=c(1, 300), ylab="Scaled Connectivity", xlab="Connectivity-Ordered Transcripts", pch=16, cex.lab=1.2, font.lab=2, axes=F)
	par(new=T)
}
axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1), labels=c('0', '0.2', '0.4', '0.6', '0.8', '1'), cex.axis=1.3, las=2, font.axis='2', cex.axis=1.2)
axis(1, at=c(0, 50, 100, 150, 200, 250, 300), labels=c('0', '50', '100', '150', '200', '250', '300'), cex.axis=1.2, font.axis='2')
labels <- c("Amino Acid Metabolism", "DNA and Chromosome Processing", "Translation", "Transporters", "Cofactor and Vitamin Metabolism", "Carbohydrate Metabolism",   "Oxidative Phosphorylation", "Fatty Acid Metabolism", "Carbon Fixation", "Transcription", "Photosynthesis")
legend("topright", legend=labels[1:11], col=mypal[1:11], pch=15, cex=1.3, bty="n")
box(lwd=2)
dev.off()
########################### End DESeq data

###################

# BioLincs
blc <- read.table(file="R_network_pipeline/BioLincs/blc.final.combined.list2", sep="\t", header=TRUE, row.names=1, quote="")

#functions <- c("Amino_Acid_metabolism", "Aminoacyl_tRNA_Biosynthesis", "Carbohydrate_Metabolism", "Carbon_Fixation", "Cell_Cycle", "Central_Carbon_Metabolism", "Chemotaxis_and_Flagellar_Assembly", "Citrate_Cycle", "Cofactor_vitamin_metabolism", "DNA/Chromosome_Processing", "Fatty_acid_metabolism", "Glycan_biosynthesis", "Glycolysis", "Hydrolase", "Isomerase", "Ligase", "Lipid_Metabolism", "Lyase", "Nitrogen_metabolism", "Nucleotide_metabolism", "Oxidative_phosphorylation", "Oxidoreductase", "Peptidase", "Photosynthesis", "Protein_RNA_Processing", "Proteorhodopsin", "Secondary_metabolism", "Secretion_System", "Signal_Transduction", "Sulfur_metabolism", "Transcription", "Transferase", "Translation", "Transport_Catabolism", "Transporters", "Unassigned", "Xenobiotic_degradation")
#functions <- c("Amino_Acid_metabolism", "Carbohydrate_Metabolism", "Carbon_Fixation", "Cell_Cycle", "Central_Carbon_Metabolism", "Chemotaxis_and_Flagellar_Assembly", "Citrate_Cycle", "Cofactor_vitamin_metabolism", "DNA/Chromosome_Processing", "Fatty_acid_metabolism", "Glycolysis", "Nitrogen_metabolism", "Nucleotide_metabolism", "Oxidative_phosphorylation", "Photosynthesis", "Protein_RNA_Processing", "Secretion_System", "Signal_Transduction", "Sulfur_metabolism", "Transcription", "Translation", "Transporters")

functions <- c("Amino_Acid_metabolism", "DNA/Chromosome_Processing", "Translation", "Transporters", "Cofactor_vitamin_metabolism", "Carbohydrate_metabolism",   "Oxidative_phosphorylation", "Fatty_acid_metabolism", "Carbon_Fixation", "Transcription", "Photosynthesis")
functions <- c("Amino_acid_metabolism", "DNA/Chromosome_Processing", "Translation", "Transporters", "Cofactor_vitamin_metabolism", "Carbohydrate_metabolism",   "Oxidative_phosphorylation", "Fatty_acid_metabolism", "Carbon_fixation", "Transcription", "Photosynthesis")

jpeg(file="blc.KEGG.connectivity2.jpg", quality=100, res=600, height=7, width=10, units="in")
for(i in 1:11) {
	clusters <- subset(blc, grepl(functions[i], blc$Final_Annote))
	clusters <- clusters[order(clusters$Connectivity.TC, decreasing=T),]
	plot(clusters$Scaled.Connectivity.TC[1:1000], col=mypal[i], ylim = c(0, 1), xlim=c(1, 300), ylab="Scaled Connectivity", xlab="Connectivity-Ordered Transcripts", pch=16, axes=F, cex.lab=1.2, font.lab=2)
	par(new=T)
}
axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1), labels=c('0', '0.2', '0.4', '0.6', '0.8', '1'), cex.axis=1.3, las=2, font.axis='2', cex.axis=1.2)
axis(1, at=c(0, 50, 100, 150, 200, 250, 300), labels=c('0', '50', '100', '150', '200', '250', '300'), cex.axis=1.2, font.axis='2')
labels <- c("Amino Acid Metabolism", "DNA and Chromosome Processing", "Translation", "Transporters", "Cofactor and Vitamin Metabolism", "Carbohydrate Metabolism",   "Oxidative Phosphorylation", "Fatty Acid Metabolism", "Carbon Fixation", "Transcription", "Photosynthesis")
legend("topright", legend=labels[1:11], col=mypal[1:11], pch=15, cex=1.3, bty="n")
box(lwd=2)
dev.off()


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


