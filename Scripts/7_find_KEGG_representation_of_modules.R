
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

canon <- read.table(file="R_network_pipeline/CANON/canon.final.combined.list2", sep="\t", header=TRUE, row.names=1, quote="")

num <- 34
result <- c()
pathways <- c("Amino_acid_metabolism",  "Carbohydrate_metabolism", 'Fatty_acid_metabolism', "DNA/Chromosome_Processing", "Aminoacyl_tRNA_Biosynthesis", "Translation", "Nucleotide_metabolism", "Carbon_fixation", "Photosynthesis", "Cell_Cycle", "Central_Carbon_Metabolism", "Citrate_cycle", "Cofactor_vitamin_metabolism", "Chemotaxis_and_Flagellar_Assembly", "Glycan_biosynthesis", "Glycolysis", "Hydrolase", "Isomerase", "Ligase", "Lipid_metabolism", "Lyase", "Nitrogen_metabolism",  "Oxidative_phosphorylation", "Oxidoreductase", "Proteorhodopsin", "Peptidase",  "Protein_RNA_Processing", "Secondary_metabolism", "Secretion_system", "Signal_transduction", "Sulfur_metabolism", "Transcription", "Transferase", "Transport_Catabolism", "Transporters", "Unassigned", "Xenobiotic_degradation")

for(i in 0:34) {
	subset <- as.character(canon$Final_Annote[which(canon$Module == (i))])
	list <- sort(c("Amino_acid_metabolism"=0, "Aminoacyl_tRNA_Biosynthesis"=0, "Carbohydrate_metabolism"=0, "Carbon_fixation"=0, "Cell_Cycle"=0, "Central_Carbon_Metabolism"=0, "Chemotaxis_and_Flagellar_Assembly"=0, "Citrate_cycle"=0, "Cofactor_vitamin_metabolism"=0, "DNA/Chromosome_Processing"=0, "Fatty_acid_metabolism"=0, "Glycan_biosynthesis"=0, "Glycolysis"=0, "Hydrolase"=0, "Isomerase"=0, "Ligase"=0, "Lipid_metabolism"=0, "Lyase"=0, "Nitrogen_metabolism"=0, "Nucleotide_metabolism"=0, "Oxidative_phosphorylation"=0, "Oxidoreductase"=0, "Peptidase"=0, "Photosynthesis"=0, "Protein_RNA_Processing"=0, "Proteorhodopsin"=0, "Secondary_metabolism"=0, "Secretion_system"=0, "Signal_transduction"=0, "Sulfur_metabolism"=0, "Transcription"=0, "Transferase"=0, "Translation"=0, "Transport_Catabolism"=0, "Transporters"=0, "Unassigned"=0, "Xenobiotic_degradation"=0))
		for(j in 1:length(subset)) {
			list[subset[[j]][1]] <- list[subset[[j]][1]] + 1
		}
	final <- as.numeric(list)
	result <- data.frame(rbind(result, final))
}
colnames(result) <- sort(pathways)
row.names(result) <- c("Unassigned", 1:34)
#write.table(result, file='output', sep='\t', quote=F)
#x <- as.matrix(read.table(file='canon.modules', header=TRUE, row.names=1, sep='\t'))

pathsub <- c("Amino_acid_metabolism",  "Carbohydrate_metabolism", "Fatty_acid_metabolism", "DNA/Chromosome_Processing", "Aminoacyl_tRNA_Biosynthesis", "Translation", "Nucleotide_metabolism", "Carbon_fixation", "Photosynthesis", "Cell_Cycle", "Central_Carbon_Metabolism", "Citrate_cycle", "Cofactor_vitamin_metabolism", "Chemotaxis_and_Flagellar_Assembly", "Glycan_biosynthesis", "Glycolysis", "Lipid_metabolism", "Nitrogen_metabolism",  "Oxidative_phosphorylation", "Protein_RNA_Processing", "Secondary_metabolism", "Secretion_system", "Signal_transduction", "Sulfur_metabolism", "Transcription", "Transporters")
x <- t(result[,pathsub])
y <- t(x[order(row.names(x)),])
#z <- t(y)

names <- sort(c("Amino Acid Met.",  "Carbohydrate Met.", "Fatty Acid Met.", "DNA/Chromosome Processing", "Aminoacyl-tRNA Biosyn.", "Translation", "Nucleotide Met.", "Carbon Fixation", "Photosynthesis", "Cell Cycle", "Central Carbon Met.", "Citrate Cycle", "Cofactor/Vitamin Met.", "Chemotaxis/Flagella", "Glycan Biosyn.", "Glycolysis", "Lipid Met.","Nitrogen Met.",  "Oxidative Phosph.", "Protein/RNA Processing", "Secondary Met.", "Secretion System", "Signal Transduction", "Sulfur Met.", "Transcription", "Transporters"))

library(ggplot2)
library(reshape2)
library(grid)
dat1 <- melt(y)
colnames(dat1) <- c("Module", "Pathway", "Transcripts")

dat2 <- dat1
zeros <- row.names(dat1)[which(dat1$Transcripts == 0)]
dat2[zeros, 'Transcripts'] <- NA

base_size=9
#tiff(file="canon_heatmap.tiff", quality=100, res=600, height=19, width=16, units="in")
png(file="canon_heatmap.png", height=19, width=16, units='in', res=400)
ggplot(dat1, aes(Pathway, Module)) + geom_tile(aes(fill = Transcripts), colour = "white") + scale_fill_gradient2(low='white', high = "dodgerblue") + theme(axis.ticks = element_blank(), axis.text.x = element_text(size = 16.5, hjust = 1, angle=30, vjust = 1, face = 'bold', colour = "black"), axis.title.x = element_blank(), axis.title.y = element_blank(), legend.title = element_text(size=16), legend.key.size= unit(0.3, "in"), legend.text=element_text(size=14), axis.text.y = element_text(face = 'bold', colour = "black", size=17)) + geom_text(aes(fill = dat2$Transcripts, label = dat2$Transcripts), colour='grey30') + scale_x_discrete(labels=names)
dev.off()

base_size=9
jpeg(file="canon.heatmap.functions.jpg", quality=100, res=600, height=19, width=14, units="in")
ggplot(dat1, aes(Pathway, Module)) + geom_tile(aes(fill = Transcripts), colour = "white") + scale_fill_gradient2(low='white', high = "dodgerblue") + theme(axis.ticks = element_blank(), axis.text.x = element_text(size = 16.5, hjust = 1, angle=30, vjust = 1, face = 'bold', colour = "black"), axis.title.x = element_blank(), axis.title.y = element_blank(), legend.title = element_text(size=16), legend.key.size= unit(0.3, "in"), legend.text=element_text(size=14), axis.text.y = element_text(face = 'bold', colour = "black", size=17)) + geom_text(aes(fill = dat2$Transcripts, label = dat2$Transcripts), colour='grey30') + scale_x_discrete(labels=names)
dev.off()
############
#################################################
################# End
#################################################

########## Load big summary file
# BioLincs

blc <- read.table(file="R_network_pipeline/BioLincs/blc.final.combined.list2", sep="\t", header=TRUE, row.names=1, quote='')

num <- 18
result <- c()
pathways <- c("Amino_acid_metabolism",  "Carbohydrate_metabolism", 'Fatty_acid_metabolism', "DNA/Chromosome_Processing", "Aminoacyl_tRNA_Biosynthesis", "Translation", "Nucleotide_metabolism", "Carbon_fixation", "Photosynthesis", "Cell_Cycle", "Central_Carbon_Metabolism", "Citrate_cycle", "Cofactor_vitamin_metabolism", "Chemotaxis_and_Flagellar_Assembly", "Glycan_biosynthesis", "Glycolysis", "Hydrolase", "Isomerase", "Ligase", "Lipid_metabolism", "Lyase", "Nitrogen_metabolism",  "Oxidative_phosphorylation", "Oxidoreductase", "Proteorhodopsin", "Peptidase",  "Protein_RNA_Processing", "Secondary_metabolism", "Secretion_system", "Signal_transduction", "Sulfur_metabolism", "Transcription", "Transferase", "Transporters", "Unassigned", "Xenobiotic_degradation")

for(i in 0:18) {
	subset <- as.character(blc$Final_Annote[which(blc$Module == (i))])
	list <- sort(c("Amino_acid_metabolism"=0, "Aminoacyl_tRNA_Biosynthesis"=0, "Carbohydrate_metabolism"=0, "Carbon_fixation"=0, "Cell_Cycle"=0, "Central_Carbon_Metabolism"=0, "Chemotaxis_and_Flagellar_Assembly"=0, "Citrate_cycle"=0, "Cofactor_vitamin_metabolism"=0, "DNA/Chromosome_Processing"=0, "Fatty_acid_metabolism"=0, "Glycan_biosynthesis"=0, "Glycolysis"=0, "Hydrolase"=0, "Isomerase"=0, "Ligase"=0, "Lipid_metabolism"=0, "Lyase"=0, "Nitrogen_metabolism"=0, "Nucleotide_metabolism"=0, "Oxidative_phosphorylation"=0, "Oxidoreductase"=0, "Peptidase"=0, "Photosynthesis"=0, "Protein_RNA_Processing"=0, "Proteorhodopsin"=0, "Secondary_metabolism"=0, "Secretion_system"=0, "Signal_transduction"=0, "Sulfur_metabolism"=0, "Transcription"=0, "Transferase"=0, "Translation"=0, "Transporters"=0, "Unassigned"=0, "Xenobiotic_degradation"=0))
		for(j in 1:length(subset)) {
			list[subset[[j]][1]] <- list[subset[[j]][1]] + 1
		}
	final <- as.numeric(list)
	result <- data.frame(rbind(result, final))
}
colnames(result) <- sort(pathways)
row.names(result) <- c("Unassigned", 1:18)

#write.table(result, file='output', sep='\t', quote=F)

pathsub <- c("Amino_acid_metabolism",  "Carbohydrate_metabolism", "Fatty_acid_metabolism", "DNA/Chromosome_Processing", "Aminoacyl_tRNA_Biosynthesis", "Translation", "Nucleotide_metabolism", "Carbon_fixation", "Photosynthesis", "Cell_Cycle", "Central_Carbon_Metabolism", "Citrate_cycle", "Cofactor_vitamin_metabolism", "Chemotaxis_and_Flagellar_Assembly", "Glycan_biosynthesis", "Glycolysis", "Lipid_metabolism", "Nitrogen_metabolism",  "Oxidative_phosphorylation", "Protein_RNA_Processing", "Secondary_metabolism", "Secretion_system", "Signal_transduction", "Sulfur_metabolism", "Transcription", "Transporters")
x <- t(result[,pathsub])
y <- t(x[order(row.names(x)),])

names <- sort(c("Amino Acid Met.",  "Carbohydrate Met.", "Fatty Acid Met.", "DNA/Chromosome Processing", "Aminoacyl-tRNA Biosyn.", "Translation", "Nucleotide Met.", "Carbon Fixation", "Photosynthesis", "Cell Cycle", "Central Carbon Met.", "Citrate Cycle", "Cofactor/Vitamin Met.", "Chemotaxis/Flagella", "Glycan Biosyn.", "Glycolysis", "Lipid Met.","Nitrogen Met.",  "Oxidative Phosph.", "Protein/RNA Processing", "Secondary Met.", "Secretion System", "Signal Transduction", "Sulfur Met.", "Transcription", "Transporters"))

#x <- as.matrix(read.table(file='blc.modules', header=TRUE, row.names=1, sep='\t'))
library(ggplot2)
library(reshape2)
library(grid)
dat1 <- melt(y)
colnames(dat1) <- c("Module", "Pathway", "Transcripts")

dat2 <- dat1
zeros <- row.names(dat1)[which(dat1$Transcripts == 0)]
dat2[zeros, 'Transcripts'] <- NA

base_size=9
jpeg(file="blc.heatmap.functions.jpg", quality=100, res=600, height=15, width=14, units="in")
ggplot(dat1, aes(Pathway, Module)) + geom_tile(aes(fill = Transcripts), colour = "white") + scale_fill_gradient2(low='white', high = "dodgerblue") + theme(axis.ticks = element_blank(), axis.text.x = element_text(size = 16.5, hjust = 1, angle=30, vjust = 1, face = 'bold', colour = "black"), axis.title.x = element_blank(), axis.title.y = element_blank(), legend.title = element_text(size=16), legend.key.size= unit(0.3, "in"), legend.text=element_text(size=14), axis.text.y = element_text(face = 'bold', colour = "black", size=17)) + geom_text(aes(fill = dat2$Transcripts, label = dat2$Transcripts), colour='grey30') + scale_x_discrete(labels=names)
dev.off()
#################################################
################# End
#################################################






sum <- apply(result, 1, sum)
result[,39] <- sum
order <- result[order(result[,39], decreasing=TRUE),]

mypal = labels2colors(1:37)
library(ggplot2)
library(reshape2)
dat1 <- melt(t(result[1:37]))
colnames(dat1) <- c("KEGG", "Module", "Transcripts")
#jpeg(file="CANON_Modules.jpg", quality=100, res=500, height=8, width=8, units="in")
ggplot(dat1, aes(Module, Transcripts, fill=KEGG), guide=guide_legend(title="size")) + geom_bar(stat="identity") + coord_flip() + scale_fill_hue(l=50, c=60) + theme(axis.text.y  = element_text(colour="black", size=10, face="bold"), axis.text.x = element_text(colour="black", face="bold"), legend.title = element_text(size=12, face="bold"), legend.text = element_text(size = 10)) + scale_fill_manual(values=mypal)
dev.off()

#############################################################################################################################
############ Load raw count data and normalize by total counts in a given timepoint [column] [this is identical to script 1].
numsamples <- 35
Alltaxa.counts <- read.table(file="R_network_pipeline/CANON/Combined_Canon", sep="\t", header=TRUE, row.names=1)
numgenes <- dim(Alltaxa.counts)[1]
Data.Counts <- Alltaxa.counts[1:numgenes,1:numsamples]
Total.Counts <- as.numeric(apply(Alltaxa.counts[1:numgenes,1:numsamples], 2, sum))
Norm=scale(Data.Counts, scale=Total.Counts, center=FALSE)
Norm_Final <- t(Norm)
#time <- as.numeric(Alltaxa.counts[1,])

########### Load module .module file generated in script 1. 
Module_Colors <- read.table(file="R_network_pipeline/CANON/canon.modules.TC.10_2014", sep="\t", header=TRUE, row.names=1)
mod <- data.frame(sort(table(Module_Colors[,2]), decreasing=T))

#modules <- row.names(data.frame(sort(table(mergedColors), decreasing=T)[1:38]))
# Define modules (by color) for which eigengenes are to be calculated (color codes are identical to those used in script 1). 
modules <- c("steelblue1", "purple", "orange", "springgreen4", "red", "magenta", "steelblue4", "brown", "cyan", "turquoise3", "violetred4", "tomato1", "violet", "tan4", "turquoise4", "wheat3", "slategray2", "yellow", "tan1", "sienna2", "darkturquoise", "saddlebrown", "orangered1", "firebrick1", "orchid4", "royalblue3", "tomato4", "olivedrab1", "olivedrab4", "springgreen3", "skyblue", "chocolate3", "darkcyan", "aquamarine", "coral", "darkgoldenrod", "blueviolet")
eig <- moduleEigengenes(Norm_Final, as.character(Module_Colors[,2]), excludeGrey=TRUE)
data <- eig$eigengenes
data1 <- eig$averageExpr
#times = as.numeric(x[1,])
timelist = list(time)
#data1 <- cbind(data, times)
colnames(data) <- gsub("ME", "", colnames(data))
#data <- data.frame(scale(data))
### Set up output for eigengene data; the user will need to adjust the plotting parameters
jpeg(file="eigengenes.new.jpg", quality=100, res=300, height=10, width=7, units="in")
par(mfrow=c(6,3), mar=c(2, 1, 2, 1))
for(i in 1:15) {
#for(i in 1:(dim(data)[2]-1)) {
	color <- modules[i]
	#print(color)
	#ag <- aggregate(data[,color], timelist, mean)
	#axis = c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24)
	#axis = c(1, 3, 5, 6, 7, 9, 11, 12, 13, 14, 15, 17, 19, 21, 22, 23)
	axis <- c(1:35)
	#row.names(ag) <- axis
	title <- paste(c("Module ", i, "; ", color, "; ", "n=", mod[modules[i],]), collapse="")
	print(title)
	plot(axis, data[,modules[i]], col=color, type="l", ylab=NA, lwd=2, main=title, cex.main=1.5, axes=F, ylim=c(-0.26, 0.26))
	rect(xleft = c(0, 4.5, 10.5, 15.5, 22.5, 30.5), xright=c(1.5, 7.5, 12.5, 19.5, 26.5, 33.5), ybottom=c(-1, -1, -0.5, -0.5, -0.5, -0.5), ytop=c(1, 1, 1, 1, 1, 1), col='grey90', border=NA)
	par(new=T)
	plot(axis, data[,modules[i]], col=color, type="l", ylab=NA, lwd=2, main=title, cex.main=1.5, axes=F, ylim=c(-0.26, 0.26))
	axis(1, at = c(0.1, 5, 10, 15, 20, 25, 30, 35), labels = c(0, 5, 10, 15, 20, 25, 30, 35), cex.axis=1.4)
	lines(seq(0, 0, length.out=35), pch=22, lty=5)
	box(lty=1)
	#par(new=T)
	#plot(data[,i], col=color, type="l", ylab=NA, lwd=1)
	#par(new=T)
	#plot(times, data[,i], col=color, ylab=NA, ylim=c(-0.2, 0.9))
	#par(new=T)
}
dev.off()


############################################################################################################
########################## Calculate Taxa Representation Across the Modules Identified #####################
############################################################################################################

############### Set color palette using RColorBrewer
library(RColorBrewer)
mypal <- c("green4", brewer.pal(9, "Set1"))

###################################################
num <- 37
result <- c()
names <- c("Ostreococcus", "SAR11", "SAR116", "SAR86", "GII Euryarchaeota", "SAR92", "Flavobacteria", "Roseobacter", "ARCTIC96-BD19", "SAR406")
for(i in 0:num) {

	if(i == 0) {
	set <- subset(Module_Colors, grepl("grey", Module_Colors[,2]))
	subset <- as.character(row.names(set))
	}
	else {
	subset <- row.names(Module_Colors)[which(Module_Colors[,1] == i)]
	}
	string <- strsplit(subset, "Cluster")
	list <- c("Ostreo_"=0, "SAR11_" = 0, "SAR116_" = 0, "SAR86_"=0, "eury_"=0, "SAR92_"=0, "Flavo_"=0, "Roseo_"=0, "Arctic_"=0, "SAR406_"=0)
		for(j in 1:length(string)) {
			list[string[[j]][1]] = list[string[[j]][1]] + 1
		}
	final <- as.numeric(list)
	result <- data.frame(rbind(result, final))
}
colnames(result) <- names
row.names(result) <- c("Unassigned", 1:37)
sum <- apply(result, 1, sum)
result[,11] <- sum
order <- result[2:38,][order(result[,11], decreasing=TRUE),]

library(ggplot2)
library(reshape2)
dat1 <- melt(t(result[1:10]))
colnames(dat1) <- c("Taxon", "Module", "Transcripts")
dat1 <- dat1[order(dat1$Taxon, decreasing=FALSE),]
dat1$Module[1:350] <- 3

jpeg(file="CANON_Modules.jpg", quality=100, res=500, height=2, width=10, units="in")
ggplot(dat1, aes(Module, Transcripts, fill=Taxon), guide=guide_legend(title="size")) + geom_bar(stat="identity") + coord_flip() + scale_fill_hue(l=50, c=60) + theme(axis.text.y  = element_text(colour="black", size=10, face="bold"), axis.text.x = element_text(colour="black", face="bold"), legend.title = element_text(size=12, face="bold"), legend.text = element_text(size = 10)) + scale_fill_manual(values=mypal)
dev.off()




################# New plot for revised manuscript with all modules combined

result <- c()
names <- c("Ostreococcus", "SAR11", "SAR116", "SAR86", "GII Euryarchaeota", "SAR92", "Flavobacteria", "Roseobacter", "ARCTIC96-BD19", "SAR406")
subset <- row.names(Module_Colors)
string <- strsplit(subset, "Cluster")
list <- c("Ostreo_"=0, "SAR11_" = 0, "SAR116_" = 0, "SAR86_"=0, "eury_"=0, "SAR92_"=0, "Flavo_"=0, "Roseo_"=0, "Arctic_"=0, "SAR406_"=0)


mypal <- c("darkolivegreen2", "magenta", "cyan", "green4", "purple", "red", "yellow",  "steelblue4", "blue", "orange")

for(j in 1:length(string)) {
	list[string[[j]][1]] = list[string[[j]][1]] + 1
}

result <- data.frame(list)
row.names(result) <- names
Module <- seq(from=2, to=2, length.out=10)
dat1 <- cbind(result, Module, names)
colnames(dat1) <- c('Count', 'Module', 'Names')
dat1 <- dat1[order(dat1$Count, decreasing=TRUE),]
dat1$Module <- c(1:10)

jpeg(file="CANON_Modules.jpg", quality=100, res=500, height=8, width=8, units="in")
ggplot(dat1, aes(Module, Count, fill=Names)) + geom_bar(stat="identity") + coord_flip() + theme(axis.text.y  = element_text(colour="black", size=10, face="bold"), axis.text.x = element_text(colour="black", face="bold", size=12), legend.title = element_text(size=12, face="bold"), legend.text = element_text(size = 10), legend.position='bottom', panel.background = element_blank()) + scale_fill_hue(l=50, c=60)  + scale_fill_manual(values=mypal)
dev.off()


 
Module_Colors <- read.table(file="R_network_pipeline/BioLincs/blc.module.10_2014", sep="\t", header=TRUE, row.names=1)

result <- c()
names <- c("Prochlorococcus", "SAR11", "SAR116", "SAR86", "Roseobacter", "SAR406", "SAR324")
subset <- row.names(Module_Colors)
string <- strsplit(subset, "Cluster")
list <- c("Pro_"=0, "SAR11_" = 0, "SAR116_" = 0, "SAR86_"=0, "Roseo_"=0, "SAR406_"=0, "SAR324_" = 0)

mypal <- c("green4", "purple", "red", "yellow", "cyan", "steelblue4",  "blue", "magenta",  "orange")

for(j in 1:length(string)) {
	list[string[[j]][1]] = list[string[[j]][1]] + 1
}

result <- data.frame(list)
row.names(result) <- names
Module <- seq(from=2, to=2, length.out=7)
dat1 <- cbind(result, Module, names)
colnames(dat1) <- c('Count', 'Module', 'Names')
dat1 <- dat1[order(dat1$Count, decreasing=TRUE),]
dat1$Module <- c(1:7)

jpeg(file="CANON_Modules.jpg", quality=100, res=500, height=7, width=8, units="in")
ggplot(dat1, aes(Module, Count, fill=Names)) + geom_bar(stat="identity") + coord_flip() + theme(axis.text.y  = element_text(colour="black", size=10, face="bold"), axis.text.x = element_text(colour="black", face="bold", size=12), legend.title = element_text(size=12, face="bold"), legend.text = element_text(size = 10), legend.position='bottom', panel.background = element_blank()) + scale_fill_hue(l=50, c=60)  + scale_fill_manual(values=mypal)
dev.off()

