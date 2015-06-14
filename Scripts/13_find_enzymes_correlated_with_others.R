


####### CANON Glyoxylate cycle enzymes analysis
diel <- read.table(file="R_network_pipeline/CANON/All.canon.hra.out", sep="\t", header=TRUE, row.names=1)
diel.transcript <- (diel$Perm.FDR <= 0.1 & diel$Regression.FDR <= 0.1)
diel.2 <- cbind(diel, diel.transcript)

canon <- read.table(file="R_network_pipeline/CANON/canon.combined.summary2", sep="\t", header=TRUE, row.names=1)
final.canon <- merge(diel.2, canon, by="row.names", all=T)
row.names(final.canon) <- final.canon$Row.names

row.names(cor)[which(cor[,row.names(perox)[4]] > 0.9)]


############### Peroxidases
row.names(final.canon)[which(
perox <- subset(final.canon, grepl('perox', final.canon$Annotation))


counts <- canon[, 25:59]
cor <- cor(t(counts))

for(i in 1:length(row.names(perox))) {



taxa.spec <- canon[(canon$Organism == list[i]),25:59]
counts <- apply(taxa.spec, 2, sum)
mat <- scale(taxa.spec, scale=counts, center=F)
taxa.mat <- data.frame(cbind(taxa.mat, t(mat)))
}
cor <- (cor(taxa.mat))^5

