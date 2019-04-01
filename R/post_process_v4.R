path_resultat=commandArgs(TRUE)[1]
setwd(path_resultat)


library(phyloseq)
library(vegan)



#Données brutes

otu_table=read.table("OTU_distribution_phyloseq.txt",h=T,sep="\t", row.names=1) # Modif 24/05/2016
otu_table[,-1]-> otu_table # remove seed sequences

do.call('rbind',strsplit(as.character(otu_table[,"LCA.taxonomy"]),";", fixed=TRUE))-> taxonomy 
rownames(taxonomy)=rownames(otu_table)
labels_taxonomy=c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
colnames(taxonomy) <- labels_taxonomy[1:ncol(taxonomy)]

otu_table[,which(colSums(otu_table[,1:(ncol(otu_table)-4)])>0)]->otu_table # Remove empty samples
otu_taxonomy= phyloseq(otu_table( otu_table, taxa_are_rows = TRUE), tax_table(taxonomy))

write.table(round(estimate_richness(otu_taxonomy),2), "../.Richness_diversity.tmp",sep="\t", quote=F, col.names=T, row.names=T)




jpeg("Heatmap_OTUs.jpg")
plot_heatmap(otu_taxonomy, low = "#66CCFF", high = "#000033", na.value = "white") 
dev.off()



jpeg("Rarefaction_curves_OTUs.jpg")
rarecurve(t(otu_table), step=100,  xlab = "Sequences", ylab = "OTUs")
dev.off()


jpeg("Richness_Diversity_OTUs.jpg")
plot_richness(otu_taxonomy,  measures=c("Chao1", "ACE", "Shannon"))
dev.off()

# 



#Données normalisées

otu_table_norm=read.table("../OTU_distribution_tax_normalized_LCA.txt",h=T,sep="\t", row.names=1)

otu_table_norm[,1:(ncol(otu_table_norm)-1)]-> otu_table_norm

if(is.data.frame(otu_table_norm)) # si 1 colonne -> vector
{

otu_table_norm[which(rowSums(otu_table_norm)>0),]-> otu_table_norm
otu_norm=phyloseq(otu_table(otu_table_norm, taxa_are_rows=T))

write.table(round(estimate_richness(otu_norm),2), "../.Richness_diversity_norm.tmp",sep="\t", quote=F, col.names=T, row.names=T)

jpeg("Rarefaction_curves_OTUs_normalized.jpg")
rarecurve(t(otu_table_norm), step=100,  xlab = "Sequences", ylab = "OTUs")
dev.off()


jpeg("Richness_Diversity_OTUs_normalized.jpg")
plot_richness(otu_norm,  measures=c("Chao1", "ACE", "Shannon"))
dev.off()


jpeg("Heatmap_OTUs_normalized.jpg")
plot_heatmap(otu_norm, low = "#66CCFF", high = "#000033", na.value = "white")
dev.off()


}

rm(path_resultat)
rm(labels_taxonomy)
rm(taxonomy)

save.image()


# derniers changements : le 28/07/2017
