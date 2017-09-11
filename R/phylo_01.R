# 15/07/2017 correction multiples bugs
library(picante)
library(MASS)
library(cluster)  ## ajout 22/10/2014

profile_name = commandArgs(TRUE)[1]	### nom du profil


########## inputs

path_otuDistrib = commandArgs(TRUE)[2]	### otu_distrib_tmp pour tous les otus
path_phylogeny=commandArgs(TRUE)[3]	### phylo d'un profil
path_otuFile=commandArgs(TRUE)[4]	### otus d'un profil

otuDistrib = read.table(path_otuDistrib, header=T, sep="\t") #15/6/2017
otuDistrib[,-dim(otuDistrib)[2]]->otuDistrib #15/6/2017 dernière colonne ("NA")
phylogeny=read.tree(path_phylogeny)	
otuFile = read.table(path_otuFile)


########### outputs

profile_fig = commandArgs(TRUE)[5]	### nom du profil pour figure
otuPhylosor_File = commandArgs(TRUE)[6]	### nom du profil pour shared_branches
otuUnifrac_File = commandArgs(TRUE)[7]	### nom du profil pour uniq_branches
pd_File = commandArgs(TRUE)[8]
mpd_File = commandArgs(TRUE)[9]
mntd_File = commandArgs(TRUE)[10]
nri_File = commandArgs(TRUE)[11]
nti_File = commandArgs(TRUE)[12]

profile_output = commandArgs(TRUE)[13]

profile_unifrac_cluster_fig = commandArgs(TRUE)[14]   ## ajout 22/10/2014
profile_unifrac_pcoa_fig = commandArgs(TRUE)[15]   ## ajout 22/10/2014


resultat=matrix(,1,8)


##################################################

otuVector = as.vector(otuFile)  ## récupérer les noms des OTUs affiliés à un arbre phylo
t_otuVector = t(otuVector)	## transposée
otuNumber = length(t_otuVector)	## nombre d'otus

#otuMatrix=matrix(1,1,length(t_otuVector))	## générer une matrice [1, nombre d'OTUs], avec des uns comme valeurs (abondance = 1)
#colnames(otuMatrix)= t_otuVectort_mat	### renommer les colonnes de la matrice avec les noms des OTUs

################# matrice avec abondances réelles

subset_otuDistrib = subset(otuDistrib, otuDistrib[,1] %in% t_otuVector)	

#15/07/2017
# mat=data.matrix(subset_otuDistrib) # !!!
# dimnames(mat)<-NULL
mat=subset_otuDistrib[,-1] # Colonne des OTUs 
mat=mat[,1:(dim(mat)[2]-2)]
t_mat=t(mat)

# rownames(t_mat)=as.vector(colnames(subset_otuDistrib)[-1])
colnames(t_mat)=as.vector(subset_otuDistrib[,1])

otuMatrix=t_mat
#15/07/2017

#############

otuTree=drop.tip(phylogeny, which(!phylogeny$tip.label %in% t_otuVector)) ## extraire le sous arbre des otus

patristic_otuTree = cophenetic(otuTree)		## distances patristiques du sous arbre des otus
mean_patristicOtu = mean(patristic_otuTree)	## max des distances patristiques 
mean_Branch_otuTree=mean(otuTree$edge.length) ## max branch (edge) length

mean_pd_otuTree = round(mean(pd(otuMatrix, otuTree, include.root=T)[,1]), digits=4) ## phylogenetic distance
pd_otuTree = pd(otuMatrix, otuTree, include.root=T)[1]

sink(pd_File)
print("The sum of the total phylogenetic branch length for one or multiple samples")
print (pd_otuTree)
sink()

mean_mpd_otuTree = round(mean(mpd(otuMatrix, patristic_otuTree, abundance.weighted=T), na.rm = T), digits=4)	## mean pairwise distance
mpd_otuTree = mpd(otuMatrix, patristic_otuTree, abundance.weighted=T)

sink(mpd_File)
print("The mean pairwise distance separating taxa in a community")
print (mpd_otuTree)
sink()

mean_mntd_otuTree = round(mean(mntd(otuMatrix, patristic_otuTree, abundance.weighted=T), na.rm= T),digits=4)	## Mean nearest taxon distance
mntd_otuTree = mntd(otuMatrix, patristic_otuTree, abundance.weighted=T)

sink(mntd_File)
print("The mean nearest taxon distance for taxa in a community")
print (mntd_otuTree)
sink()

# nri clustering of phylogeny from root to leaves; reveals non random patternes in community phylogenetic structure
# The higher nri is (positive values), more phylogenetically clustered are the analyzed communities; more sensitive at tree scale 
# The lower nri is (negative values), more phylogenetically overdispersed are the analyzed communities. 

ses.mpd = ses.mpd(otuMatrix, patristic_otuTree, null.model="taxa.labels", runs=50)
nri_values = as.vector(ses.mpd[,6])
nri_values_minus = -1 * nri_values
#mean_nri_otu = mean(ses.mpd(otuMatrix, patristic_otuTree, null.model="taxa.labels", runs=50)[,6], na.rm = T)	
mean_nri_otu = round(mean (nri_values_minus, na.rm = T),digits=4)

sink(nri_File)
print("The net related index. Reveals non random patterns in community phylogenetic structure")
print (ses.mpd)
sink()

# nti : clustering of the terminal nodes; measures the branch tip clustering of a species of a community
# The higher nti is (positive values), more phylogenetically clustered are the analyzed communities; more sensitive at local scale
ses.mntd = ses.mntd(otuMatrix, patristic_otuTree, null.model="taxa.labels", runs=50)
nti_values = as.vector(ses.mntd[,6])
nti_values_minus = -1 * nti_values
#mean_nti_otu = mean(ses.mntd(otuMatrix, patristic_otuTree, null.model="taxa.labels", runs=50)[,6], na.rm = T)	
mean_nti_otu = round(mean(nti_values_minus, na.rm = T),digits=4)

sink(nti_File)
print("The nearest taxon index. Measures the branch tip clustering for taxa of a community")
print (ses.mntd)
sink()

############## fig

svg(profile_fig)
zoom(phylogeny, t_otuVector,  subtree=T, col="blue", no.margin=F, root.edge=T, show.node.label = F, cex=0.4)
dev.off()

############# community comparison 
### phylosor : Fraction of branch-length shared between two communities
### Unweighted UniFrac, a phylogenetic beta diversity metric of the unique (non-shared) fraction of total phylogenetic diversity (branch-length) between two communities

otuPhylosor = phylosor(otuMatrix, otuTree)
mean_otuPhylosor = mean(phylosor(otuMatrix, otuTree), na.rm = T)

sink(otuPhylosor_File)
print("The fraction of branch-length shared between two communities")
print (otuPhylosor)
sink()

otuUnifrac = unifrac(otuMatrix, otuTree)
mean_unifrac = round(mean(unifrac(otuMatrix, otuTree), na.rm = T),digits=4)

######## ajout 22/10/2014

unifrac.clusters = hclust(unifrac(otuMatrix[rowSums(otuMatrix)>0,], otuTree)) # DD : Modif 22/07/2015
svg(profile_unifrac_cluster_fig)   
plot (unifrac.clusters, xlab="Samples", sub="(Agglomeration method used: complete)")#16/06/2017
dev.off()

unifrac.pcoa = pcoa(unifrac(otuMatrix[rowSums(otuMatrix)>0,], otuTree), correction= "none", rn = NULL)  # DD : Modif 22/07/2015
labels = rownames(unifrac.pcoa$vectors)
xlabel= paste("PC1 = ", round(unifrac.pcoa$values[1,2]*100, digits=2),"%")
ylabel= paste("PC2 = ", round(unifrac.pcoa$values[2,2]*100, digits=2),"%")

xlim_lower = min(unifrac.pcoa$vectors[,1]) - 0.1;
xlim_upper = max(unifrac.pcoa$vectors[,1]) + 0.1;

ylim_lower = min(unifrac.pcoa$vectors[,2]) - 0.1;
ylim_upper = max(unifrac.pcoa$vectors[,2]) + 0.1;

svg(profile_unifrac_pcoa_fig) 
plot(unifrac.pcoa$vectors[,1], unifrac.pcoa$vectors[,2], xlab = xlabel, ylab = ylabel, pch=20, main="PCoA Unifrac metrics", ylim=c(ylim_lower, ylim_upper), xlim=c(xlim_lower, xlim_upper))
text(unifrac.pcoa$vectors[,1], unifrac.pcoa$vectors[,2]+.02,labels)
points(unifrac.pcoa$vectors[,1], unifrac.pcoa$vectors[,2])
dev.off()	

sink(otuUnifrac_File)
print("The unweighted UniFrac, the phylogenetic beta diversity metric of the unique (non-shared) fraction of total phylogenetic diversity (branch-length) between two communities")
print (otuUnifrac)
print (unifrac.pcoa)
sink()

######## fin ajout 22/10/2014

#################### printing results

resultat[1,1] = profile_name
resultat[1,2] = otuNumber
resultat[1,3] = mean_pd_otuTree
resultat[1,4] = mean_mpd_otuTree
resultat[1,5] = mean_mntd_otuTree
resultat[1,6] = mean_nri_otu
resultat[1,7] = mean_nti_otu
resultat[1,8] = mean_unifrac


write.matrix(resultat, profile_output, sep="\t")




