## Modification du 21/7

library(picante)


clade_name = commandArgs(TRUE)[1]
path_phylogeny=commandArgs(TRUE)[2]	### phylo d'un profil
path_otuFile=commandArgs(TRUE)[3]	### otus d'un profil
profile_fig=commandArgs(TRUE)[4]

path_otuDistrib = commandArgs(TRUE)[5]   ### ajout 1/9/2014

phylogeny = read.tree(path_phylogeny)	
otuFile = read.table(path_otuFile)
otuDistrib = read.table(path_otuDistrib, header=T, sep="\t") 
otuDistrib[,-dim(otuDistrib)[2]]->otuDistrib #15/6/2017 dernière colonne ("NA")

#mntd_CladeTree=commandArgs(TRUE)[5]
#depth_deepest=commandArgs(TRUE)[6]

#file = commandArgs(TRUE)[5]    ### modif 1/9/2014
file = commandArgs(TRUE)[6]    ### modif 1/9/2014

nri_File = commandArgs(TRUE)[7]    ### ajout 2/9/2014

############

otuVector = as.vector(otuFile)  ## récupérer les noms des OTUs affiliés à un clade
t_otuVector = t(otuVector)	## transposée
otuNumber = length(t_otuVector)	## nombre d'otus

#otuMatrix=matrix(1,1,length(t_otuVector))	## générer une matrice [1, nombre d'OTUs], avec des uns comme valeurs (abondance = 1)  ### modif 1/9/2014
#colnames(otuMatrix)= t_otuVector	### renommer les colonnes de la matrice avec les noms des OTUs   ### modif 1/9/2014

################# ajout 1/9/2014
### matrice avec abondances réelles   

subset_otuDistrib = subset(otuDistrib, otuDistrib[,1] %in% t_otuVector)	

# 15/6/2017
# mat=data.matrix(subset_otuDistrib)
# dimnames(mat)<-NULL
# mat=mat[,-1]
# t_mat=t(mat)
# 
# rownames(t_mat)=as.vector(colnames(subset_otuDistrib)[-1])
# colnames(t_mat)=as.vector(subset_otuDistrib[,1])
# 
# otuMatrix=t_mat

mat=subset_otuDistrib[,-1] # Colonne des OTUs 
mat=mat[,1:(dim(mat)[2]-2)]
t_mat=t(mat)
colnames(t_mat)=as.vector(subset_otuDistrib[,1])

otuMatrix=t_mat
#15/07/2017


### récupérer le nombre d'échantillons dans lequel on retouve ce clade
 
VectorSample = as.vector(colSums(mat)) ## vecteur des sommes des colonnes (échantillons) de la matrice des OTUs de clade
l_VectSample = length (VectorSample)  ### taille du vecteur
for (i in 1: l_VectSample) {    ### si somme=0 (un échantillon ne contient aucun OTU du clade) la taille du vecteur est décrémenter de 1
  if (VectorSample[i] == 0) {
    l_VectSample = l_VectSample -1
  }
}

################# fin ajout 1/9/2014

#############

CladeTree=drop.tip(phylogeny, which(!phylogeny$tip.label %in% t_otuVector)) ## extraire le sous arbre des otus
patristic_CladeTree = cophenetic(CladeTree)

## Calcul de la somme de colonnes pour calculer un seul mntd 21/07/2015
abondance=matrix(0,1, dim(otuMatrix)[2]) 
abondance[1,]=colSums(otuMatrix) 
colnames(abondance)=colnames(otuMatrix)
rownames(abondance)="Clade"


mntd_CladeTree = round(mntd(abondance, patristic_CladeTree, abundance.weighted=T), digits=4) 
mntd_CladeTree_pa = round(mntd(abondance, patristic_CladeTree, abundance.weighted=F), digits=4) 

## 21/07/2015
  
depth_deepest = round(sum(CladeTree$edge.length)/length(CladeTree$tip.label), digits=4)

################################### ajout 1/9/2014

ses.mpd = ses.mpd(otuMatrix, patristic_CladeTree, null.model="taxa.labels", runs=500)
nri_values = as.vector(ses.mpd[,6])
nri_values_minus1 = -1 * nri_values
nri_values_minus = na.omit(nri_values_minus1)
mean_nri_otu = round(mean (nri_values_minus),digits=4)

sink(nri_File)
print("The net related index. Reveals non random patterns in community phylogenetic structure")
print (ses.mpd)
sink()

########################################  fin ajout 1/9/2014

sink(file)
print (l_VectSample)  ### modif 1/9/2014
print (mntd_CladeTree)   ### modif 1/9/2014
print(mntd_CladeTree_pa) ## 21/07/2015
print (depth_deepest)
sink()

#### fig

svg(profile_fig)
zoom(phylogeny, t_otuVector,  subtree=T, col="red", no.margin=F, root.edge=T, show.node.label = F, cex=0.4)
dev.off()

###########


