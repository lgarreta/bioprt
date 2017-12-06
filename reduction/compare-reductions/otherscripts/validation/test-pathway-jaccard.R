#!/usr/bin/Rscript

library ("clusteval")

# Jaccard function
jaccard <- function (a,b) {
	sizeInter = length (intersect (a,b))
	sizeUnion = length (union (a,b))

	distance =  1 - (sizeInter / sizeUnion)

	return (distance)
}

#----------------------------------------------------------
# LOAD DATA
#----------------------------------------------------------
# Load previous results of  PCA scores and 
# clustering of PDB structures
scores <- read.table ("pathways-lnx-win.scores")
structures <- read.table ("structures-groups-win.tbl")

#----------------------------------------------------------
# STRUCTURAL CLUSTERING
#----------------------------------------------------------
# Get the individual results by group
c1 = subset (structures, structures==1)
c2 = subset (structures, structures==2)
c3 = subset (structures, structures==3)
c4 = subset (structures, structures==4)

# Get the names of the proteins by group
listStructs = list (c1,c2,c3,c4)
cNames = mapply (rownames, listStructs)

#----------------------------------------------------------
# KMEANS CLUSTERING
#----------------------------------------------------------
# Run kmeans clustering of PCAs
kmeanscluster <- kmeans(scores, 4)
Kclusterprotein1 <- (kmeanscluster$cluster)

# Get kmeans results by resulting groups
kc = as.data.frame (Kclusterprotein1)
k1 = subset (kc, kc$Kclusterprotein1==1)
k2 = subset (kc, kc$Kclusterprotein1==2)
k3 = subset (kc, kc$Kclusterprotein1==3)
k4 = subset (kc, kc$Kclusterprotein1==4)

# Get the names of the proteins by groupd
listKmeans = list (k1,k2,k3,k4)
kNames = mapply (rownames, listKmeans)


# Compute the Jaccard distance between
# groups obtained by clustering structures
# and groups resulting by PCA analysis
kMatrix = matrix (nrow=4,ncol=4)
for (j in 1:4) {
	structCluster = cNames [[j]]
	for (i in 1:4) {
		kmeansCluster = kNames [[i]]
		kMatrix [j,i] = jaccard (structCluster,kmeansCluster)
	}
}

cat ("KMEANs vs Structural:\n")
cat ("Jaccard distances between structural and kmeans clustering of PCAs\n")
rownames (kMatrix) <- c("c1","c2","c3","c4")
colnames (kMatrix) <- c("c1'","c2'","c3'","c4'")
print (round (kMatrix,2))

#----------------------------------------------------------
# HIERARCHICAL CLUSTERING
#----------------------------------------------------------
HierCluster <- dist(scores, method = "euclidean")
Hclusterprotein1 <- hclust(HierCluster, method="ward.D")
groupsHclusterprotein1 <- cutree(Hclusterprotein1, k=4)

# Get kmeans results by resulting groups
hc = as.data.frame (groupsHclusterprotein1)
h1 = subset (hc, hc$groupsHclusterprotein1==1)
h2 = subset (hc, hc$groupsHclusterprotein1==2)
h3 = subset (hc, hc$groupsHclusterprotein1==3)
h4 = subset (hc, hc$groupsHclusterprotein1==4)

# Get the names of the proteins by groupd
listHclust = list (h1,h2,h3,h4)
hNames = mapply (rownames, listHclust)


# Compute the Jaccard distance between
# groups obtained by clustering structures
# and groups resulting by PCA analysis
hMatrix = matrix (nrow=4,ncol=4)
for (j in 1:4) {
	structCluster = cNames [[j]]
	for (i in 1:4) {
		hclustCluster = hNames [[i]]
		hMatrix [j,i] = jaccard (structCluster,hclustCluster)
	}
}

cat ("\nHCLUST vs Structural:\n")
cat ("Jaccard distances between structural and hierarchical clustering of PCAs\n")
rownames (hMatrix) <- c("c1","c2","c3","c4")
colnames (hMatrix) <- c("c1'","c2'","c3'","c4'")
print (round (hMatrix,2))





