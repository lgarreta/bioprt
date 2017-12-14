#!/home/mmartinez/bin/Rscript

# Takes a dir with proteins and calculates a cluster over them, taking representative proteins.
library (bio3d)

args <- commandArgs (TRUE)
proteinsDir <- args[1];

# We define the number of groups that we want to retrieve
numGroups <- args[2];

files = list.files(proteinsDir)
pdb2files=paste (proteinsDir, "/",files, sep="")

pdbs = pdbaln (pdb2files)

core = core.find (pdbs)
#plot (core)

col = rep("black", length(core$volume))
col[core$volume < 2] = "pink"
col[core$volume < 1] = "red"

xyz = pdbfit (pdbs)
rd = rmsd (xyz)

hist(rd, breaks = 20, xlab = "RMSD (Å)")
hc.rd <- hclust(as.dist(rd))

#grps <- cutree(hc.rd, h = numGroups)
grps <- cutree(hc.rd, k = numGroups)
###write(t(grps), file = paste ("groups", numGroups, sep="_"))
###grps 

###pdf ( paste ("clust", numGroups, sep="_"), width=20)
hc <- hclust(as.dist(rd))
##- draw dendrogram 
####hclustplot(hc, k=2)
###dev.off()


###pdf ("histogram.pdf", width=20)
#plot(hc.rd, labels = pdbs$id, ylab = "RMSD", main = "RMSD Cluster Dendrogram")
###dev.off()


# Get a list with the groups of proteins

# First, combine the information from pdf files with their groups
grpsPdbs = cbind (pdb2files, group=grps)


# Then convert to dataframe to be able to have subsets (queries)
frgrp = as.data.frame (grpsPdbs)

# We name the columns from the dataframe
colnames (frgrp) = c ("pdb", "grp")

# The groups
lslvl = levels(frgrp$"grp")

i = 1
classifByGroups = list()
for(l in lslvl) {
    gi = frgrp [ which (frgrp$"grp"==l),]
    classifByGroups[[i]] = gi[1,]
    i = i + 1
}

##print(classifByGroups)

#Create folder to set the results
clusterOutputDirName = "output_cluster"
if (!file.exists(clusterOutputDirName))
{
	dir.create(clusterOutputDirName)
}

processorName = basename(dirname(dirname(proteinsDir)))
blockName = basename(proteinsDir)
outputFileName = paste(processorName, blockName, sep="_")

classifByGroups.df = as.data.frame(do.call(cbind, classifByGroups))

write(t(classifByGroups.df), file = paste (clusterOutputDirName, "/", outputFileName, ".output", sep=""))


