#!/usr/bin/Rscript

# Reduces a protein folding trajectory by three
# methods: medoids, static and random

library (bio3d)
library (parallel)
library (cluster)

source ("libs/createDir.R")
source ("libs/splitFilesToBins.R")

nCPUS = 4
USAGE = "USAGE: medoids.R <input dir> <output dir> <size chunks>\n"
options (width=400)
tmpDir = "/dev/shm"
#--------------------------------------------------------------
# Main function
#--------------------------------------------------------------
main <- function (args) {
	#args = c("shm/2YCC", "40")
	if (length (args) < 1) {
		cat (USAGE)
		quit ()
	}
	inputPathname  = args [1]
	sizeMedoids    = strtoi (args [2])
	outputDir      = getwd ()
	if (length (args) > 2) 
		outputDir      = args [3]

	cat (">>> Splitting Files into bins...\n")
	binDirList = splitFilesToBins (inputPathname, sizeMedoids, tmpDir)

	#cat (">>> Calculating medoids...\n")
	#medoidsReduction (binDirList, outputDir)

	cat (">>> Static reduction...\n")
	staticReduction (binDirList, outputDir)
}

#--------------------------------------------------------------
# Get the center structure for each bin from the "binDirList"
# and write results to the "outputDir"
#--------------------------------------------------------------
staticReduction <- function (binDirList, outputDir) {
	staticDir =  sprintf ("%s/%s", outputDir, "static")
	createDir (staticDir)
	for (binDir in binDirList) {
		pdbs = list.files (binDir)
		n = length (pdbs)
		centerPdb = pdbs [ceiling (n / 2)]
		pdbFilename = sprintf ("%s/%s", binDir, centerPdb)
		print (sprintf ("cp %s, %s", pdbFilename, staticDir)) 
	}
}

#--------------------------------------------------------------
# Calculate the medoids for each bin of the "binDirList"
# and write results to the "outputDir"
#--------------------------------------------------------------
medoidsReduction <- function (binDirList, outputDir) {
	medoidsList <- mclapply (X=binDirList, FUN=medoidsFromBin, mc.cores=nCPUS )

	medoidsDir =  sprintf ("%s/%s", outputDir, "medoids")
	createDir (medoidsDir)

	for (medoid in medoidsList) {
		system (sprintf ("cp %s %s/%s", medoid, outputDir, "medoids"))
	}
}
#--------------------------------------------------------------
#--------------------------------------------------------------
medoidsFromBin <- function (binDir) {
	outputGroups   = paste (binDir, ".groups", sep="")
	outputMedoids  = paste (binDir, ".medoids", sep="")

	cat ("Medoids from..", binDir, "\n")
	pdbObjects <- getPDBFiles (binDir)

	results = clusteringPDBs (pdbObjects)

	write.table (file=outputGroups, results$groups)
	write.table (file=outputMedoids, results$medoids)

	medoidFilename = list.files (binDir) [results$medoids]
	return (sprintf ("%s/%s", binDir, medoidFilename))
}

#--------------------------------------------------------------
clusteringPDBs <- function (pdbObjects) {
	rmsdDistances <- pairwiseDistancesRMSDs (pdbObjects)

	pamPDBs <- pam (rmsdDistances, 1, diss=T)
	groups  <- pamPDBs$clustering
	medoids <- pamPDBs$medoids

	return (list(groups=groups,medoids=medoids))
}


#--------------------------------------------------------------
# Calculate pairwise RMSDs
#--------------------------------------------------------------
pairwiseDistancesRMSDs <- function (pdbObjects) {
	pdbs = pdbObjects$pdbs
	pdb = pdbs [[1]]
	CAs <- atom.select (pdb, elety="CA", verbose=FALSE)
	firstPdb = pdb$xyz[CAs$xyz]

	n = length (pdbObjects$pdbs)

	# Calculate matrix of coordinates xyz
	xyzMatrix <- matrix (firstPdb,nrow=1)
	for (pdb in pdbs [2:n]) {
		CAs = atom.select (pdb, elety="CA", verbose=FALSE)
		pdb = pdb$xyz[CAs$xyz]
		xyzMatrix = rbind (xyzMatrix, n1=pdb)
	}
	rownames (xyzMatrix) <- names (pdbObjects$pdbs)

	# Calculate RMSDs
	xyz <- fit.xyz (fixed = firstPdb, mobile = xyzMatrix, ncore=nCPUS)

	rmsdDistances <- as.dist (rmsd (xyz, ncore=nCPUS))

	return (rmsdDistances)
}

#--------------------------------------------------------------
# Load pdb files to pdb objects
#--------------------------------------------------------------
getPDBFiles <- function (pathname) {
	# Extracts to an pathname 
	if (grepl ("gz", pathname)==T) {
		stemName = strsplit (pathname, split="[.]")[[1]][1]
		untar (pathname, compressed=T, exdir=stemName)
		inputDir = stemName 
	} else 
		inputDir = pathname

	pdbNames     = list.files (inputDir)
	pdbNamesFull = sapply (pdbNames, function (x) paste (inputDir, x, sep="/"))
	n = length (pdbNamesFull)
	native = pdbNamesFull [[n]]

	# Load PDB Objects
	nativeObject <<- read.pdb2 (native)
	pdbObjects <<- mclapply (X=pdbNamesFull, FUN=read.pdb2, mc.cores=nCPUS )
	
	return (list (n=n, native=nativeObject, pdbs=pdbObjects))
}

#--------------------------------------------------------------
# Print to text file the input object
#--------------------------------------------------------------
log <- function (object, dwfilename) {
	sink (filename)
	print (object)
	sink()
}
#--------------------------------------------------------------
# Call main function
#--------------------------------------------------------------

args = commandArgs (TRUE)
main (args)
