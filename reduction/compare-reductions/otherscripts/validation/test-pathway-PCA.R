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
	#args = c("shm/1FCA1")
	if (length (args) < 1) {
		cat (USAGE)
		quit ()
	}
	inputDir  = args [1]
	outputDir      = getwd ()
	if (length (args) > 2) 
		outputDir      = args [2]

	cat ("\nLoading PDBs...\n")
	pdbObjects <- getPDBFiles (inputDir)
	n = length (pdbObjects$pdbs)

	cat ("\nCalculating RMSDs\n")
	rmsdDistances <- calculateRMSD (pdbObjects)

	mds <- cmdscale (rmsdDistances)

	pdf (file=sprintf ("%s/%s", outputDir,"pca.pdf"))
		colors = rainbow (n,start=.2,end=.1)
		colors2 = heat.colors (n, alpha=1)
		plot (mds, type="p", pch=21, col=colors,bg=colors)
	dev.off ()
}
#--------------------------------------------------------------
# Calculate the RMSD between two protein structures
#--------------------------------------------------------------
calculateRMSD <- function (pdbObjects) {
	n = length (pdbObjects$pdbs)
	pdbRef     <- pdbObjects$pdbs [[n]]
	matRMSDs = matrix (nrow=n, ncol=n)

	ref    <- pdbRef
	refCAs <- atom.select (ref, elety="CA", verbose=FALSE)
	refx   <- ref$xyz [refCAs$xyz]
	for (i in 1:n) {
		trg    <- pdbObjects$pdbs [[i]]
		trgCAs <- atom.select (trg, elety="CA", verbose=FALSE)
		trgx   <- trg$xyz [trgCAs$xyz]
		
		rm       <- rmsd (trgx, refx, fit=T)
		matRMSDs [i,] = rm
	}

	print (round (matRMSDs,2))
	return (as.dist (matRMSDs))
}
#--------------------------------------------------------------
# Calculate pairwise RMSDs
#--------------------------------------------------------------
pairwiseDistancesRMSDs <- function (pdbObjects) {
	pdbs = pdbObjects$pdbs
	n = length (pdbObjects$pdbs)
	pdb = pdbs [[n]]
	CAs <- atom.select (pdb, elety="CA", verbose=FALSE)
	lastPdb = pdb$xyz[CAs$xyz]


	# Calculate matrix of coordinates xyz
	xyzMatrix <- matrix (lastPdb,nrow=1)
	for (pdb in pdbs [2:n]) {
		CAs = atom.select (pdb, elety="CA", verbose=FALSE)
		pdb = pdb$xyz[CAs$xyz]
		print (">>> Dims")
		print (dim (xyzMatrix))
		print (length (pdb))
		xyzMatrix = rbind (xyzMatrix, n1=pdb)
	}
	rownames (xyzMatrix) <- names (pdbObjects$pdbs)

	# Calculate RMSDs
	xyz <- fit.xyz (fixed = lastPdb, mobile = xyzMatrix, ncore=nCPUS)

	rmsdDistances <- as.dist (rmsd (lastPdb, xyz, ncore=nCPUS))

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
	nativeObject <- read.pdb2 (native)
	pdbObjects <- mclapply (X=pdbNamesFull, FUN=read.pdb2, mc.cores=nCPUS )
	
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
