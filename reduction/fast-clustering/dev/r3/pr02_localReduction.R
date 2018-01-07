#!/usr/bin/Rscript
#!/home/mmartinez/bin/Rscript

# Release 1.0:

#----------------------------------------------------------
# Make a fast clustering of protein conformations from a 
# trayectory by doing first a fast local clustering 
# INPUT:  
#   <inputDir>       An input directory with the trayectory files 
#   <outputDir>      An output directory with the results (representative files)
#   <RMSD threshold> Threshold for local reduction comparisons that uses RMSD
#   <Size Bin>       Number of structures for each bin of the partiioned trajectory
#   <N Cores>        Number of cores to use for paralelizing the procedure
#
# OUTPUT:            An output dir with the clustering results for each bin
#----------------------------------------------------------
USAGE="USAGE: reduction.R <inputDir> <outputDir> <RMSD Threshold> <num cores>\n" 

library (bio3d)
library (parallel)
library (cluster)
options (width=300)

#THRESHOLD = 1.3
#NCORES= 1
#----------------------------------------------------------
# Main function
#----------------------------------------------------------
main <- function () {
	args <- commandArgs (TRUE)
	#args = c("io/out1000/outbins", "io/out1000/outrepr", "1.5", "1")
	print (args)
	if (length (args) < 4){
		cat (USAGE)
		quit (status=1)
	}

	INPUTDIR  = args [1] 
	OUTPUTDIR = args [2]
	THRESHOLD = as.numeric (args [3])
	NCORES    = as.numeric (args [4])

	createDir (OUTPUTDIR)

	binPathLst = list.files (INPUTDIR, pattern="bin", full.names=T)
	results=mclapply (binPathLst, reduceLocal, outputDir=OUTPUTDIR, threshold=THRESHOLD, mc.cores=NCORES)
	print (binPathLst)
	#for (binPath in binPathLst) 
	#	reduceLocal (binPath, outputDir, THRESHOLD)
}

#----------------------------------------------------------
# Reduction function to reduce a single bin
#----------------------------------------------------------
reduceLocal <- function (inputBinPath, outputDir, threshold) {
	cat ("\n>>> Local Reducing ", inputBinPath )
	# Create the output dir for representatives
	outputBinPath = (paste (getwd(), outputDir, basename (inputBinPath), sep="/"))
	createDir (outputBinPath)		

	pdbsBinLst <<- getPDBFiles (inputBinPath)
	# Fast clustering for bin, writes representatives to outputBinPath
	localClustering (inputBinPath, outputBinPath, threshold, pdbsBinLst)
}

#----------------------------------------------------------
# Fast clustering following hobbohm algorith
# Write the links to the representatives in the output dir
# Return a list with the representative as the first pdb in the group
#----------------------------------------------------------
localClustering <- function (inputBinPath, outputBinPath, threshold, pdbsBinLst) {
	# Get proteins from input dir
	inputProteinsLst   = paste (inputBinPath, rownames (pdbsBinLst), sep="/")

	representativePDBsLst = c ()
	lstGroups = list ()
	n = length (inputProteinsLst)
	for (k in 1:n) {
		protein = basename (inputProteinsLst[[k]]) 
		groupNumber = getGroupNumberByRmsd (protein, lstGroups, threshold, pdbsBinLst)

		if (groupNumber == -1) {
			group = list (protein)
			cmm = sprintf ("ln -s %s/%s/%s %s/%s", getwd(), inputBinPath, protein, outputBinPath, basename(protein))
			cat ("\n>>> ",cmm,"\n")
			system (cmm)
			representativePDBsLst = append (representativePDBsLst, basename (protein))
			lstGroups = append (lstGroups, list (group), after=0)
		} else 
			lstGroups [[groupNumber]] = append (lstGroups[[groupNumber]],protein)
	} 
	return (representativePDBsLst)
}	

#----------------------------------------------------------
# calculates the protein closest group 
#----------------------------------------------------------
getGroupNumberByRmsd <- function (protein, lstGroups, threshold, pdbsBinLst) {
	groupNumber = -1
	counter = 1
	for (group in lstGroups) {
		localProtein       = group [[1]]
		localPdbObject     <<- pdbsBinLst [basename (localProtein),]
		referencePdbObject <<- pdbsBinLst [basename (protein),]

		referencePdbObject = fit.xyz (localPdbObject, referencePdbObject)
		rmsdValue = rmsd (localPdbObject, referencePdbObject, fit=F)
		#rmsdValue2 = rmsdPartial (localPdbObject, referencePdbObject, 0.8, 10)
		#cat ("\nRMSD:", rmsdValue, " > ", localProtein, " GROUP: ", protein)

		if (rmsdValue < threshold) {
			groupNumber = counter
			break
		}
		counter = counter + 1	
	}	
	return (groupNumber)
}

#--------------------------------------------------------------
# Load pdb files to pdb objects 
#--------------------------------------------------------------
getPDBFiles <- function (inputDir) {
	cat ("\n>>> Loading PDB files to objects from: ", inputDir, "\n")
	pdbNames <- list.files (inputDir, full.names=T)

	# Load PDB Objects
	pdbObjects <- mclapply (X=pdbNames, FUN=read.pdb2, mc.cores=1)
	names (pdbObjects) <- basename (pdbNames)
	pdbObjects <- preparePdbObjects (pdbObjects)
	return (pdbObjects)
}

#--------------------------------------------------------------
#--------------------------------------------------------------
preparePdbObjects <- function (pdbObjects) {
	n <- length (pdbObjects)

	# Calculate matrix of coordinates xyz
	for (k in 1:n) {
		pdb =  pdbObjects [[k]] 
		CAs = atom.select (pdb, elety="CA", verbose=FALSE)
		pdb = pdb$xyz[CAs$xyz]
		lastPdb = pdb

		if (k==1) xyzMatrix = rbind (pdb)
		else      xyzMatrix = rbind (xyzMatrix, n1=pdb)

	}
	rownames (xyzMatrix) <- names (pdbObjects)

	# Calculate RMSDs
	xyzMatriz <- fit.xyz (fixed = lastPdb, mobile = xyzMatrix, ncore=1)
	rownames (xyzMatrix) <- names (pdbObjects)

	return (xyzMatrix)
}

#----------------------------------------------------------
# Create dir, if it exists the it is renamed old-XXX
#----------------------------------------------------------
createDir <- function (newDir) {
	checkOldDir <- function (newDir) {
		name  = basename (newDir)
		path  = dirname  (newDir)
		if (dir.exists (newDir) == T) {
			oldDir = sprintf ("%s/old-%s", path, name)
			if (dir.exists (oldDir) == T) {
				checkOldDir (oldDir)
			}

			file.rename (newDir, oldDir)
		}
	}

	checkOldDir (newDir)
	system (sprintf ("mkdir %s", newDir))
}

#--------------------------------------------------------------
#--------------------------------------------------------------
main () 
	
