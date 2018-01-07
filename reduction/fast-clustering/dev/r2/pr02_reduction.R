#!/usr/bin/Rscript
#!/home/mmartinez/bin/Rscript

# Release 4.0:
#  - RMSD for partial clustering
#  - TMscore for full clustering
#  - Command line parameters: inDir, outDir, RMSD, cores 
#  - Uses a external program (TMscore)

# Release 3.0
#  - Load all PDB files instead bin by bin
#  - Create groups in memory instead in disk

#----------------------------------------------------------
# Make a fast clustering of protein conformations from a 
# trayectory by doing first a fast clustering, and then 
# doing a detailed clustering with the representatives 
# of the groups resulting in the fast clustering
# Input: inputDir filename with the protein conformations
#        outputDir filename to write the results
# Output: Two files: "groups.txt" and "representatives.txt"
#----------------------------------------------------------
USAGE="USAGE: reduction.R <inputDir> <outputDir> <RMSD Threshod> <num cores>\n" 

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
	results=mclapply (binPathLst, reduceSingleBin, outputDir=OUTPUTDIR, threshold=THRESHOLD, mc.cores=NCORES)
	#for (binPath in binPathLst) 
	#	reduceSingleBin (binPath, outputDir, THRESHOLD)
}

#----------------------------------------------------------
# Reduction function to reduce a single bin
#----------------------------------------------------------
reduceSingleBin <- function (binPath, outputDir, threshold) {
	cat ("\n>>> Reducing ", binPath )
	# Create the output dir for representatives
	clusDir = (paste (getwd(), outputDir, basename (binPath), sep="/"))
	createDir (clusDir)		

	pdbsBinLst <<- getPDBFiles (binPath)
	# Fast clustering for bin, writes representatives to clusDir
	representativePDBNames <<- partialClustering (binPath, clusDir, threshold, pdbsBinLst)

	# Get Medoid from clustir and write to output dir
	#fullClustering (clusDir, outputDir, representativePDBNames, pdbsBinLst)
	fullClustering (clusDir, outputDir, representativePDBNames, pdbsBinLst)
	cat ("\n")
}

	
#----------------------------------------------------------
# Clustering around medoids. Return one medoid for all inputs
#----------------------------------------------------------
fullClustering <- function (inputDir, outputDir, representativePDBNames, pdbsBinLst) {
	cat ("\n>>> fullClustering...", inputDir )
	if (length (representativePDBNames) < 2)
		medoid = 1
	else {
		#rmsdDistanceMatrix   <<- calcRmsdDistMatrix (representativePDBNames, pdbsBinLst)
		binDir = paste (outputDir, basename (inputDir), sep="/")
		rmsdDistanceMatrix <<- getTMDistanceMatrix (binDir, pdbsBinLst)

		pamPDBs <<- pam (rmsdDistanceMatrix, k=1, diss=F)
		medoid <<- pamPDBs$id.med
	}
					
	medoidPath <<- representativePDBNames [medoid]
	cmm <<- sprintf ("ln -s %s/%s %s/%s", inputDir, medoidPath, outputDir, medoidPath)   
	system (cmm)
	scan ()

	return (medoidPath)
}

#--------------------------------------------------------------
# Calculate pairwise RMSDs
#--------------------------------------------------------------
getTMDistanceMatrix <- function (binDir, pdbsBinLst) {
	system (paste ("TMscore-distance-matrix.py", binDir))
	distanceFilename = paste (binDir, ".dist", sep="")
	matrixTM = read.table (distanceFilename, header=T)
	distanceMatrixTM =  as.dist (matrixTM)
	return (distanceMatrixTM)
}

#--------------------------------------------------------------
# Calculate pairwise RMSDs
#--------------------------------------------------------------
calcRmsdDistMatrix <- function (pdbNames, pdbsBinLst) {
	pdbObjects = pdbsBinLst [pdbNames,]
	n <- length (pdbObjects)

	rmsdDistances <- as.dist (rmsd (pdbObjects, ncore=1), diag=T)

	return (rmsdDistances)
}

#----------------------------------------------------------
# Fast clustering following hobbohm algorith
# Write the links to the representatives in the output dir
# Return a list with the representative as the first pdb in the group
#----------------------------------------------------------
partialClustering <- function (inputDir, outputDir, threshold, pdbsBinLst) {
	cat ("\n>>> Partial clustering... ", inputDir)
	# Get proteins from input dir
	inputProteinsLst   = paste (inputDir, rownames (pdbsBinLst), sep="/")
	representativePDBsLst = c ()
	lstGroups = list ()
	n = length (inputProteinsLst)
	for (k in 1:n) {
		protein = basename (inputProteinsLst[[k]]) 
		groupNumber = getGroupNumberByRmsd (protein, lstGroups, threshold, pdbsBinLst)

		if (groupNumber == -1) {
			group = list (protein)
			cmm = sprintf ("ln -s %s/%s/%s %s/%s", getwd(), inputDir, protein, outputDir, basename(protein))
			#cat ("\n>>> ",cmm,"\n")
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

#----------------------------------------------------------
# Evaluates the RMSD for each atom positions until it reaches
# a min number of atoms (minAminos) with RMSD below a threshold 
#----------------------------------------------------------
rmsdPartial <- function (loc,ref, rmsdThreshold, minAminos) {
	print ("-------------------------------------------------------")
	n = length (ref)
	k = 0
	rmsdAvr = 0
	for (i in seq(1,n,3)) {
		refBlock  = ref [i: (i+2)]
		locBlock  = loc [i: (i+2)]
		#refBlock  = ref 
		#locBlock  = loc 

		rmsdValue = rmsd (refBlock, locBlock, fit=F)
		cat (">>>> ", rmsdValue, " ", k, "\n")
		cat ("loc: "); print (locBlock) 
		cat ("ref: "); print (refBlock)

		if (rmsdValue < 20) {
			k = k + 1
			rmsdAvr = rmsdAvr + rmsdValue
		}
		if (k >= minAminos)
			break;
	}
	return (rmsdAvr/(k))
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
	
