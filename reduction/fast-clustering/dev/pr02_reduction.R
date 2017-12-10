#!/usr/bin/Rscript
#!/home/mmartinez/bin/Rscript

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

THRESHOLD = 1.3
NCORES    = 1
#----------------------------------------------------------
# Main function
#----------------------------------------------------------
main <- function () {
	args <- commandArgs (TRUE)
	#args = c("outbins", "res", "1.3", "1")
	if (length (args) < 4){
		cat (USAGE)
		quit (status=1)
	}

	inputDir  = args [1] 
	outputDir = args [2]
	THRESHOLD = as.numeric (args [3])
	NCORES    = as.numeric (args [4])

	# Process each bin from input dir
	binPathsLst = list.files (inputDir, pattern="bin", full.names=T)
	for (binPath in binPathsLst) {
		cat ("\n>>> Reducing ", binPath )
		# Create the output dir for representatives
		clusDir = (paste (getwd(), outputDir, basename (binPath), sep="/"))
		createDir (clusDir)

		# Fast clustering for bin, writes representatives to clusDir
		partialClustering (binPath, clusDir)

		# Get Medoid from clustir and write to output dir
		fullClustering (clusDir, outputDir)
	}
	cat ("\n")
}

#----------------------------------------------------------
# Fast clustering following hobbohm algorith
# Write the links to the representatives in the output dir
# Return a list with the representative as the first pdb in the group
#----------------------------------------------------------
partialClustering <- function (inputDir, outputDir) {
	cat (">>> Partial clustering...")
	# Get proteins from input dir
	inputProteinsLst   = list.files (inputDir, pattern=".pdb", full.names=T)

	# Assign the first protein as first group
	lstGroups = list (inputProteinsLst[[1]])

	# Assign the other proteins to groups
	for (protein in inputProteinsLst[2:length(inputProteinsLst)]) {
		groupNumber = getGroupNumberByRmsd (protein, lstGroups, THRESHOLD)

		if (groupNumber != -1) {
			group = lstGroups [[groupNumber]]
			lstGroups [[groupNumber]] = append (lstGroups[[groupNumber]],protein)
		} else {
			group = list (protein)
			cmm = sprintf ("ln -s %s/%s %s/%s", getwd(), protein, outputDir, basename(protein))
			#print (cmm)
			system (cmm)
			lstGroups = append (lstGroups, list (group))
		}
	} 
	return (lstGroups)
}	

#----------------------------------------------------------
# Clustering around medoids. Return one medoid for all inputs
#----------------------------------------------------------
fullClustering <- function (inputDir, outputDir) {
	cat (">>> fullClustering...")
	pdbNames = list.files (inputDir, full.names=T)
	if (length (pdbNames) == 1)
		medoid = 1
	else {
		#rmsdDistanceMatrix   <<- getRmsdDistanceMatrixWithAlginments (inputDir)
		rmsdDistanceMatrix   <<- getRmsdDistanceMatrix (inputDir)

		pamPDBs <- pam (rmsdDistanceMatrix, 1, diss=T)
		medoid = pamPDBs$medoid
	}

	medoidPath <- pdbNames [medoid]
	cmm = sprintf ("ln -s %s %s/%s", medoidPath, outputDir, basename (medoidPath))   
	system (cmm)

	return (medoidPath)
}

#----------------------------------------------------------
# calculates the protein closest group 
#----------------------------------------------------------
getGroupNumberByRmsd <- function (protein, lstGroups, threshold) {		
	counter = 1
	for (group in lstGroups) {
		localProtein = group [[1]]
		r = calculateRMSD (protein, localProtein)
		if (r < threshold)
			return (counter)
		counter = counter + 1	
	}	
	return (-1)
}
#----------------------------------------------------------
# calculate the rmsd of two proteins.
#----------------------------------------------------------
calculateRMSD <- function (proteinTarget, proteinReference) {
	target <- read.pdb2 (proteinTarget, rm.alt=FALSE, verbose=FALSE)
	reference <- read.pdb2 (proteinReference, rm.alt=FALSE, verbose=FALSE)
	value = rmsd (target$xyz, reference$xyz, fit=TRUE)
	return (value)
}

#----------------------------------------------------------
# Write the groups' representatives to the ouput dir
#----------------------------------------------------------
writeRepresentatives <- function (lstGroups, outputDir) {
	for (g in lstGroups) {
	  	pdb = g[[1]]
		cmm = sprintf ("ln -s %s %s/%s", pdb, outputDir, basename (pdb))
		#print (cmm)
		system (cmm)
	}
}
#----------------------------------------------------------
# Write the groups formatted as lines of text
#----------------------------------------------------------
writeGroups <- function (lstGroups, outputDir=NULL) {
  cat (">>> Writing representatives...\n")
  outFile = paste (outputDir, "representatives.txt", sep="/")
  sink (outFile)
	for (g in lstGroups) {
	  cat (g[[1]])
	  cat ("\n")
	}
  sink ()

  cat (">>> Writing groups...\n")
  outFile = paste (outputDir, "groups.txt", sep="/")
  sink (outFile)
	for (g in lstGroups) {
	  cat (unlist (g))
	  cat ("\n")
	}
  sink ()
  cat (">>> End writing \n")
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
# Calculate pairwise RMSDs
#--------------------------------------------------------------
getRmsdDistanceMatrix <- function (inputDir) {
	pdbObjects = getPDBFiles (inputDir)

	pdb <<- pdbObjects [[1]]
	CAs <- atom.select (pdb, elety="CA", verbose=FALSE)
	firstPdb <<- pdb$xyz[CAs$xyz]

	n <<- length (pdbObjects)

	# Calculate matrix of coordinates xyz
	xyzMatrix <- matrix (firstPdb,nrow=1)
	for (pdb in pdbObjects [2:n]) {
		CAs = atom.select (pdb, elety="CA", verbose=FALSE)
		pdb = pdb$xyz[CAs$xyz]
		xyzMatrix = rbind (xyzMatrix, n1=pdb)
	}
	rownames (xyzMatrix) <- names (pdbObjects)
	mat <<-xyzMatrix

	# Calculate RMSDs
	xyz <- fit.xyz (fixed = firstPdb, mobile = xyzMatrix, ncore=1)

	rmsdDistances <<- as.dist (rmsd (xyz, ncore=1), diag=T)

	return (rmsdDistances)
}

#--------------------------------------------------------------
# Load pdb files to pdb objects, makes firstly an alignment
# and return a distance matrix
#--------------------------------------------------------------
getRmsdDistanceMatrixWithAlginments <- function (inputDir) {
	pdbNames     = list.files (inputDir, full.names=T)
	pdbObjects     <<- pdbaln (pdbNames, NCORES)
	rmsdDistances  <<- rmsd (pdbObjects, fit=T)
}

#--------------------------------------------------------------
# Load pdb files to pdb objects 
#--------------------------------------------------------------
getPDBFiles <- function (inputDir) {
	pdbNames     = list.files (inputDir, full.names=T)

	# Load PDB Objects
	pdbObjects <<- mclapply (X=pdbNames, FUN=read.pdb2, mc.cores=NCORES)
	return (pdbObjects)
}

#--------------------------------------------------------------
#--------------------------------------------------------------
main () 
	
