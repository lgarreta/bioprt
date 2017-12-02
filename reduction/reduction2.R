#!/usr/bin/Rscript
#!/home/mmartinez/bin/Rscript

#----------------------------------------------------------
# It clusters protein conformations from a trayectory by
# doing first a fast clustering, and then doing a detailed 
# clustering with the representatives of the groups resulting
# in the fast clustering

# Input: inputDir filename with the protein conformations
#        outputDir filename to write the results
# Output: Two files: "groups.txt" and "representatives.txt"
#----------------------------------------------------------

suppressPackageStartupMessages (library (bio3d))
THRESHOLD = 1.3
options (width=300)

#----------------------------------------------------------
# Main function
#----------------------------------------------------------
main <- function () {

	#args <- commandArgs (TRUE)
	args <- c("in", "out")

	inputDir = args [1]
	outputDir = args [2]

	# Get proteins from input dir
	inputProteinsLst   = paste (inputDir, list.files (inputDir, pattern=".pdb"), sep="/")

	# Assign the first protein as first group
	lstGroups = list (inputProteinsLst[[1]])

	# Assign the other proteins to groups
	for (protein in inputProteinsLst[2:length(inputProteinsLst)]) {
		groupNumber = getGroupNumberByRmsd (protein, lstGroups, THRESHOLD)
		cat (paste (">>> Asigning group number ", groupNumber, "to", protein, "\n"))

		if (groupNumber != -1) {
			group = lstGroups [[groupNumber]]
			lstGroups [[groupNumber]] = append (lstGroups[[groupNumber]],protein)
		} else {
			group = list (protein)
			lstGroups = append (lstGroups, list (group))
		}
	} 

	cat ("\n")
	writeGroups (lstGroups, outputDir)
}	
	
#----------------------------------------------------------
# Write the groups formatted as lines of text
#----------------------------------------------------------
writeGroups <- function (lstGroups, outDir=NULL) {
  cat (">>> Writing representatives...\n")
  outFile = paste (outDir, "representatives.txt", sep="/")
  sink (outFile)
	for (g in lstGroups) {
	  cat (g[[1]])
	  cat ("\n")
	}
  sink ()

  cat (">>> Writing groups...\n")
  outFile = paste (outDir, "groups.txt", sep="/")
  sink (outFile)
	for (g in lstGroups) {
	  cat (unlist (g))
	  cat ("\n")
	}
  sink ()
  cat (">>> End writing \n")
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

#---------------------------------------------------------
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

###############################################################################
main () 
	
