#!/usr/bin/Rscript

library (bio3d)
library (parallel)

MC_CORES=4

# Creates a .pdf plot for a protein trajectory o pathway
# The input is either a compressed (.tgz) trajectory or
# a directory name with the PDBs files inside it.
USAGE="plot-pathway.R <input pathway dir> [output dir]\n"
#--------------------------------------------------------------
# Main function
#--------------------------------------------------------------
main <- function () {
	args = commandArgs (TRUE)
	if (length (args) < 1) {
		cat (USAGE)
		quit ()
	}
	pathname  = args [1]
	outputDir = getwd ()
	if (length (args) == 2)
		outputDir = args [2]

	filenames = getInOutNames (pathname, outputDir)
	pathname = filenames$inputDir
	outputFile = filenames$outputFile

	# Extract or load filename to calculate RMSDs
	cat ("\nLoading PDBs...\n")
	files = getPDBFiles (pathname)
	
	cat ("\nCalculating RMSDs...")	
	rmsdValues = parallerRmsdPathway (files$n, files$native, files$pdbs)

	cat ("\nWriting output file ", outputFile, "\n")
	plotPathway (rmsdValues, outputFile)
}

#--------------------------------------------------------------
# Calculate the RMSD between two protein structures
#--------------------------------------------------------------
parallerRmsdPathway <- function (n, native, pdbs) {
	rmsdValues = mclapply (pdbs, calculateRMSD, native, 
						   mc.preschedule=T, mc.set.seed=T, mc.cores=30)
	return (rmsdValues)
}
	
rmsdPathway <- function (n, native, pdbs) {
	rmsdValues = c()
	for (i in 1:n) {
		rmsd = calculateRMSD (pdbs [i], native)
		rmsdValues = append (rmsdValues, rmsd)
	}
	return (rmsdValues)
}

#--------------------------------------------------------------
# Calculate the RMSD between two protein structures
#--------------------------------------------------------------
calculateRMSD <- function (pdbNameTarget, pdbNameRef) {
	#--- Obtain the stem pathname of the compared protein (pdbNameRef)
	proteinTargetFilename = unlist (strsplit (pdbNameTarget, "\\."))[1]

	target <- read.pdb (pdbNameTarget, rm.alt=FALSE, verbose=FALSE)
	reference <- read.pdb (pdbNameRef, rm.alt=FALSE, verbose=FALSE)

	targetCAs <- atom.select (target, elety="CA", verbose=FALSE)
	referenceCAs <- atom.select (reference, elety="CA", verbose=FALSE)

	#--- Calculte the RMSD fitting the two proteins (coordinate superpostion)
	r1 = rmsd (target$xyz[targetCAs$xyz], reference$xyz[referenceCAs$xyz], fit=TRUE)

	return  (r1)
}

#-------------------------------------------------------------
# Creata a XY plot from the RMSD values of each conformation
#-------------------------------------------------------------
plotPathway <- function (rmsdValues, outputFile) {
	pdf (outputFile, width=14)
		n = length(rmsdValues)
		rd = rmsdValues[1:(n-1)]
		time = 1:(n-1)
		plot(time, rd, typ = "l", ylab = "RMSD", xlab = "Frame No.")
		points (lowess(time,rd, f=2/10), typ="l", col="red", lty=2, lwd=2)
		#steps = n / 21
		#xPoints = seq (0,n, ceiling (steps))
		#axis (side=1, xPoints)
	dev.off ()
}

#--------------------------------------------------------------
# Get the PDB files from either a compressed file or a dir
#--------------------------------------------------------------
getPDBFiles <- function (pathname) {
	# Extracts to an inputDir 
	if (grepl ("gz", pathname)==T) {
		stemName = strsplit (pathname, split="[.]")[[1]][1]
		inputDir = stemName 
		untar (pathname, compressed=T, exdir=inputDir)
	} else 
		inputDir = pathname

	inputFiles = list.files (inputDir)
	inputFilesFull = sapply (inputFiles, function (x) paste (inputDir, x, sep="/"))
	n = length (inputFilesFull)
	native = inputFilesFull [[n]]
	outfile = paste (inputDir, ".pdf", sep="")

	return (list (n=n, native=native, pdbs=inputFilesFull, outfile=outfile))
}

#--------------------------------------------------------------
#--------------------------------------------------------------
getInOutNames <- function (pathname, outputDir) {
	if (grepl ("gz", pathname)==T) 
		inputDir = strsplit (pathname, split="[.]")[[1]][1]
	else  
		inputDir = pathname

	outputFile = sprintf ("%s/%s.pdf", outputDir, basename (inputDir))
	return (list (inputDir=inputDir, outputFile=outputFile))
}
#--------------------------------------------------------------
# Call main function
#--------------------------------------------------------------
main ()
