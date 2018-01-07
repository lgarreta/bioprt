#!/usr/bin/Rscript
#!/home/mmartinez/bin/Rscript

# Release 1.0:

#----------------------------------------------------------
# Makes a detailed global clustering of protein conformations 
# from the representatives resulting from the local clustering 
# INPUT:  inputDir filename with the protein conformations
#         outputDir filename to write the results
#
# OUTPUT: Medoids for each bin cluster and their distance matrix
#----------------------------------------------------------
USAGE="USAGE: reduction.R <inputDir> <outputDir> <num cores>\n" 

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
	if (length (args) < 3){
		cat (USAGE)
		quit (status=1)
	}

	INPUTDIR  = args [1] 
	OUTPUTDIR = args [2]
	NCORES    = as.numeric (args [3])

	#createDir (OUTPUTDIR)
	print ("\n FULL")
	print (OUTPUTDIR)
	binPathLst = list.files (INPUTDIR, pattern="bin", full.names=T)
	results=mclapply (binPathLst, reduceGlobal, outputDir=OUTPUTDIR, mc.cores=NCORES)
	#for (inputBinPath in binPathLst) 
	#	reduceGlobal (inputBinPath, outputDir, THRESHOLD)
}

#----------------------------------------------------------
# Reduction function to reduce a single bin
#----------------------------------------------------------
reduceGlobal <- function (inputBinPath, outputDir) {
	cat ("\n>>> Global Reducing ", inputBinPath )
	# Create the output dir for representatives
	clusDir = (paste (getwd(), outputDir, basename (inputBinPath), sep="/"))
	#createDir (clusDir)
	clusDir = inputBinPath

	# Fast clustering for bin, writes representatives to clusDir
	representativePDBPaths <- list.files (inputBinPath, pattern=".pdb", full.names=F)
	cat ("\n>>>>>>>>>>>>>>>>>>>>>>>>>>\n")
	print (inputBinPath)
	print (outputDir)
	print (clusDir)
	print (representativePDBPaths)
	cat ("\n>>>>>>>>>>>>>>>>>>>>>>>>>>")

	# Get Medoid from clustir and write to output dir
	fullClustering (clusDir, outputDir, representativePDBPaths)
	cat ("\n")
}

#----------------------------------------------------------
# Clustering around medoids. Return one medoid for all inputs
#----------------------------------------------------------
fullClustering <- function (inputDir, outputDir, representativePDBPaths) {
	cat ("\n>>> fullClustering...", inputDir, "\n" )
	cat ("\n>>> fullClustering...", outputDir, "\n" )
	if (length (representativePDBPaths) < 2)
		medoid = 1
	else {
		binDir = paste (inputDir, basename (inputDir), sep="/")
		binDir = inputDir
		print (binDir)
		TMscoreDistanceMatrix <<- getTMDistanceMatrix (binDir, outputDir)

		pamPDBs <<- pam (TMscoreDistanceMatrix, k=1, diss=F)
		medoid <<- pamPDBs$id.med
	}
					
	medoidName <<- representativePDBPaths [medoid]
	cmm <<- sprintf ("ln -s %s//%s//%s %s//%s", getwd (), inputDir, medoidName, outputDir, medoidName)   
	cat (paste ("\n", cmm, "\n"))
	system (cmm)
	scan ()

	return (medoidName)
}

#--------------------------------------------------------------
# Calculate pairwise RMSDs
#--------------------------------------------------------------
getTMDistanceMatrix <- function (inputBinDir, outputBinDir) {
	cat ("\nTMD: ", outputBinDir, "\n")
	system (paste ("TMscore-distance-matrix.py", inputBinDir, outputBinDir))
	distanceFilename = paste (outputBinDir, "/", basename (inputBinDir), ".dist", sep="")
	matrixTM = read.table (distanceFilename, header=T)
	distanceMatrixTM =  as.dist (matrixTM)
	return (distanceMatrixTM)
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
	
