#--------------------------------------------------------------
# Load pdb files to pdb objects, makes firstly an alignment
# and return a distance matrix
#--------------------------------------------------------------
getRmsdDistanceMatrixWithAlginments <- function (inputDir) {
	pdbNames     = list.files (inputDir, full.names=T)
	pdbObjects     <- pdbaln (pdbNames, NCORES)
	rmsdDistances  <- rmsd (pdbObjects, fit=T)
}
