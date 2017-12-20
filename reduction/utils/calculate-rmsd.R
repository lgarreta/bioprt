
#----------------------------------------------------------
# calculate the rmsd of two proteins.
#----------------------------------------------------------
calculateRMSD <- function (proteinTarget, proteinReference, pdbsBin) {
	target <- read.pdb2 (proteinTarget, rm.alt=FALSE, verbose=FALSE)
	reference <- read.pdb2 (proteinReference, rm.alt=FALSE, verbose=FALSE)
	value = rmsd (target$xyz, reference$xyz, fit=TRUE)
	return (value)
}
