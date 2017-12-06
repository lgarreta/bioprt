#!/usr/bin/Rscript

# Calculate the RMSD of two proteins.
# It fits the second one --coordinte superposition.
# Proteins are entered as command line arguments
# USAGE: rmsd.R <proteinTarget> <proteinReference> [outputFilename]

suppressPackageStartupMessages (library (bio3d))

args <- commandArgs (TRUE)

proteinTarget <- args[1]; 
proteinReference <- args[2]

#--- Obtain the stem filename of the compared protein (proteinReference)
proteinTargetFilename = unlist (strsplit (proteinTarget, "\\."))[1]
rmsdFilename = paste (proteinTargetFilename, "", ".rmsd", sep="")

if (length (args) == 3)
	rmsdFilename <- args [3]

target <- read.pdb2 (proteinTarget, rm.alt=FALSE, verbose=FALSE)
reference <- read.pdb2 (proteinReference, rm.alt=FALSE, verbose=FALSE)
lenAtomsTarget <- length (target$atom[,1])
lenAtomsReference <- length (reference$atom[,1])

residueRange <- NULL

if (length (args) == 4)
	residueRange <- paste (args [3], ":", args[4])

target.ind <- atom.select (target, elety="CA", resno=residueRange, verbose=FALSE)
reference.ind <- atom.select (reference, elety="CA", resno=residueRange, verbose=FALSE)

#--- Calculte the RMSD fitting the two proteins (coordinate superpostion)
r1 = rmsd (target$xyz[target.ind$xyz], reference$xyz[reference.ind$xyz], fit=TRUE)
write (r1, "")

