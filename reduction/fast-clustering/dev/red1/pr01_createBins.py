#!/usr/bin/python

USAGE="\
Split the files of an input directory in bins according to the\n\
size of the bin. The bins are put in an output directory\n\
USAGE: createBins.py <inputDir> <outputDir>\n"

import os, sys, math
SIZEBIN  = 100   # Number of files for each bin
SIZEFILL = 6      # Number of zeros to fill bin names

#--------------------------------------------------------
# Main function to be called from command line
#--------------------------------------------------------
def main (args):
    if len (args) < 3:
        print USAGE
        sys.exit (1)

    inputDir  = "%s/%s" % (os.getcwd (), args [1])
    outputDir = "%s/%s" % (os.getcwd (), args [2])

    createDir (outputDir)

    createBins (inputDir, outputDir, SIZEBIN, SIZEFILL)

#--------------------------------------------------------
# Creates bins from files in an input dir 
# outputDir: destiny dir for bins
# binSize is the number of file by bin
# binFill is the prefix for each bin filename
#--------------------------------------------------------
def createBins (inputDir, outputDir, binSize, binFill):
	inputFiles  = getSortedFilesDir (inputDir, ".pdb")
	binList = splitBins (inputFiles, binSize)

	for k,lst in enumerate (binList):
		binNumber = k+1
		print ">>> Creating bin %s..." % binNumber
		binDirname = "%s/%s%s" % (outputDir, "bin", str (binNumber).zfill (binFill))
		os.mkdir (binDirname)
		for filename in lst:
			sourceFilename  = "%s/%s" % (inputDir, filename)
			destinyFilename = "%s/%s" % (binDirname, filename)
			os.symlink (sourceFilename, destinyFilename)

#--------------------------------------------------------
# Creates a list of sublist where each sublist correspond to
# the files of each bin
#--------------------------------------------------------
def splitBins (inputFiles, binSize):
	nSeqs = len (inputFiles)
	#nBins = nSeqs / binSize
	nBins = int (math.ceil (1.0*nSeqs / binSize))

	binList = []
	for k in range (nBins):
		start = k*binSize
		end   = start + binSize 
		if k < binSize-1:
			binList.append (inputFiles [start:end])
		else:
			binList.append (inputFiles [start:])

	return binList

#------------------------------------------------------------------
# Utility to create a directory safely.
# If it exists it is renamed as old-dir 
#------------------------------------------------------------------
def createDir (dir):
	def checkExistingDir (dir):
		if os.path.lexists (dir):
			headDir, tailDir = os.path.split (dir)
			oldDir = os.path.join (headDir, "old-" + tailDir)
			if os.path.lexists (oldDir):
					checkExistingDir (oldDir)

			os.rename (dir, oldDir)
	checkExistingDir (dir)
	os.system ("mkdir %s" % dir)
#--------------------------------------------------------------------
# Get the files containing the pattern from a inputDir 
#--------------------------------------------------------------------
def getSortedFilesDir (inputDir, pattern=""):
	files  = [x for x in os.listdir (inputDir) if pattern in x ]
	return sorted (files)
#--------------------------------------------------------------------
# Call main with input parameter
#--------------------------------------------------------------------
if __name__ == "__main__":
    main (sys.argv)
