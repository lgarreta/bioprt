#!/usr/bin/python
import os, sys

"""
 Given a trayectory it uses a fast clustering algorith to
 reduce the trayectory to only the main representatives.
 INPUT:  
   <inputDir>  An input directory with the trayectory files 
   <outputDir> An output directory with the results (representative files)
"""
USAGE  = "\nReduces a trayectory using a fast clustering"
USAGE += "\nUSAGE   : pr00_main.py <inputDir> <outputDir> <RMSD> <SizeBin> <nCores>"
USAGE += "\nExample : pr00_main.py in1000 out1000 1.5 100 4"

# Constant values
RMSDTHRESHOLD = 1.5   # Threshold for RMSD comparison between structures"
NCORES	      = 4     # Number of cores for multiprocessing
SIZEBIN       = 1000  # Number of files for each bin

def main (args):
	if len (args) < 6:
		print USAGE
		sys.exit (1)

	inputDir      = args [1]
	outputDir     = args [2]
	RMSDTHRESHOLD = float (args [3])
	SIZEBIN       = int (args [4])
	NCORES        = int (args [5])

	outputDirBins = "%s/outbins" % outputDir
	outputDirRepr = "%s/outrepr" % outputDir

	createDir (outputDir)

	print "Parameters: "
	print "\t Input dir: ", inputDir
	print "\t Output dir bins: ", outputDirBins
	print "\t Output dir representatives: ", outputDirRepr
	print "\t RMSD THRESHOLD: ", RMSDTHRESHOLD
	print "\t NUM CORES : ", NCORES
	print "\n"

	# Split full trajectory in bins (blocks of 1000 pdbs)
	cmm ="pr01_createBins.py %s %s %s" % (inputDir, outputDirBins, SIZEBIN)
	os.system (cmm) 

	# Get Representatives for each bin
	cmm = "pr02_reduction.R %s %s %s %s" % (outputDirBins, outputDirRepr, RMSDTHRESHOLD, NCORES)
	os.system (cmm)

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

#------------------------------------------------------------------
# Call main with input parameter
#------------------------------------------------------------------
if __name__ == "__main__":
	main (sys.argv)
