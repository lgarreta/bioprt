#!/usr/bin/python
import os, sys

#------------------------------------------------------------------
# Given a trayectory it uses a fast clustering algorith to
# reduce the trayectory to only the main representatives.
# INPUT:  An input directory where the trayectory files are located
# OUTPUT: An output directory with the representative files
#------------------------------------------------------------------
USAGE="\
Reduce a trayectory using a fast clustering\n\
USAGE:  pr00_main.py <inputDir>\n"

RMSDTHRESHOLD = 1.5
NCORES	= 4

def main (args):
	if len (args) < 2:
		print USAGE
		sys.exit (1)

	inputDir = args [1]
	outputDirBins = "outbins"
	outputDirRepr = "outrepr"

	print "Parameters: "
	print "\t Input dir: ", inputDir
	print "\t Output dir bins: ", outputDirBins
	print "\t Output dir representatives: ", outputDirRepr
	print "\t RMSD THRESHOLD: ", RMSDTHRESHOLD
	print "\t NUM CORES : ", NCORES
	print "\n"

	# Split full trajectory in bins (blocks of 1000 pdbs)
	cmm ="pr01_createBins.py %s %s" % (inputDir, outputDirBins)
	os.system (cmm)

	# Get Representatives for each bin
	createDir (outputDirRepr)
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
