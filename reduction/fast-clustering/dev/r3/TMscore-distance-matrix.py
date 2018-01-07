#!/usr/bin/python
"""
Calculate the distance (TM-score) matrix for structures within an input dir.
It wrappers the binary program "TMscore" that takes two PDB files
and prints the results to the screen. This output is preprocessed by extracting
the tm-score value.
INPUT  : 
  <input dir>  The name of the input dir with the protein structures
  <output dir> The name of the output dir
OUTPUT : 
  A file named as the input dir plus the suffix ".dist" containing the 
  matrix of dissimilarites.
"""
USAGE="calculate-dist-matriz-tmscore.py <input dir> <outpud dir>"

import os, sys
import subprocess

def main (args):
	if len (args) < 3:
		print USAGE
		sys.exit (1)

	inputDir  = args [1]
	outputDir = args [2]
	pdbLst    = [inputDir+"/"+x for x in os.listdir (inputDir) if ".pdb" in x] 
	pdbLst.sort()

	n = len (pdbLst)
	matrix = [[""]]
	for i in range (n):
		pdb = pdbLst [i]
		#matrix[0].append (pdb[13:16])
		matrix[0].append (os.path.basename (pdb))
	for i in range (n):
		pdb1 = pdbLst [i]
		#matrix.append ([pdb1[13:16]])
		matrix.append ([os.path.basename (pdb1)])
		for j in range (n):
			pdb2 = pdbLst [j]
			value = 1-getTMScore (pdb1, pdb2)
			matrix [i+1].append (round(value,2))

	#outDistanceFilename = os.path.dirname (inputDir) + "/" + os.path.basename (inputDir) + ".dist"
	outDistanceFilename = outputDir + "/" + os.path.basename (inputDir) + ".dist"
	print "\n>>> Calculating TM score..", inputDir, outDistanceFilename
	sys.stdout = open (outDistanceFilename, "w")
	printMatrix (matrix)

#--------------------------------------------------------------------
#--------------------------------------------------------------------
def printMatrix (matrix):
	n = len (matrix)
	for i in range (n):
		m = len (matrix[i])
		for j in range (m):
			print matrix [i][j],"\t",
		print ""

#--------------------------------------------------------------------
#--------------------------------------------------------------------
def getTMScore (pdb1, pdb2):
	output = runProgram (["TMscore %s %s" % (pdb1, pdb2)])
	value  = preprocesOutputTM (output) 
	return  value
#--------------------------------------------------------------------
#--------------------------------------------------------------------
def preprocesOutputTM (output):
	value = -1
	for line in output.split("\n"):
		if len (line.strip()) == 0:
			continue

		if line.split()[0] == "TM-score":
			value = float (line.split()[2])
			break
	return value

#--------------------------------------------------------------------
# call a external programm that returns a value running on "workingDir"
#--------------------------------------------------------------------
def runProgram (listOfParams):
	try:
		#value = subprocess.Popen (listOfParams, cwd=workingDir, stdout=subprocess.PIPE, shell=True).communicate()[0]
		value = subprocess.Popen (listOfParams, stdout=subprocess.PIPE, shell=True).communicate()[0]
		return value
	except:
		print("runProgram Error:" + str(listOfParams))


main (sys.argv)
