#!/usr/bin/python

"""
  r0.8
  Executes distributively a command over N servers using remote ssh.
  - It takes as inputs the dir of the input files and the dir where it
  distributes (links) the files in N dirs for the N servers.
  - Then it connects to each server with ssh and execute the command.
  - It assumes that all the servers are in a common filesystem so the
  output subdirectories are common for all servers. (copied lgarreta)
""" 

USAGE="python run_multicore.py <matches input dir> <clustering output dir>"

import os, sys, time, shutil
from multiprocessing import Pool 

#--------------------------------------------------------------------
# Main
#--------------------------------------------------------------------
def main (args):
	if len (args) < 3:
		print USAGE
		sys.exit (0)

	inputDir = args [1]
	outputDir = args [2]
	nCpus = int (args [3])

        shutil.rmtree(outputDir, ignore_errors=True) # Nota: Se pueden borrar resultados importantes de otras ejecuciones

	# Name of servers with their number of cores
        now = time.strftime("%c")
        print ("Current time (ini) %s"  % now )
	listOfSublists = splitFilesIntoSublists (inputDir, nCpus)
	listOfSubirs   = createSubdirsFromSublists (listOfSublists, outputDir, nCpus)
	executeClusterCommands (outputDir, listOfSubirs, nCpus)
	print "End---------------------------"
        now = time.strftime("%c")
        print ("Current time (end) %s"  % now )

#----------------------------------------------------------
# Split the files into sublist 
#----------------------------------------------------------
def splitFilesIntoSublists (inputDir, nCpus):
	nPCs = nCpus
	filteredFiles = filter (lambda x: ".pdb" in x, os.listdir (inputDir))
        #print "filteredFiles"
	inputFiles = [inputDir+"/"+x for x in filteredFiles]
	nFiles     = len (inputFiles)

	# Split into sublists
	nFilesSublist  = nFiles / nPCs
	#print "nFiles " + str(nFiles) + " nPCs " + str(nPCs)
        #print "Sublists size: " + str(nFilesSublist)
	listOfSublists = []
	ini = 0

	for i in  range(nCpus):
		ini = i * nFilesSublist
		end = (i+1) * nFilesSublist
		if i == nPCs - 1:
			end = nFiles

		sublist = inputFiles [ini:end]
		listOfSublists.append (sublist)
	#print listOfSublists
	return listOfSublists	

#----------------------------------------------------------
# Create subdirectories with the links to the files
#----------------------------------------------------------
def createSubdirsFromSublists (listOfSublists, outputDir, nCpus):
	os.system ("mkdir -p %s" % outputDir)
	curdir = os.getcwd ()
	listOfSubirs = []
	for n, sublist in enumerate (listOfSublists):
		dirName = "%s/%s/proc_%s" % (curdir, outputDir, str (n).zfill (5))
		listOfSubirs.append (dirName)
		#print dirName
		cmm = "mkdir %s" % dirName
		os.system (cmm)
		cmmList = []
 
		for matchFile in sublist:
			cmm = "ln -s %s/%s %s" % (curdir, matchFile, dirName)
			cmmList.append (cmm)
		#print "---------------------------------------------------------------------------------------------"                
		pool = Pool (processes=nCpus)	
		pool.map (os.system, cmmList)
		pool.close ()
		pool.join ()	

	return listOfSubirs

#----------------------------------------------------------
# Execute remote commands into subdirs
#----------------------------------------------------------
def executeClusterCommands (outputDir, listOfSubirs, nCpus):
        command = "python apply_cluster_by_blocks.py %s 1000"

	cmmList = []
	for n, dirName in enumerate (listOfSubirs):	
		#pcName, pcUnits = listOfPCs[n][0], listOfPCs[n][1]
		#strCmm = command % (dirName, "/dev/shm/out",  pcUnits)
                strCmm = command % (dirName)

		#print (">>>", strCmm)
		cmmList.append (strCmm)

	pool = Pool (processes=nCpus)	
	pool.map (os.system, cmmList)
	pool.close ()
	pool.join ()
	

	#input (">>>")		
		
#----------------------------------------------------------
# Main
#----------------------------------------------------------
if __name__ == "__main__":
    main (sys.argv)

