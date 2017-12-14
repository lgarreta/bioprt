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
USAGE="python apply_cluster_by_blocks.py <sub_dir> <max_files_by_dir>"

import os, sys, math
PERCENTAGE = 0.1

#--------------------------------------------------------------------
# Main
#--------------------------------------------------------------------
def main (args):
	if len (args) < 2:
		print USAGE
		sys.exit (0)

	dirName = args [1]
	maxFiles = int (args [2])

	listOfSublists = splitFilesIntoSublists (dirName, maxFiles)
        listOfSubirs   = createSubdirsFromSublists (listOfSublists, dirName+"/blocks")
        executeClusterCommands (listOfSubirs)
        #print listOfSubirs

#----------------------------------------------------------
# Split the files into sublist 
#----------------------------------------------------------
def splitFilesIntoSublists (inputDir, maxFiles):

	filteredFiles = filter (lambda x: ".pdb" in x, os.listdir (inputDir))
	inputFiles = [inputDir+"/"+x for x in filteredFiles]
	nFiles     = len (inputFiles)

	# Split into sublists
	nFilesSublist  = int(math.ceil(float(nFiles) / maxFiles))

	listOfSublists = []
	ini = 0

	for i in  range(nFilesSublist):
		ini = i * maxFiles
		end = (i+1) * maxFiles
		if i == nFilesSublist - 1:
			end = nFiles
		sublist = inputFiles [ini:end]
		listOfSublists.append (sublist)
	return listOfSublists	

#----------------------------------------------------------
# Create subdirectories with the links to the files
#----------------------------------------------------------
def createSubdirsFromSublists (listOfSublists, outputDir):
	os.system ("mkdir -p %s" % outputDir)
	curdir = os.getcwd ()

	listOfSubirs = []
	for n, sublist in enumerate (listOfSublists):
		dirName = "%s/%s/dir%s" % (curdir, outputDir, str (n).zfill (5))
                dirName = "%s/block_%s" % (outputDir, str (n).zfill (5))
		listOfSubirs.append (dirName)
		#print dirName
		cmm = "mkdir -p %s" % dirName
		os.system (cmm)
		cmmList = []
 
		for matchFile in sublist:
			cmm = "ln -s %s/%s %s" % (curdir, matchFile, dirName)
                        cmm = "ln -s %s %s" % ( matchFile, dirName)
			#cmmList.append (cmm)
                        os.system (cmm)

	return listOfSubirs

#----------------------------------------------------------
# Execute remote commands into subdirs
#----------------------------------------------------------
def executeClusterCommands (listOfSubirs):
        command = "./cluster.r %s %s"

	cmmList = []
	for n, dirName in enumerate (listOfSubirs):
		filteredFiles = filter (lambda x: ".pdb" in x, os.listdir (dirName))
	        inputFiles = [dirName+"/"+x for x in filteredFiles]
		
                numGroupsByCluster = int(math.ceil(len(inputFiles) * PERCENTAGE))
                strCmm = command % (dirName, str(numGroupsByCluster))

		print ("El comando a ejecutar finalmente.............>>>", strCmm)
		cmmList.append (strCmm)
                os.system (strCmm)

#----------------------------------------------------------
# Main
#----------------------------------------------------------
if __name__ == "__main__":
    main (sys.argv)
