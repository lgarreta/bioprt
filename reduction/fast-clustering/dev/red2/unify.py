#!/usr/bin/python

import os, sys

#--------------------------------------------------------------------
# Main
#--------------------------------------------------------------------
def main (args):
	loadFilesToProcess()


#----------------------------------------------------------
# Split the files into sublist 
#----------------------------------------------------------
def loadFilesToProcess():
	filteredFiles = filter (lambda x: ".output" in x, (os.listdir (os.getcwd () + "/output_cluster")))
	filteredFiles = map(lambda x: "output_cluster" + "/" +  x , filteredFiles)

	for i in  filteredFiles:
		extactResultFiles(i)
	return filteredFiles

def extactResultFiles(path):
	f = open(path, "r")
	os.system ("mkdir -p %s" % "output_reduced")
	lines = f.readlines()
    
	for i, val in  enumerate(lines):
		if i % 2 == 0:
			print "["+val.strip()+"]"
			cmd = "cp %s %s" % (val.strip(), "output_reduced")
			os.system (cmd)
        f.close()

#----------------------------------------------------------
# Main
#----------------------------------------------------------
if __name__ == "__main__":
    main (sys.argv)