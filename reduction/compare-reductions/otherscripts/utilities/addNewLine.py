#!/usr/bin/python

import os, sys

args = sys.argv
inputDir = args [1]

filenameList = os.listdir (inputDir)

outputDir = "new-%s" % inputDir
os.system ("mkdir %s" %outputDir)
for filename in filenameList: 
    inFilename = "%s/%s" % (inputDir, filename)
    inFile      = open  (inFilename)

    outFilename = "%s/%s" % (outputDir, filename)
    outFile     = open (outFilename, "w")
    for line in inFile:
        outFile.write (line.strip())
        outFile.write ("\n")

        
