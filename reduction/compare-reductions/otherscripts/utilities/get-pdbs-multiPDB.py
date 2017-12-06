#!/usr/bin/python

"""
 Get from a multipdb (separated by "END") the range of Pdbs,
 complete aminos, and make with the resultinf files a .tgz file. i
 At the end put the native 2FK4-nle.pdb
"""
"""
LOG:
	2013
	Nov/27: Modified to include the native as "...-r.pdb"
"""


import sys, os
import shutil

args = sys.argv

if len (args) < 4:
	print "USAGE: %s %s" % (args[0], "<filename> <label> <start> <end>")
	sys.exit(0)

filename = args [1]
label = args [2]
nStart = int (args [3])
nEnd   = int (args [4])

outName = "%s-%s-%s" % (label, nStart, nEnd)

n = 0
file = open (filename)

os.mkdir (outName)
os.chdir (outName)

buffer = []
for line in file:
	if not "END" in line:
		buffer.append (line)
	else:
		if n < nStart:
			n = n+1
			buffer = []
			continue
		if n >= nEnd:
			break

		pdbName = "%s-%06d.pdb" % (label, n)
		pdbName = "%s-%06d.pdb" % (label, n)
		newFile = open (pdbName, "w")
		newFile.writelines (buffer)
		newFile.close()

		buffer = []
			
		if n % 100 == 0:
			print ">>>", pdbName

		n = n+1

lastPdbName = "%s-%06d-r.pdb" % (label, 999999)
shutil.copy ("../2F4K-nonle.pdb", lastPdbName)

