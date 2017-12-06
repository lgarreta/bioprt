#!/usr/bin/python

"""
 Get from a multipdb (separated by "END") the range of Pdbs,
 complete aminos, and make with the resultinf files a .tgz file.
"""
"""
LOG:
	2013
	Nov/27: Modified to include the native as "...-r.pdb"
"""


import sys, os
import shutil
import gzip

args = sys.argv

if len (args) < 3:
	print "USAGE: %s %s" % (args[0], "<filename> <label>")
	sys.exit(0)

filename = args [1]
label = args [2]
outName = "%s" % (label)

n = 0
file = gzip.open (filename)

os.mkdir (outName)
os.chdir (outName)

buffer = []
for line in file:
	if not "END" in line:
		buffer.append (line)
	else:
		pdbName = "%s-%06d.pdb.gz" % (label, n)
		newFile = gzip.open (pdbName, "w")
		newFile.writelines (buffer)
		newFile.close()

		buffer = []
			
		if n % 100 == 0:
			print ">>>", pdbName

		n = n+1
