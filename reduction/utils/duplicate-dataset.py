#!/usr/bin/python

import os, sys

# Given an input dir and an increment 
# it takes the files of the input dir and adds the increment to 
# the number of the file and makes a new copy of it
dir   = sys.argv [1]
inc   = int (sys.argv [2])
files = os.listdir (dir)
files.sort()

def sp (filename, increment): 
    prefix      = filename.split("-")[0]
    number      = filename.split("-")[1].split(".")[0]; 
    newNumber   = int (number)+increment;
    newFilename = prefix + "-" + str(newNumber).zfill(6)+".pdb"
    return newFilename

for i in files: 
    cmm="cp %s/%s %s/%s" % (dir, i, dir, sp (i, inc));
    print cmm; 
