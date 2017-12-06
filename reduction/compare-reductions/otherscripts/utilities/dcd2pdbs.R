#!/usr/bin/Rscript

library (bio3d)

args = commandArgs (TRUE)
#pathname  = args [1]
pathname  = "hivp.dcd"
outpath   = "hivp"

trj = read.dcd (pathname)
n = nrow (trj)

for (i in 1:n) {
	name = sprintf ("hivp/pdb%03d.pdb", i)
	write.pdb (file=name, xyz=trj[i,])
}

