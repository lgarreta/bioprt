#!/usr/bin/Rscript

args = commandArgs (TRUE)
inDir = args [1]
inDir = "in/tables/full"

fullFiles = paste (inDir, list.files (inDir), sep="/")
print (fullFiles)

fullTable = data.frame ()
for (f in fullFiles) {
	tbl = read.table (f)
	fullTable = rbind (fullTable,tbl)
}
write.table (fullTable, file="clustering-structures.tbl")


