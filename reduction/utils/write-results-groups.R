
#----------------------------------------------------------
# NOT USED
# Write the groups' representatives to the ouput dir
#----------------------------------------------------------
writeRepresentatives <- function (lstGroups, outputDir) {
	for (g in lstGroups) {
	  	pdb = g[[1]]
		cmm = sprintf ("ln -s %s %s/%s", pdb, outputDir, basename (pdb))
		#print (cmm)
		system (cmm)
	}
}
#----------------------------------------------------------
# NOT USED
# Write the groups formatted as lines of text
#----------------------------------------------------------
writeGroups <- function (lstGroups, outputDir=NULL) {
  cat (">>> Writing representatives...\n")
  outFile = paste (outputDir, "representatives.txt", sep="/")
  sink (outFile)
	for (g in lstGroups) {
	  cat (g[[1]])
	  cat ("\n")
	}
  sink ()

  cat (">>> Writing groups...\n")
  outFile = paste (outputDir, "groups.txt", sep="/")
  sink (outFile)
	for (g in lstGroups) {
	  cat (unlist (g))
	  cat ("\n")
	}
  sink ()
  cat (">>> End writing \n")
}




