#--------------------------------------------------------------
# Split (symlinks) files to bins of "sizeBin"
#--------------------------------------------------------------
splitFilesToBins <- function (inputDir, sizeBin, outputDir) {
	tmpDir = sprintf ("%s/%s", outputDir, "tmpBins")
	createDir (tmpDir)

	files <<- normalizePath (paste (inputDir, "/", list.files (inputDir),sep=""))
	bins <- split(files, ceiling(seq_along(files)/sizeBin))
	nBins = length (bins)

	binDirList = list ()
	for (k in 1:nBins) {	
		binFiles = bins [[k]]

		binDir = sprintf ("%s/bin-%.5d", tmpDir, k)
		cat ("\n>>> Bin dir ", binDir,"\n")
		binDirList = append (binDirList, binDir)
		createDir (binDir)

		for (f in binFiles) {
			system (sprintf ("ln -s %s %s/%s", f, binDir, basename(f)))
		}
	}
	return (binDirList)
}
