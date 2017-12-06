#!/usr/bin/Rscript
# Given an input k it creates partitions of k elements and
# calculates the medoids for each one.

nCPUS = 4
USAGE = "USAGE: clusters.R <input filename|input dir> <output filename>\n"
options (width=400)
#--------------------------------------------------------------
# Main function
#--------------------------------------------------------------
main <- function (args) {
	#args = c("in/hivp", "in/hivp-tbl.log")
	args = c ("in", "out", 10)
	if (length (args) < 3) {
		cat (USAGE)
		quit ()
	}
	inputDir  = args [1]
	outputDir = args [2]
	sizeBin   = strtoi (args [3])

	nBins = splitFilesToBins (inputDir, sizeBin, outputDir)
}

#--------------------------------------------------------------
# Split (symlinks) files to bins of "sizeBin"
#--------------------------------------------------------------
splitFilesToBins <- function (inputDir, sizeBin, outputDir) {
	files <<- normalizePath (paste (inputDir, "/", list.files (inputDir),sep=""))
	bins <- split(files, ceiling(seq_along(files)/sizeBin))
	nBins = length (bins)

	for (k in 1:nBins) {	
		binFiles = bins [[k]]

		binDir = sprintf ("bin-%.5d", k)
		system (sprintf ("mkdir %s/%s", outputDir, binDir))
		for (f in binFiles) {
			cat ("\n", ">>>>", f, "\n")
			cmm = sprintf ("ln -s %s %s/%s/%s", f, outputDir, binDir, basename(f))
			print (cmm)
			system (cmm)
		}
	}
	return (nBins)
}

#--------------------------------------------------------------
# Call main function
#--------------------------------------------------------------

args = commandArgs (TRUE)
main (args)
