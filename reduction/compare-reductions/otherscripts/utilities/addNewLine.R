
args = commandArgs (T)
inputDir = args [1]

filenameList = list.files (inputDir)

newDir = sprintf ("new-%s", inputDir)
system (sprintf ("mkdir %s", newDir))
for (f in filenameList) {
	file    = read.csv (sprintf ("%s/%s", inputDir, f))
	newFilename = sprintf ("%s/%s", newDir,f)

	write.csv (newFilename, file)
}

