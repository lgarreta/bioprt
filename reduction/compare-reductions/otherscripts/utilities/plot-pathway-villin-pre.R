#!/usr/bin/Rscript

# Plot properties for the villin preliminary study
# Make a XY plot of a pathway given as a parameter
# INPUT: The whole values, the pathway name, the average window 
# OUTPUT: A pdf file named as the protein pathway
# In the same figure, it is plotted the averaged ploted (n=51)

################################################################
# Log
# Nov 01/2013: Hidded the average line, legend, and name.
# Statict values file "vldm.values"
################################################################

usage <- function () {
	print ("Make a XY plot of a pathway given as a parameter")
	print ("Usage: program <values file> <name> <average window>")
	print ("INPUT: The values of the measurenments on the pathway conformations")
	print ("OUTPUT: A pdf file named as the protein pathway")
	print ("In the same figure, it is plotted the averaged ploted (n=51)")
	quit()
}

##############################################################################
# MAIN
# INPUT: (1) name of the pathway file (2) Average window
###############################################################################
mainR <- function () {
	args = commandArgs (TRUE)
	args = c("vldm.values", "villin", 100)

	## Take and assign command line args (or assign by default)
	if (length (args) == 0)
		usage ();

	valuesFilename = args [1]
	pathName       = args [2]
	averageWindow  = as.integer (args [3])

	#pathName = unlist (strsplit (valuesFilename, split="-"))[1]

	## load the values and create the average file
	values     = read.table (valuesFilename, sep="",header=T,row.names=1)
	#pathValues = getPathwayValues (values, pathName)
	pathValues = values

	## Plot for each property with its average
	for (property in colnames (pathValues)) {
		#outFilename = sprintf ("property-plot-avr%s-%s-%s.pdf", averageWindow, pathName, property)
		outFilename = sprintf ("%s-pre-pca-property-%s.pdf", pathName, property)  # When average is plotted

		propValues    = pathValues [,property]
		log (property, outFilename)
		propAvrValues  = averageFilter (propValues, averageWindow) 

		plotXYProperty (property, propValues, propAvrValues, outFilename)
	}
}

###############################################################################
# Creata a XY plot of a  property
###############################################################################
plotXYProperty <- function (property, propValues, avrPropValues,  outFilename) {
	
	pdf (outFilename, width=2.8, height=2)
		par (mar=c(2.8,2.8,1,1), cex=0.7, mgp = c(1.5,0.8,0)) # Outer margins, scale fonts, pos labels (

		n = length(avrPropValues) - 1
		time =(((0:n)*50)*(1/1e12)*(1e9/1))
		name = properties [property]
	 	plot (time, propValues, xlab="Time (ns)", ylab=name, type="l", lty="dashed", col=1)
		#legend (x="top", legend=c("Original","Averaged"), lty=c(2,1), col=c(1,2))

		par (new=T)
	 	#lines (avrPropValues, col="red", lty=1)
	dev.off()
}

###############################################################################
# Ids and long names to use in plots
###############################################################################
properties = c (	
	NC = "Native Contacts (#)",
	CO = "Contact Order (%)",
	RG = "Radius of Gyration (Angstrom)",
	HB = "Hydrogen Bonds (#)",
	AS = expression (paste ("Accessible Surface Area (", Angstrom^2, ")")),
	VD = expression (paste ("Voids (", Angstrom^3, ")")),
	RM = "RMSD (Angstrom)",
	PE = "Potential Energy (Debyes)",
	DM = expression (paste ("Dipolar Moment (", Debyes^2, ")"))
	)
###############################################################################
# Print a log message with the parameter
###############################################################################
log <- function (label, params=c()) {
	cat (">>>> ", label, ": ", params, "\n") 
}

###############################################################################
# Smooth the pathway values by se self implementation of the moving averaging 
# Remove noisy values -- Better smoothing at the ends that other filtering meths
###############################################################################
averageFilter <- function (data, n) {
	n = as.integer (n/2)
	size = length (data)
	filteredData = data
	for (i in 1:size) {
		if (i < n) {
			start = 1; end = 2*i-1
		}else if (i >  (size-n)) {
			start = i - (size-i); end=size
		}else {
			start = i-n; end=i+n
		}
		filteredData [i] <- mean (filteredData [start:end])
	}
	return (filteredData)
}
##############################################################################
# Average (smooth) the whole measures of a pathway calling the "averageFilter"
# Smooth the pathway values by self implementation of the moving averaging 
###############################################################################
averagePathway <- function (pathwayNamefile, win, prefixName) {
	values = read.csv (pathwayNamefile, sep="",header=T,row.names=1)

	for (prop in names (values)) {
		propValues = values [prop][,1]
		values [prop] = averageFilter (propValues, win) 
	}
	pathName = unlist (strsplit (pathwayNamefile,split="[.]")[[1]])[1]
	outputFilename = sprintf ("%s-%s.%s", pathName, prefixName, "values")
	write.table (values, file=outputFilename, sep=" ")
	return (outputFilename)
}

##############################################################################
# Get the pathway values of a specified pathway name
##############################################################################
getPathwayValues  <- function (values, name) {
	# get the full values given partial row name
	filter <- function (path) {values [grep (path, rownames (values)),]}

	dataset = filter (name)
	
	return (dataset)
}

##############################################################################
# call to main
##############################################################################
mainR ()
