#!/usr/bin/Rscript

require(bio3d); require(graphics);

log <- function (object, filename) {
	sink (filename)
	print (object)
	sink()
}
pause <- function() {
  cat("Press ENTER/RETURN/NEWLINE to continue.")
  readLines(n=1)
  invisible()
}

# Read example trajectory file
trtfile <- system.file("examples/hivp.dcd", package="bio3d")
trj <- read.dcd(trtfile)

# Read the starting PDB file to determine atom correspondence
pdbfile <- system.file("examples/hivp.pdb", package="bio3d")
#pdb <- read.pdb(pdbfile)
pdb <- read.pdb("pdb001.pdb")

ncol(trj) == length(pdb$xyz)

# Trajectory Frame Superposition on Calpha atoms
ca.inds <- atom.select(pdb, elety = "CA")
xyz <- fit.xyz(fixed = pdb$xyz, mobile = trj, 
	           fixed.inds = ca.inds$xyz, 
	           mobile.inds = ca.inds$xyz)

log (xyz[1,], "xyz-hivp.log")

# Root Mean Square Deviation (RMSD)
rd <- rmsd(xyz[1, ca.inds$xyz], xyz[, ca.inds$xyz])
plot(rd, typ = "l", ylab = "RMSD", xlab = "Frame No.")

points(lowess(rd), typ = "l", col = "red", lty = 2, lwd = 2)

summary(rd)

# Root Mean Squared Fluctuations (RMSF)
rf <- rmsf(xyz[, ca.inds$xyz])
plot(rf, ylab = "RMSF", xlab = "Residue Position", typ="l")
pause()

# Principal Component Analysis
pc <- pca.xyz(xyz[, ca.inds$xyz])
plot(pc, col = bwr.colors(nrow(xyz)))

# Cluster in PC space
#hc <- hclust(dist(pc$z[, 1:2]))
hc <- hclust(dist(xyz))

grps <- cutree(hc, k = 4)
log (grps, "groups-hivp.log")
plot(pc, col = grps)
cat ("Writing groups...\n")
write.table (file="hivp-grps.tbl", grps)

pause()

# Cross-Correlation Analysis
cij <- dccm(xyz[, ca.inds$xyz])
plot(cij)

## view.dccm(cij, pdb, launch = TRUE)
