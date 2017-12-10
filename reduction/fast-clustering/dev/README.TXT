Programs that make trayectory reduction using a fast clustering.
The main program that call the other is the pr00_XXXX.py. 
It takes as parameter an input dir filename where the protein conformations are located. Then, it calls de "createBins" program (pr01) to split the trayectory files in blocks of 1000 conformations called bins (create links).  Next, the main program calls the reduction program (pr02) that processes each bin by firstly making a fast clustering then a full clustering.  

The fast clustering follows the hobbom algorithm by ordering the files, in our case they are really ordered by time, and taking the first as the representative for the first cluster and comparing the next ones with this representative, if they are closer (RMSD < 1.3) then they are assigned to this cluster, if not, this conformation is taken as the next representative for the next group, and so on.

The input of the program is a dir filename with the protein conformations (*.pdb) of the trajectory, and the output is an outpudir ("outrepr") with the medoids or representatives for each cluster (bin0000XXX) and the bins with the protein conformations resulting from the fast clustering.

NOTE:
In the program "pr02_reduction.R" the "partialClustering" makes the fast clustering, and the "fullClustering" makes the clustering using medoids (k-means).
The fullClustering is using a own methods (getRmsdDistanceMatrix) to create the RMSD distance matrix, but exists other "bio3d" methods to to the same (getRmsdDistanceMatrixWithAlignments). So it is needed to test if they are producing the same results, so try first with the first method, save the results, and then try with the second method, and compare the results.

The *.tgz files are the results with the "in" inputDir that contains a small set of
protein conformations. the other *.tgz are the results ("outbins", "outrepr"). Of course it is needed to test with a full trayectory.

In the main program you have to set the constant parameters for the reduction, mainly the block size (e.g. 10 for this example) and the RMSD threshold (e.g. 1.3), and the number of cores (e.g. 1). Change it to real values (e.g. 1000 for block size). You can try to use more cores (NCORES) (I have not tried this).
