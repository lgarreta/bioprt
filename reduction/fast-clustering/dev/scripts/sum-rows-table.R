m=read.table ("bin000005.dist", header=T)
for (i in 1:12) print (paste (i, sum (m[i,1:12])))
