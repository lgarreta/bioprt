require ("clusteval")
set.seed(42)
K <- 3
n <- 10
labels1 <- sample.int(K, n, replace = TRUE)
labels2 <- sample.int(K, 8, replace = TRUE)
jc = jaccard (labels1, labels2)
print (jc)
