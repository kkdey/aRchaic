

tab <- read.csv("../data/Pinhasi/IR1.q30.csv", header=FALSE)
tab1 <- tab[1:50000,]
write.csv(tab1, file = "../data/Pinhasi/IR1-subsampled.q30.csv", col.names=FALSE)


tab1 <- read.csv("../data/Pinhasi/IR1-subsampled.q30.csv", header=FALSE)
