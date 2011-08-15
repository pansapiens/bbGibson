tms<-read.table("1000_pairs_30bp_70deg_filt.dat", sep="\t", header=TRUE)
hist(tms$Tm, xlab="Tm", ylab="number of pairs", main="All-against-all melting temperatures of 1000 primer pairs")
