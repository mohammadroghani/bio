genome_length=c(18890,18874,18934,18954,18939)
global=read.csv(file ="~/bio_project/bio/output/global.csv",check.names=FALSE)
times=global
for (i in 1:5)
  for(j in 1:5)
    times[[i,j]]=-1/(1.9*10^(-3))*log(1-4/3*global[[i,j]]/((genome_length[i]+genome_length[j])/2))
times