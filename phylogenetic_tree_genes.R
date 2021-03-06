library(phangorn)
library(ape)
genes = c("NP", "L", "GP", "VP24", "VP30", "VP35", "VP40")
for(gene in genes){
  path=paste("~/bio_project/bio/output/",gene,".csv",sep="")
  matrix=read.csv(path,check.names=FALSE)
  tree=upgma(matrix)
  path=paste("~/bio_project/bio/output/upgma_",gene,".png",sep="")
  png(path)
  name=paste(gene,"UPGMA","Phylogenetic","Tree")
  plot(tree,main=name)
  dev.off()
  tree=nj(as.dist(matrix))
  path=paste("~/bio_project/bio/output/nj_",gene,".png",sep="")
  png(path)
  name=paste(gene,"NJ","Phylogenetic","Tree")
  plot(tree,main=name)
  dev.off()
}