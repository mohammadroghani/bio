library(phangorn)
library(ape)
matrix=read.csv("~/bio_project/bio/output/global.csv",check.names=FALSE)
tree=upgma(matrix)
png("~/bio_project/bio/output/global_tree.png")
plot(tree,main="Global Tree")
dev.off()