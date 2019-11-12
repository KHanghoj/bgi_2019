df<-read.table("af_gt.txt", row.names=1)
pdf("af_gt.pdf")
barplot(t(df), main = "AF based on called genotypes")
dev.off()

df<- read.table("af_gl.txt", row.names=1)
pdf("af_gl.pdf")
barplot(t(df), main = "AF based on Genotype Likelihoods")
dev.off()
