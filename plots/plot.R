df<-read.table("af_gt.txt", row.names=1)
bitmap("af_gt.png", res=300)
barplot(t(df), main = "AF based on called genotypes")
dev.off()

df<- read.table("af_gl.txt", row.names=1)
bitmap("af_gl.png", res=300)
barplot(t(df), main = "AF based on GL genotypes")
dev.off()
