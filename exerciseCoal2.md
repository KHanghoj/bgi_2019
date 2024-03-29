Copy paste the code below into an R terminal.

This code puts mutations on the simulated gene trees from last
exercise under the assumption of an infinite sites mutation model.

In addition to the mutations it also generates the sequences as a
vector of zeroes and ones.

And it also calculates the site frequency spectrum.

```
library(ape)

#functions returns list with 1) newick 2) phylo 3) tks
coal <- function(nnodes){
    nodes <- rep(0,nnodes)
    names(nodes) <- 1:nnodes
    tks <- c()
    T <- 0
    while(length(nodes)>1){
        pick <- sample(1:length(nodes),2)
        tk <- rexp(1,rate=choose(length(nodes),2))
        tks <- c(tks,tk)
        T <- T+tk
        pair <- nodes[pick]
        nodes <- c(nodes[-pick])
        inode <- paste("(", names(pair[1]), ":", T-pair[1], ",", names(pair[2]), ":", T-pair[2], ")", sep="")
        nodes <- c(nodes,T)
        names(nodes)[length(nodes)] <- inode
    }
    names(nodes) <- paste0(names(nodes[1]),";")
    phylo <- read.tree(text=names(nodes[1]))
    return(list(newick=names(nodes[1]),phyl=phylo,tks=tks))
}
##function calculates number of mutations on each branch
addmut <- function(mat,theta)
   sapply(mat,function(x) rpois(1,theta*x/2))


getSeq <- function(phyl){
    edges <- phyl$edge
    mut <- phyl$mut
    nmut <- sum(mut) ##num of segsites
    nsam <- length(phyl$tip)
    seqs <- matrix(0,ncol=nmut,nrow=nrow(edges))
    at <- 0
    for(i in 1:length(mut)){
        if(mut[i]>0)
            seqs[i,(at+1):(at+mut[i])] <- 1
        at <- at+mut[i]
    }
    ##seqs now contains the sequence mutation for each edge
    chain <- function(x){
        persam <- rep(FALSE,nmut)
        while(TRUE){
##            cat('looking for: ',x,'\n')
            wh <- which(edges[,2]==x)
            if(length(wh)==0)
                break
            persam <- persam+seqs[wh,]
            x <- edges[wh,1]
        }
        return(persam)
    }
    return(t(sapply(1:length(phyl$tip.label),chain)))
    
}
getSFS <- function(seqs){
    return(table(factor(colSums(seqs),levels=0:nrow(seqs))))
}

wrapper <- function(nnodes,theta){
    tmp <- coal(nnodes)
    tmp$phyl$"mut" <- addmut(tmp$phyl$edge.length,theta)
    tmp$"seqs" <- getSeq(tmp$phyl)
    tmp$"sfs" <- getSFS(tmp$seqs)
    return(tmp)
}

```

Try to visualize the sfs for a single tree
```
T20 <- wrapper(20,10)
barplot(T20$sfs)
```
How does it look? The SFS is normally estimated genomwide how many
sites do we have

```
ncol(T20$seqs)
```

Lets try to make more realizations:

```
barplot(rowMeans(replicate(10,wrapper(20,10)$sfs)))
barplot(rowMeans(replicate(100,wrapper(20,10)$sfs)))
barplot(rowMeans(replicate(1000,wrapper(20,10)$sfs)))
```

Let us also try to estimate theta using both watterson and tajima

```
sfs <- rowMeans(replicate(1000,wrapper(20,10)$sfs))
##watterson
sum(sfs)/sum(1/(1:19))
##pi
sum(0:20*20:0*sfs)/choose(20,2)
```
