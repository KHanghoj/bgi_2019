This is an R based exercise. Either implement your own gene genealogy
generator, or copy paste the code below in an R session.
```
library(ape)

##functions returns list with 1) newick 2) phylo 3) tks
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
```
Generate a tree with 3 nodes, and examine the contents

```
T4 <- coal(3)
names(T4)
```

It contains: 1) The tree represented in the 'newick' format. 2) The tree as a
phylo object. 3) The interarrival waiting times.

Plot four trees with 5 nodes

```
plot.phylo(coal(5)$phyl)
plot.phylo(coal(5)$phyl)
plot.phylo(coal(5)$phyl)
plot.phylo(coal(5)$phyl)
plot.phylo(coal(5)$phyl)
```

The expected interarrival times are given by 1/(k(k-1)/2). The
simulated interarrival times are given by

```
coal(5)$tks
1/choose(5:2,2)
```
How close are they? Lets try to do multiple replicates

```
rowMeans(replicate(10,coal(5)$tks))-1/choose(5:2,2)
rowMeans(replicate(100,coal(5)$tks))-1/choose(5:2,2)
rowMeans(replicate(10000,coal(5)$tks))-1/choose(5:2,2)
```