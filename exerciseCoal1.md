This is an R based exercise. Either implement your own gene genealogy
generator, or copy paste the code below in an R session. 

For this exercise you new to use the R packages `ape`. If it is not installed on your laptop you can run the following code:

```R
install.packages("ape")
```


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
Finally,	imagine	the	following	case:	a	researcher	has
estimated	the	structure	of	a	tree	for mtDNA	from a
species	sampled	in	a	single	location.	She	obtains	a	tree
looking	as	follows:

![alt text](https://github.com/KHanghoj/bgi_2019/blob/master/Tree0.png)

Based	on	the	structure	of	the	tree,	i.e.	two	groups	of	related	individuals	separated	by	long branches	down to the	root	of	the	tree,	she	concludes	that	there	must	be population subdivision	with	two	clearly differentiated	groups.	Based	on	what	you	have	learned	from the	simulations,	do	you	agree	with	this conclusion?
