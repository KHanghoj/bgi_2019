Purpose of this exercise is to simulate genetic drift using the Wright-Fisher model without mutation.

This exercise is in R. Either write your own simulator or use the code below:


```
wf_freq <- function(f,n=200,g=100,new=FALSE,...){
    if(f<1)
        f <- f*n
    if(new==T)
        plot(-1,xlim=c(1,g),ylim=c(0,1),xlab='Generations',ylab='Allele Frequency',...)
    x <- f
    for(i in 1:g){
        x <- sample(0:n,1,prob=dbinom(0:n,n,x/n))
        points(i,x/n,...)
    }
}
wf_freq(f=0.5,n=200,g=2000,new=T,col=1)
```

Try 6 different realizations with a population of 200  and overlay the results:

```
wf_freq(f=0.5,n=200,g=2000,new=T,col=1)
wf_freq(f=0.5,n=200,g=2000,new=T,col=2)
wf_freq(f=0.5,n=200,g=2000,new=T,col=3)
wf_freq(f=0.5,n=200,g=2000,new=T,col=4)
wf_freq(f=0.5,n=200,g=2000,new=T,col=5)
wf_freq(f=0.5,n=200,g=2000,new=T,col=6)
```

Try 6 different realizations with a population of 1000  and overlay the results:

```
wf_freq(f=0.5,n=200,g=2000,new=T,col=1)
wf_freq(f=0.5,n=200,g=2000,new=T,col=2)
wf_freq(f=0.5,n=200,g=2000,new=T,col=3)
wf_freq(f=0.5,n=200,g=2000,new=T,col=4)
wf_freq(f=0.5,n=200,g=2000,new=T,col=5)
wf_freq(f=0.5,n=200,g=2000,new=T,col=6)
```
