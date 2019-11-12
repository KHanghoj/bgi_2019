# Genotype likelihoods and SNP calling

** ALways give links for the files **

The exercises are a modified version of this [[INSERT LINK TO FUMAGALLI]]

In this session you will learn how to do:
* genotype calling
* allele frequency estimation
* variant (or SNP) calling


We are using the program [ANGSD](http://popgen.dk/wiki/index.php/ANGSD) (Analysis of Next Generation Sequencing Data).
More information about its rationale and implemented methods can be found [here](http://www.ncbi.nlm.nih.gov/pubmed/25420514).


## First set some paths

Briefly, you need to set the path to the software and various data that will be used.
Also, you will have to create two folders on your working directory, one for your results and one for your intermediate data.

PATHS MUST BE FIXED

```
NGS=/ricco/data/matteo/Software/ngsTools
DIR=/home/matteo/Copenhagen
DATA=/ricco/data/matteo/Data
REF=$DATA/ref.fa.gz
ANC=$DATA/anc.fa.gz

mkdir Results
mkdir Data
```


### Workflow

![stages](https://github.com/mfumagalli/Copenhagen/blob/master/Files/stages.png)

This might seem a bit overwhelming at first, but we will go through every step in this exercise.

The workflow can be divided into four steps

1. Filtering data
2. Genotype likelihoods
3. Genotype probability and Genotype Calling
4. SNP calling

### Data Filtering

In contrast to Exercise I, we will be using multiple bam (individuals) files.

Have a look at our list of BAM files:
```
cat $DATA/ALL.bams
wc -l $DATA/ALL.bams
ls $DATA/*.bams
```

The default data filtering done by ANGSD can be found using the following command

```
$NGS/angsd/angsd -bam
```

1. What is the default base quality score filter


Below we have shown the base of the ANGSD command:

```
# $NGS/angsd/angsd -b ALL.bams -ref $REF -out Results/ALL \
#        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
```

These filters will retain only uniquely mapping reads, not tagged as bad, considering only proper pairs, without trimming, and adjusting for indel/mapping (as in samtools).
`-C 50` reduces the effect of reads with excessive mismatches, while `-baq 1` computes base alignment quality as explained here ([BAQ](http://samtools.sourceforge.net/mpileup.shtml)) used to rule out false SNPs close to INDELS.

We will use an additional set of filters throughout the exercise.
```
#        -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 7 -setMaxDepth 30 -doCounts 1
```
This will remove reads with mapping quality below 20, basequalities < 20 and sites where less that 5 individuals have reads mapping.

`-setMinDepth 7 -setMaxDepth 30` are per bam filters. They ensure that we only use sites with a depth between 7 and 30.

### Genotype Likelihoods ###
We now wish to calculate the ***genotype likelihoods*** for each site at each individual.

To do so you need to specify which genotype likelihood model to use.
```
$NGS/angsd/angsd -GL
```

For this exercise, we will use `-GL 2` (GATK) that was introduced in the lecture.

Run the following command:
A possible command line to estimate genotype likelihoods might be:
```
$NGS/angsd/angsd -b $DATA/EUR.bams -ref $REF -out Results/EUR \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 7 -setMaxDepth 30 -doCounts 1 \
        -GL 2 -doGlf 4
```
where we specify:
* -GL 2: genotype likelihood model as in GATK
* -doGlf 4: output in text format

Ignore the various warning messages

** QUESTINS **

1. What are the output files?
2. What kind of information do they contain?


```
ls Results/EUR.*
```

```
less -S Results/EUR.arg
less -S Results/EUR.glf.gz
```


### Genotype probability and Genotype Calling ###

Here we will explore several ways to call genotypes from sequencing data.
We will also calculate genotypes probabilities (aka posterior probability of the genotypes) to each site for each individual.

In ANGSD, the option to call genotypes is `-doGeno`:
```
$NGS/angsd/angsd -doGeno
```

```
$NGS/angsd/angsd -b $DATA/EUR.bams -ref $REF -out Results/EUR \
    -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 7 -setMaxDepth 30 -doCounts 1 \
    -GL 2 -doGlf 1

$NGS/angsd/angsd -glf Results/EUR.glf.gz -fai $REF.fai -nInd 10 -out Results/EUR \
    -doMajorMinor 1 -doGeno 3 -doPost 2 -doMaf 1
```

**QUESTION**
1. What files will be produced when with `-doGeno 3`?
2. How many sites have at least one missing genotype (-1)?

```
zcat Results/EUR.geno.gz | grep -1 - | wc -l
```
3. Why is that?

You can control how to set missing genotype when their confidence is low with `-postCutoff`.
For instance, we can set as missing genotypes when their (highest) genotype posterior probability is below 0.95:
```
$NGS/angsd/angsd -glf Results/EUR.glf.gz -fai $REF.fai -nInd 10 -out Results/EUR \
        -doMajorMinor 1 -doGeno 3 -doPost 2 -doMaf 1 -postCutoff 0.95
```

4. How many sites do we have in total?
5. How many sites have at least one missing genotype now?



### EXERCISE 1  ###

In this exercise we will use *ANGSD* to calculate the allele frequency of a variant in the EDAR gene in 5 human populations.

The SNP of interest (rs3827760) is described in [this paper](http://www.nature.com/ncomms/2016/160519/ncomms11616/full/ncomms11616.html).


The bam lists for each population can be found here:

```
ls $DATA/*.bams
```

*Do not use `ALL.bams`*

We retain only these populations: AFR (Africans), EUR (Europeans), EAS (East Asians), LAT (Latinos), NAM (Native Americans).


#### Restricting analysis to a genomic region ####

For this exercise, we are only interested in a single genomic position.

In ANGSD we can restrict our analyses on a subset of positions of interest using the `-sites` option.
The file with these positions need to be formatted as (chromosome positions).
```
echo 2 109513601 > Data/snp.txt
```
We need to index this file in order for ANGSD to process it.
```
$NGS/angsd/angsd sites index Data/snp.txt
```

Write the code that performs the following genotype calling for EDAR variants in all populations.
Also, you can directly call genotypes without generating the genotype likelihood files, by starting from bam files directly.
As an indication, you can follow these guidelines:
- use the SAMtools genotype likelihood model
- calculate genotype posterior probabilities using a HWE-based prior
- filter out bases with a quality score less than 20
- filter our reads with a mapping quality score less than 20
- use ony sites where you have at least one sample with data (-mindInd)
- do not set any filtering based on min and max depth
- use -doMajorMinor 1 and -doMaf 1 options
- set genotypes as missing if the highest genotype probability is less than 0.50
but feel free to choose some parameters yourself.



```bash
for pop in AFR EUR EAS LAT NAM
do
    $NGS/angsd/angsd -b $DATA/${pop}.bams -out Results/${pop} -sites Data/snp.txt
done
```

**QUESTIONS**
1. What is the allele frequency for each populations?
2. Do you observe a difference between the populations?



### SNP calling ###
In the last section we will calculate the allele frequency without using called genotypes. Instead we will estimate allele frequencies taking into account the uncertainty in the data by using genotype likelihoods.

Have a look at `-doMaf`

```
$NGS/angsd/angsd -doMaf
```

To use -doMaf, we need to also assign major and minor alleles:

```
$NGS/angsd/angsd -doMajorMinor
```


A possible command line to estimate allele frequencies might be (this may take 1 min to run):
```
$NGS/angsd/angsd -b $DATA/EUR.bams -ref $REF -out Results/EUR \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -minQ 20 -minInd 5 -setMinDepth 7 -setMaxDepth 30 -doCounts 1 \
        -GL 1 -doGlf 1 \
    -doMajorMinor 1 -doMaf 1
```

Have a look at this file which contains estimates of allele frequency values.
```
zcat Results/EUR.mafs.gz | head
```

`knownEM` specifies the algorithm used to estimate the allele frequency (of the minor allele) which is given under that column.

3.
