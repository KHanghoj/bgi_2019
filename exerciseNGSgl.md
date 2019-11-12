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

TODO: PATHS MUST BE FIXED

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

#### Examples of I/O filters ####

Below are some examples of I/O filters for ANGSD that we will use today:

```
# -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1
```

These filters will retain only uniquely mapping reads, not tagged as bad, considering only proper pairs, without trimming, and adjusting for indel/mapping (as in samtools).
`-C 50` reduces the effect of reads with excessive mismatches, while `-baq 1` computes base alignment quality as explained here ([BAQ](http://samtools.sourceforge.net/mpileup.shtml)) used to rule out false SNPs close to INDELS.

We can also add filters related to the quality of the NGS data:
```
# -minMapQ 20 -minQ 20
```
This will remove reads with mapping quality below 20, basequalities < 20 and

Criteria to exclude sites with low/high amounts of data can also be set:
```
# -minInd 5 -doCounts
# -setMinDepth 7 -setMaxDepth 30
```
Sites where less that 5 individuals have reads mapping.
This last filter ensures that we only analyse sites with a minimum depth of 7 and maxmimum depth of 30.

### Genotype Likelihoods ###
We now wish to calculate the ***genotype likelihoods*** for each site at each individual.

![stages1](https://github.com/mfumagalli/Copenhagen/blob/master/Files/stage1.png)

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

** QUESTIONS **

1. What are the output files?
2. What kind of information do they contain?
3. How many columns are present in the `glf.gz` file? Is that excepted?


```
ls Results/EUR.*
```

```
less -S Results/EUR.arg
less -S Results/EUR.glf.gz
```


### Genotype probability and Genotype Calling ###

Here, we will explore several ways to call genotypes from sequencing data.
We will also calculate genotypes probabilities (aka posterior probability of the genotypes) for each site for each individual.

![stages2](https://github.com/mfumagalli/Copenhagen/blob/master/Files/stage2.png)

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
1. How would the output files differ between `-doGeno 2` and `-doGeno 3`?
2. What kind of prior have we used?
3. How many sites in total and how many have at least one missing genotype (-1)?


```
zcat Results/EUR.geno.gz | wc -l
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
6. What has happened?

### EXERCISE 1  ###

In this exercise we will use *ANGSD* to calculate the allele frequency of a variant in the EDAR gene in 5 human populations.

The SNP of interest (rs3827760) is described in [this paper](http://www.nature.com/ncomms/2016/160519/ncomms11616/full/ncomms11616.html). Briefly,

TODO


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

#### Build the ANGSD command ####

Write the code that performs the following genotype calling for EDAR variants in all populations.
You can directly call genotypes by adding the neccesary arguments. This way we do not need to first generate `.glf.gz` file.


As an indication, you can follow these guidelines:
- use the SAMtools genotype likelihood model
- calculate genotype posterior probabilities using a HWE-based prior (allele frequencies)
- Add the base filters shown in the I/O section (without min/max depth)
- filter out bases with a quality score less than 20
- filter our reads with a mapping quality score less than 20
- use ony sites where you have at least one sample with data (-mindInd)
- use -doMaf 1 (major and minor is fixed)
- Skip triallelic sites
- set genotypes as missing if the highest genotype probability is less than 0.50
- Make sure to set major allele to ancestral ($ANC) hint: -doMajorMinor


```bash
for pop in AFR EUR EAS LAT NAM
do
    $NGS/angsd/angsd -b $DATA/${pop}.bams -out Results/EDAR.${pop} -sites Data/snp.txt ...
done
```

**QUESTIONS**
1. What is the allele frequency for each populations based on the genotype calls?
2. Plot the allele frequencies for all the populations for this snp
3. Do we observe between population differences?


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

`knownEM` column specifies the algorithm used to estimate the allele frequency (of the minor allele) which is given under that column.

It is evident that these sites are monomorphic. In order to to perform **SNP calling**, we can use two different approaches:

```
        -minMaf         0.000000        (Remove sites with MAF below)
        -SNP_pval       1.000000        (Remove sites with a pvalue larger)
```

`-minMaf` is a hard filter of the estimated allele frequency

`-SNP_pval` is drawn from a chisq distribution based on a likelihood ratio test.

### QUESTIONS ###

1. How many variants with an minor allele frequency `>=5%`
2. How many variants with an p-value `<=1e-6` hint: `snp pval 1e-6`


### EXERCISE II ###
Now we have explored how to obtain allele frequencies without calling genotypes and how to obtain genotype probabilities. Lets recalculate allele frequencies of the EDAR variant for each population using the likelihood approach.

1. What is the allele frequency for each populations?
2. Are these estimates any different from those obtained from the genotype calling procedure?
