```
for POP in AFR EUR EAS LAT NAM
do
        echo $POP
        $NGS/angsd/angsd -b $DATA/$POP.bams -ref $REF -anc $ANC -out Results/EDAR.$POP \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
                -minMapQ 20 -minQ 20 -minInd 1 \
                -GL 1 -doMajorMinor 5 -doMaf 1 -skipTriallelic 1 \
                -doGeno 3 -doPost 1 -postCutoff 0.50 \
                -sites Data/snp.txt
done
```


```
for POP in AFR EUR EAS LAT NAM; do echo -n "$POP " ;zcat Results/EDAR.${POP}.geno.gz | cut -f 5- | Rscript -e 'message(mean(scan("stdin", quiet=TRUE))/2)'; done
```


```
for POP in AFR EUR EAS LAT NAM; do echo -n "$POP " ;zcat Results/EDAR.${POP}.mafs.gz | sed '1d' | awk '{print $(NF-1)}'; done
```
