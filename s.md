```
for POP in AFR EUR EAS LAT NAM
do
        echo $POP
        $ANGSD/angsd -b $DATA/$POP.bams -ref $REF -anc $ANC -out Results/EDAR.$POP \
                -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
                -minMapQ 20 -minQ 20 -minInd 1 \
                -GL 1 -doMajorMinor 5 -doMaf 1 -skipTriallelic 1 \
                -doGeno 3 -doPost 1 -postCutoff 0.50 \
                -sites Data/snp.txt
done
```
