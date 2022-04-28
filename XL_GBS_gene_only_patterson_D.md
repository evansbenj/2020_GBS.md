# Make a bed file with only coordinates of genes on known L or S chrs:
locally:
```
cat XLv9.2_xenbase_annotations.gff | grep 'gene' | cut -f1,4,5 > gene_only_bed.bed
```
# Now use excel to add a buffer of 10,000 bp on either side
