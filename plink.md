# Concatenate files from each chr
```
module load StdEnv/2023  gcc/12.3 bcftools/1.19
```
```
bcftools concat DB_new_Chr1L_genotyped.vcf DB_new_Chr1S_genotyped.vcf DB_new_Chr2L_genotyped.vcf DB_new_Chr2S_genotyped.vcf DB_new_Chr3L_genotyped.vcf DB_new_Chr3S_genotyped.vcf DB_new_Chr4L_genotyped.vcf DB_new_Chr4S_genotyped.vcf DB_new_Chr5L_genotyped.vcf DB_new_Chr5S_genotyped.vcf DB_new_Chr6L_genotyped.vcf DB_new_Chr6S_genotyped.vcf DB_new_Chr7L_genotyped.vcf DB_new_Chr7S_genotyped.vcf DB_new_Chr8L_genotyped.vcf DB_new_Chr8S_genotyped.vcf DB_new_Chr9_10L_genotyped.vcf DB_new_Chr9_10S_genotyped.vcf -O z -o clivii_unfiltered_allChrs.vcf.gz
```
# Plink

Associations between SNPs and a phenotype (such as phenotypic sex) are easily calculated using plink.

```
module load nixpkgs/16.09 plink/1.9b_5.2-x86_64
module load StdEnv/2020
module load r/4.0.2
plink --vcf temp.vcf.gz --recode --const-fid 0 --chr-set 37 no-y no-xy no-mt --allow-extra-chr --out myplink
plink --file myplink --pheno sex_phenotype --assoc --allow-no-sex
```
where the "sex_phenotype" file is a tab-delimited file that looks like this:
```
0	sample1	1
0	sample2	2
0	sample3	1
0	sample4	2
0	sample5	2
0	sample6	1
```
The first column is the family ID (just zeros here).  The second column is the sample name - this is the same as in the vcf file.  The third column is the phenotype - should use 1 and 2, NOT 1 and 0, because 0 *might* be interpreted as a missing phenotye.

Also this flag "--const-fid 0" sets the family id to zero and tells plink to use the vcf sample name as the sample ID irrespective of whether there is an underscore in the name.

This flag "--chr-set 36" allows extra chrs.  They will be numbers in the order they are encountered in the vcf file (I think - this will need to be confirmed...)
module load nixpkgs/16.09 plink/1.9b_5.2-x86_64

# Example commandlines
```
module load nixpkgs/16.09 plink/1.9b_5.2-x86_64

plink --vcf pygmaeus_filtered_removed_allchrs.vcf.gz --recode --const-fid 0 --chr-set 36 no-y no-xy no-mt --allow-extra-chr --out pygmaeus_filtered_removed_allchrs.vcf.gz_myplink

plink --file pygmaeus_filtered_removed_allchrs.vcf.gz_myplink --pheno sex_phenotype --assoc --allow-no-sex --allow-extra-chr

mv plink.assoc pygmaeus_filtered_removed_allchrs.vcf.gz_plink.assoc
```

# get rid of rows with no P value
```R
setwd("./")
library(ggplot2)
dat<-read.table("./pygmaeus_filtered_removed_allchrs.vcf.gz_plink.assoc",header=TRUE)
newdat <- dat[!is.na(dat$P),]
write.table(newdat, file = "pygmaeus_filtered_removed_allchrs.vcf.gz_plink_noNAs.assoc", sep = "\t", quote = FALSE)
,row.names = TRUE, col.names = NA)
```
or using unix (much faster):
```
awk '$9 == "NA" { next } { print }' XB_WGS_not_filted_allchrs.vcf.gz_plink.assoc > XB_WGS_not_filtered_allchrs.vcf.gz_plink_noNAs.assoc
```
And further thin by getting rid of non-significant positions:
```
awk '$9 > 0.01 { next } { print }' allsites_Scafs.g.vcf.gz_genotyped.vcf.gz_myplink_plink_noNAs.assoc > allsites_Scafs.g.vcf.gz_genotyped.vcf.gz_myplink_plink_noNAs_Pval_gt_0.01.assoc
```
If the file is still large (it was for the XB WGS data) you can split the file by filing the first occurance of Chr5L like this:

```
grep -m1 'Chr5L' file_name
```

And then using head or sed to split it up (need to add the column header to the second half)...

to print a range of lines (eg lines 5-8 in this example):
```
sed -n 5,8p file
```
