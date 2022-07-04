# Concatenate files from each chr
```
module load StdEnv/2020  gcc/9.3.0 bcftools/1.13
```
```
bcftools concat DB_new_chr1L_genotyped.vcf DB_new_chr1S_genotyped.vcf DB_new_chr2L_genotyped.vcf DB_new_chr2S_genotyped.vcf DB_new_chr3L_genotyped.vcf DB_new_chr3S_genotyped.vcf DB_new_chr4L_genotyped.vcf DB_new_chr4S_genotyped.vcf DB_new_chr5L_genotyped.vcf DB_new_chr5S_genotyped.vcf DB_new_chr6L_genotyped.vcf DB_new_chr6S_genotyped.vcf DB_new_chr7L_genotyped.vcf DB_new_chr7S_genotyped.vcf DB_new_chr8L_genotyped.vcf DB_new_chr8S_genotyped.vcf DB_new_chr9_10L_genotyped.vcf DB_new_chr9_10S_genotyped.vcf -O z -o clivii_unfiltered_allchrs.vcf.gz
```
# Plink

Associations between SNPs and a phenotype (such as phenotypic sex) are easily calculated using plink.

```
module load nixpkgs/16.09 plink/1.9b_5.2-x86_64
module load StdEnv/2020
module load r/4.0.2
plink --vcf temp.vcf.gz --recode --const-fid 0 --chr-set 36 no-y no-xy no-mt --out myplink
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

