# Diversity of sex-associated sites in boumbaensis

Unfortunately we do not have the parents of the boumbaensis family. This means we can't focus on maternal and paternal variation. Instead we can quantify diversity of sex-associated positions in daughters and sons. The sex with the higher diversity should be the heterogametic sex.

In this directory:
```
/home/ben/projects/rrg-ben/ben/2022_GBS_lotsofxennies/individual_gvcfs_by_species/2022_boum/mapped_to_XLv10_concatscaf
```
I did this using vcftools:
For ten sons:
```
zcat combined_Chr8L.g.vcf.gz_genotyped.vcf.gz | vcftools --vcf - --site-pi --positions SNP_list.txt --indv boum_M_M10_GCATACGGAG_sorted.bam --indv boum_M_M1_GTACCATGG_sorted.bam --indv boum_M_M2_CATCACGT_sorted.bam --indv boum_M_M3_CGTGGACAAT_sorted.bam --indv boum_M_M4_TCTTGG_sorted.bam --indv boum_M_M5_CTCAGCGTAG_sorted.bam --indv boum_M_M6_GACACT_sorted.bam --indv boum_M_M7_ATATGCGGT_sorted.bam --indv boum_M_M8_TACCGTCAT_sorted.bam --indv boum_M_M9_CTGTACCT_sorted.bam --out nucleotide_diversity_males
```
and five daughters:
```
zcat combined_Chr8L.g.vcf.gz_genotyped.vcf.gz | vcftools --vcf - --site-pi --positions SNP_list.txt --indv boum_F_F1_AGTTAGG_sorted.bam --indv boum_F_F2_CCATTCAGT_sorted.bam --indv boum_F_F3_TCAGTACG_sorted.bam --indv boum_F_F4_GGTTAGCT_sorted.bam --indv boum_F_F5_GTGCAT_sorted.bam --out nucleotide_diversity_females
```
