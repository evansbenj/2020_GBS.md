# Diversity of sex-associated sites in boumbaensis

Unfortunately we do not have the parents of the boumbaensis family. This means we can't focus on maternal and paternal variation. Instead we can quantify diversity of sex-associated positions in daughters and sons. The sex with the higher diversity should be the heterogametic sex.

In this directory:
```
/home/ben/projects/rrg-ben/ben/2022_GBS_lotsofxennies/individual_gvcfs_by_species/2022_boum/mapped_to_XLv10_concatscaf
```
load vcftools:
```
module load StdEnv/2020 vcftools/0.1.16
```
I did this using vcftools:
For ten sons:
```
zcat combined_Chr8L.g.vcf.gz_genotyped.vcf.gz | vcftools --vcf - --site-pi --positions SNP_list.txt --indv boum_M_M10_GCATACGGAG_sorted.bam --indv boum_M_M1_GTACCATGG_sorted.bam --indv boum_M_M2_CATCACGT_sorted.bam --indv boum_M_M3_CGTGGACAAT_sorted.bam --indv boum_M_M4_TCTTGG_sorted.bam --indv boum_M_M5_CTCAGCGTAG_sorted.bam --indv boum_M_M6_GACACT_sorted.bam --indv boum_M_M7_ATATGCGGT_sorted.bam --indv boum_M_M8_TACCGTCAT_sorted.bam --indv boum_M_M9_CTGTACCT_sorted.bam --out boum_nucleotide_diversity_males
```
and five daughters:
```
zcat combined_Chr8L.g.vcf.gz_genotyped.vcf.gz | vcftools --vcf - --site-pi --positions SNP_list.txt --indv boum_F_F1_AGTTAGG_sorted.bam --indv boum_F_F2_CCATTCAGT_sorted.bam --indv boum_F_F3_TCAGTACG_sorted.bam --indv boum_F_F4_GGTTAGCT_sorted.bam --indv boum_F_F5_GTGCAT_sorted.bam --out boum_nucleotide_diversity_females
```


# Pygmaeus
directory:
```
/home/ben/projects/rrg-ben/ben/2022_GBS_lotsofxennies/individual_gvcfs_by_species/2022_pygmaeus/mapped_to_XLv10_concatscaf_unfiltered
```
```
vcftools --vcf DB__Chr8L_out.vcf --site-pi --positions SNP_list.txt --indv pyg_fem_Z23338_TGCA_sorted.bam --indv pyg_fem_Z23340_CAGA_sorted.bam --indv pyg_fem_Z23341_AACT_sorted.bam --indv pyg_fem_Z23342_GCGT_sorted.bam --indv pyg_fem_Z23343_CGAT_sorted.bam --indv pyg_fem_Z23344_GTAA_sorted.bam --indv pyg_fem_Z23345_AGCG_sorted.bam --indv pyg_fem_Z23346_GATG_sorted.bam --indv pyg_fem_Z23347_TCAG_sorted.bam --indv pyg_fem_Z23348_TGCGA_sorted.bam --indv pyg_fem_Z23351_CTAGG_sorted.bam --indv pyg_fem_Z23352_ACAAA_sorted.bam --indv pyg_fem_Z23354_AGCCG_sorted.bam --indv pyg_fem_Z23355_GTATT_sorted.bam --indv pyg_fem_Z23357_ACCGT_sorted.bam --indv pyg_fem_Z23358_GCTTA_sorted.bam --indv pyg_fem_Z23360_AGGAT_sorted.bam --indv pyg_fem_Z23361_ATTGA_sorted.bam --indv pyg_fem_Z23363_CCTAG_sorted.bam --indv pyg_fem_Z23364_GAGGA_sorted.bam --indv pyg_fem_Z23367_TAATA_sorted.bam --indv pyg_fem_Z23369_TCGTT_sorted.bam --indv pyg_fem_Z23370_GGTTGT_sorted.bam --indv pyg_fem_Z23371_CCACGT_sorted.bam --indv pyg_fem_Z23372_TTCAGA_sorted.bam --indv pyg_fem_Z23373_TAGGAA_sorted.bam --indv pyg_fem_Z23374_GCTCTA_sorted.bam --indv pyg_fem_Z23375_CCACAA_sorted.bam --indv pyg_fem_Z23377_GAGATA_sorted.bam --indv pyg_fem_Z23378_ATGCCT_sorted.bam --indv pyg_fem_Z23379_AGTGGA_sorted.bam --indv pyg_fem_Z23380_ACCTAA_sorted.bam --indv pyg_fem_Z23381_ATATGT_sorted.bam --indv pyg_fem_Z23382_ATCGTA_sorted.bam --indv pyg_fem_Z23383_CATCGT_sorted.bam --indv pyg_fem_Z23384_CGCGGT_sorted.bam --indv pyg_fem_Z23385_CTATTA_sorted.bam --indv pyg_fem_Z23386_GCCAGT_sorted.bam --out pygm_nucleotide_diversity_females
```

```
vcftools --vcf DB__Chr8L_out.vcf --site-pi --positions SNP_list.txt --indv pyg_mal_Z23337_CTCG_sorted.bam --indv pyg_mal_Z23339_ACTA_sorted.bam --indv pyg_mal_Z23349_CGCTT_sorted.bam --indv pyg_mal_Z23350_TCACG_sorted.bam --indv pyg_mal_Z23353_TTCTG_sorted.bam --indv pyg_mal_Z23356_CTGTA_sorted.bam --indv pyg_mal_Z23359_GGTGT_sorted.bam --indv pyg_mal_Z23362_CATCT_sorted.bam --indv pyg_mal_Z23365_GGAAG_sorted.bam --indv pyg_mal_Z23366_GTCAA_sorted.bam --indv pyg_mal_Z23368_TACAT_sorted.bam --indv pyg_mal_Z23376_CTTCCA_sorted.bam --out pygm_nucleotide_diversity_males
```
# fischbergi
```
/home/ben/projects/rrg-ben/ben/2022_GBS_lotsofxennies/individual_gvcfs_by_species/2023_fischbergi/mapped_to_XLv10_unfiltered
```
```
vcftools --vcf DB__Chr3L_out.vcf --site-pi --positions SNP_list.txt --indv ./fisch_Z23798_F_TAGCATGG_sorted --indv ./fisch_Z23799_F_TAGGCCAT_sorted --indv ./fisch_Z23774_F_CGCGATAT_TATTTTT_sorted --indv ./fisch_Z23777_F_AAAAGTT_CGCCTTAT_sorted --indv ./fisch_Z23778_F_GAATTCA_sorted --indv ./fisch_Z23779_F_GAACTTG_sorted --indv ./fisch_Z23780_F_GGACCTA_sorted --indv ./fisch_Z23782_F_AACGCCT_ACAGGGA_sorted --indv ./fisch_Z23783_F_AATATGG_sorted --indv ./fisch_Z23784_F_ACGTGTT_sorted --indv ./fisch_Z23785_F_ATTAATT_sorted --indv ./fisch_Z23787_F_CATAAGT_sorted --indv ./fisch_Z23788_F_CGCTGAT_sorted --indv ./fisch_Z23789_F_CGGTAGA_sorted --indv ./fisch_Z23790_F_CTACGGA_sorted --indv ./fisch_Z23792_F_TAGCGGA_sorted --indv ./fisch_Z23793_F_TCGAAGA_sorted --indv ./fisch_Z23794_F_TCTGTGA_sorted --indv ./fisch_Z23795_F_TGCTGGA_sorted --out fisc_nucleotide_diversity_females
```
```
vcftools --vcf DB__Chr3L_out.vcf --site-pi --positions SNP_list.txt --indv ./fisch_Z23786_M_ATTGGAT_sorted --indv ./fisch_Z23791_M_GCGGAAT_sorted --indv ./fisch_Z23781_M_AACCGAGA_GTCGATT_sorted --indv ./fisch_Z23796_M_ACGACTAG_sorted --indv ./fisch_Z23771_M_GTTGAA_TGGTACGT_sorted --indv ./fisch_Z23772_M_TAACGA_TCTCAGTG_sorted --indv ./fisch_Z23773_M_TGGCTA_sorted --indv ./fisch_Z23775_M_CTTGCTT_sorted --indv ./fisch_Z23776_M_ATGAAAG_sorted --out fisc_nucleotide_diversity_males
```

# allofraseri family 1
```
/home/ben/projects/rrg-ben/ben/2022_GBS_lotsofxennies/individual_gvcfs_by_species/2022_allofraseri/mapped_to_XLv10_concatscaf
```
```
vcftools --vcf DB__Chr7L_out.vcf --site-pi --positions SNP_list.txt  --indv allo_fem1_CAS265119_CGCTT_sorted.bam --indv allo_fem1_CAS265165_TCAG_sorted.bam --indv allo_fem1_Z23697_AGCG_sorted.bam --indv allo_fem1_Z23698_GTAA_sorted.bam --indv allo_fem1_Z23699_CGAT_sorted.bam --indv allo_fem1_Z23705_TGCA_sorted.bam --indv allo_fem1_Z23707_CTCG_sorted.bam --indv allo_fem1_Z23708_TTCCTGGA_sorted.bam --indv allo_fem1_Z23709_TATCGGGA_sorted.bam --indv allo_fem1_Z23710_GTGAGGGT_sorted.bam --indv allo_fem1_Z23711_GGATTGGT_sorted.bam --indv allo_fem1_Z23712_GCTGTGGA_sorted.bam --indv allo_fem1_Z23720_CGCCTTAT_sorted.bam --indv allo_fem1_Z23724_TGCAAGGA_sorted.bam --indv allo_fem1_Z23725_TAGGCCAT_sorted.bam --indv allo_fem1_Z23726_TAGCATGG_sorted.bam --indv allo_fem1_Z23727_ACGACTAG_sorted.bam --indv allo_fem1_Z23728_TGCTGGA_sorted.bam --indv allo_fem1_Z23730_TCGAAGA_sorted.bam --indv allo_fem1_Z23731_TAGCGGA_sorted.bam --indv allo_fem1_Z23732_GCGGAAT_sorted.bam --indv allo_fem1_Z23733_CTACGGA_sorted.bam --indv allo_fem1_Z23734_CGGTAGA_sorted.bam --indv allo_fem1_Z23736_CATAAGT_sorted.bam --indv allo_fem1_Z23737_ATTGGAT_sorted.bam --indv allo_fem1_Z23740_AATATGG_sorted.bam --indv allo_fem2_Z23702_CAGA_sorted.bam --indv allo_fem2_Z23704_ACTA_sorted.bam --indv allo_fem2_Z23714_CGTGTGGT_sorted.bam --indv allo_fem2_Z23718_ACAGGGA_sorted.bam --indv allo_fem2_Z23721_CGCGATAT_sorted.bam --out allo_fam1_fam2_nucleotide_diversity_females
```
```
vcftools --vcf DB__Chr7L_out.vcf --site-pi --positions SNP_list.txt --indv allo_mal1_CAS265071_TCACG_sorted.bam --indv allo_mal1_CAS265072_CTAGG_sorted.bam --indv allo_mal1_CAS265118_TGCGA_sorted.bam --indv allo_mal1_Z23696_GATG_sorted.bam --indv allo_mal1_Z23722_TCTCAGTG_sorted.bam --indv allo_mal1_Z23723_TGGTACGT_sorted.bam --indv allo_mal1_Z23729_TCTGTGA_sorted.bam --indv allo_mal1_Z23735_CGCTGAT_sorted.bam --indv allo_mal1_Z23738_ATTAATT_sorted.bam --indv allo_mal1_Z23739_ACGTGTT_sorted.bam --indv allo_mal1_Z23741_AACGCCT_sorted.bam --indv allo_mal1_Z23742_GTCGATT_sorted.bam --indv allo_mal1_Z23743_GGACCTA_sorted.bam --indv allo_mal1_Z23749_GAACTTG_sorted.bam --indv allo_mal2_Z23700_GCGT_sorted.bam --indv allo_mal2_Z23701_AACT_sorted.bam --indv allo_mal2_Z23715_CGCGGAGA_sorted.bam --indv allo_mal2_Z23716_CCATGGGT_sorted.bam --indv allo_mal2_Z23717_ACGTGGTA_sorted.bam --indv allo_mal2_Z23719_AACCGAGA_sorted.bam --out allo_fam1_fam2_nucleotide_diversity_males
```

