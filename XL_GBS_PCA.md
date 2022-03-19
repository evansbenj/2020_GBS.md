# Final XL file

In this directory:
```
/home/ben/projects/rrg-ben/ben/2020_GBS_muel_fish_allo_cliv_laev/raw_data/cutaddapted_by_species_across_three_plates/laevis/vcfs_after_filtering_and_removal_with_XG
```

# For PCA, get rid of non-XL samples and divide file up by subgenome:
```
bcftools concat DB_newnew_chr1L_out.vcf_filtered.vcf.gz_filtered_removed.vcf DB_newnew_chr1S_out.vcf_filtered.vcf.gz_filtered_removed.vcf DB_newnew_chr2L_out.vcf_filtered.vcf.gz_filtered_removed.vcf DB_newnew_chr2S_out.vcf_filtered.vcf.gz_filtered_removed.vcf DB_newnew_chr3L_out.vcf_filtered.vcf.gz_filtered_removed.vcf DB_newnew_chr3S_out.vcf_filtered.vcf.gz_filtered_removed.vcf DB_newnew_chr4L_out.vcf_filtered.vcf.gz_filtered_removed.vcf DB_newnew_chr4S_out.vcf_filtered.vcf.gz_filtered_removed.vcf DB_newnew_chr5L_out.vcf_filtered.vcf.gz_filtered_removed.vcf DB_newnew_chr5S_out.vcf_filtered.vcf.gz_filtered_removed.vcf DB_newnew_chr6L_out.vcf_filtered.vcf.gz_filtered_removed.vcf DB_newnew_chr6S_out.vcf_filtered.vcf.gz_filtered_removed.vcf DB_newnew_chr7L_out.vcf_filtered.vcf.gz_filtered_removed.vcf DB_newnew_chr7S_out.vcf_filtered.vcf.gz_filtered_removed.vcf DB_newnew_chr8L_out.vcf_filtered.vcf.gz_filtered_removed.vcf DB_newnew_chr8S_out.vcf_filtered.vcf.gz_filtered_removed.vcf DB_newnew_chr9_10L_out.vcf_filtered.vcf.gz_filtered_removed.vcf DB_newnew_chr9_10S_out.vcf_filtered.vcf.gz_filtered_removed.vcf -o all_XL_withXG.vcf.gz -O z   
```
```
bcftools view -Ou --samples "^RT5_Botsw_GGATTGGT_cuttrim_sorted.bam,XGL713_123.fq.gz_sorted.bam,amnh17260_Nigeria_GTGAGGGT_cuttrim_sorted.bam,BJE2897_Lendu_TTCCTGGA_cuttrim_sorted.bam,BJE3252_Cameroon_TATCGGGA_cuttrim_sorted.bam" all_XL_withXG.vcf.gz -o all_XL_onlyXL.vcf.gz -O z
```
```
bcftools view all_XL_onlyXL.vcf.gz --regions chr1L,chr2L,chr3L,chr4L,chr5L,chr6L,chr7L,chr8L,chr9_10L -o all_XL_onlyXL_only_Lsubgenome.vcf.gz -O z
```
```
bcftools view all_XL_onlyXL.vcf.gz --regions chr1S,chr2S,chr3S,chr4S,chr5S,chr6S,chr7S,chr8S,chr9_10S -o all_XL_onlyXL_only_Ssubgenome.vcf.gz -O z
```
