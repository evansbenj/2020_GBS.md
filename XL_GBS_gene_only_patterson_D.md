# Make a bed file with only coordinates of genes on known L or S chrs:
locally:
```
cat XLv9.2_xenbase_annotations.gff | grep 'gene' | grep 'chr' | cut -f1,4,5 > gene_only_bed.bed
```
# Now use excel to add a buffer of 10,000 bp on either side

# Make beds plus buffer for L and S subgenomes only
```
grep 'L' gene_only_plusminus_10000bp_bed.bed > gene_only_plusminus_10000bp_Lsubgenome_only.bed
gene_only_bed ben$ grep 'S' gene_only_plusminus_10000bp_bed.bed > gene_only_plusminus_10000bp_Ssubgenome_only.bed
```
# Extract SNPs from each subgenome vcf
Path:
```
/home/ben/projects/rrg-ben/ben/2020_GBS_muel_fish_allo_cliv_laev/raw_data/cutaddapted_by_species_across_three_plates/laevis/vcfs_after_filtering_and_removal_with_XG
```

```
module load StdEnv/2020  gcc/9.3.0 bcftools/1.13

bcftools view -R gene_only_plusminus_10000bp_Lsubgenome_only.bed allchr_filtered_removed_subgenomeLonly.vcf.gz -Oz -o allchr_filtered_removed_subgenomeLonly_genicwith10kbuffer_only.vcf.gz

bcftools view -R gene_only_plusminus_10000bp_Ssubgenome_only.bed allchr_filtered_removed_subgenomeSonly.vcf.gz -Oz -o allchr_filtered_removed_subgenomeSonly_genicwith10kbuffer_only.vcf.gz
```
# Extract individual chrs
```
bcftools view -r chr1L allchr_filtered_removed_subgenomeLonly_genicwith10kbuffer_only.vcf.gz -Oz -o allchr_filtered_removed_subgenomeLonly_genicwith10kbuffer_only_chr1Lonly.vcf.gz
```
