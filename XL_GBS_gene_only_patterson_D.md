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
# Phase
```
sbatch Beagle.sh  ../vcfs_after_filtering_and_removal_with_XG/allchr_filtered_removed_subgenomeLonly_genicwith10kbuffer_only_chr1Lonly.vcf.gz
```
and so on for each chr

# Make geno files
```
python /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_depth_3sigmas/final_data_including_sites_with_lots_of_missing_data/genomics_general/VCF_processing/parseVCF.py -i ../vcfs_after_filtering_and_removal_with_XG/allchr_filtered_removed_subgenomeSonly_genicwith10kbuffer_only_chr9_10Sonly.vcf.gz -o chr9_10S__genicwith10kbuffer_only.geno.gz
```
and so on for each chr
