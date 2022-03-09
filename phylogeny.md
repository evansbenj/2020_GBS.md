# Phylogeny

This is a great tool to convert vcf files to fasta/nexus/phylip files, including IUPAC symbols:
```
https://github.com/edgardomortiz/vcf2phylip
```
in this directory:
```
/home/ben/projects/rrg-ben/ben/2020_GBS_muel_fish_allo_cliv_laev/raw_data/cutaddapted_by_species_across_three_plates/laevis/bin/vcf2phylip
```
with this command:
```
python vcf2phylip.py -i ../../vcfs_after_filtering_and_removal_with_XG/DB_newnew_chr9_10L_out.vcf_filtered.vcf.gz_filtered_removed.vcf  -o XGL713_123.fq.gz_sorted.bam -n --min-samples-locus 50
```
This makes a nexus file out of sites that have data for at least 50 individuals (out of 95 total, including the XG outgroup, which is also defined with the -o flag)
