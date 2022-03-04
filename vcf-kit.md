In this directory on graham:
```
/home/ben/projects/rrg-ben/ben/2020_GBS_muel_fish_allo_cliv_laev/raw_data/cutaddapted_by_species_across_three_plates/laevis/bin
```
Make a nj tree:
```
module load gcc bwa samtools bcftools blast muscle primer3 python/3
create a virtual environment and activate it
virtualenv --no-download vcf_env
source vcf_env/bin/activate
vcf_env/bin/vk phylo tree nj ../allchr_filtered_removed.vcf > allchr_filtered_removed_with_XG.tre
```
