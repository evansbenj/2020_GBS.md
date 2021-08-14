# Create geno files

For each chromosome:
```
module load StdEnv/2020
module load scipy-stack/2020b
module load python/3.8.2
python /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_depth_3sigmas/final_data_including_sites_with_lots_of_missing_data/genomics_general/VCF_processing/parseVCF.py -i DB_chr1L_out.vcf -o chr1L.geno.gz
```
