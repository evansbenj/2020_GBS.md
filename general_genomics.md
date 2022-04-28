# Phase vcf files
```
#!/bin/sh
#SBATCH --job-name=beagle
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=3:00:00
#SBATCH --mem=128gb
#SBATCH --output=beagle.%J.out
#SBATCH --error=beagle.%J.err
#SBATCH --account=def-ben

# sbatch Beagle.sh chr

module load java

java -Xmx12g -jar /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_dept
h_3sigmas/final_data_including_sites_with_lots_of_missing_data/twisst/beagle.18May20.d20.jar gt=${1} out=${1}_phased.vcf.gz impu
te=true
```

# Create geno files

For each chromosome:
```
module load StdEnv/2020
module load scipy-stack/2020b
module load python/3.8.2
python /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_depth_3sigmas/final_data_including_sites_with_lots_of_missing_data/genomics_general/VCF_processing/parseVCF.py -i DB_chr1L_out.vcf_phased.vcf.gz.vcf.gz -o chr1L.geno.gz



python /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_depth_3sigmas/final_data_including_sites_with_lots_of_missing_data/genomics_general/VCF_processing/parseVCF.py -i DB_chr1S_out.vcf_phased.vcf.gz.vcf.gz -o chr1S.geno.gz
```
and so on for each chr


Then it is necessary to swap any astrisks with Ns:
```
gunzip chr18.geno.gz
sed -i 's/\*/N/g' chr18.geno 
gzip -c chr18.geno > chr18.geno.gz
```

for autosomes:
```
#!/bin/sh
#SBATCH --job-name=abba
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=8gb
#SBATCH --output=abba.%J.out
#SBATCH --error=abba.%J.err
#SBATCH --account=def-ben

#sbatch ABBABABA.sh chr P1 P2 P3 O
#sbatch ABBABABA.sh chr1L blue green yellow others

#populations
#blue green yellow oout

module load StdEnv/2020
module load scipy-stack/2020b
module load python/3.8.2

echo python3 /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_depth_3sigmas/final_data_including_sites_with_lots_of_missing_data/genomics_general/ABBABABAwindows.py -g ../vcfs_after_filtering_and_removal/${1}.geno.gz -f phased -o ../vcfs_after_filtering_and_removal/${1}_${2}_${3}_${4}_${5}.csv -w 2000000 -m 100 -s 2000000 -P1 ${2} -P2 ${3} -P3 ${4} -O ${5} -T 10 --minData 0.5 --popsFile pops_all.txt --writeFailedWindows --windType coordinate

python3 /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_depth_3sigmas/final_data_including_sites_with_lots_of_missing_data/genomics_general/ABBABABAwindows.py -g ../vcfs_after_filtering_and_removal/${1}.geno.gz -f phased -o ../vcfs_after_filtering_and_removal/${1}_${2}_${3}_${4}_${5}.csv -w 2000000 -m 100 -s 2000000 -P1 ${2} -P2 ${3} -P3 ${4} -O ${5} -T 10 --minData 0.5 --popsFile pops_all.txt --writeFailedWindows --windType coordinate
```
sometimes out of memory errors arose.  This could be partially dealt with by echoing the sbatch command and then copying and pasting the commandline after loading modules. 
