# Calculate proportion of mapped reads

This script is located here:
```
/home/ben/projects/rrg-ben/ben/2022_GBS_lotsofxennies/ben_scripts
```
and is modified from suggestions here:
```
https://sarahpenir.github.io/bioinformatics/awk/calculating-mapping-stats-from-a-bam-file-using-samtools-and-awk/
```
It will calculate the proportion of unmapped reads of all bam files in a directory and save it in a file called "mapped_flagstat.txt"
```
#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=8gb
#SBATCH --output=samtools_mapped_flagstat.%J.out
#SBATCH --error=samtools_mapped_flagstat.%J.err
#SBATCH --account=rrg-ben

module load StdEnv/2023  gcc/12.3 samtools/1.20

echo "Proportion_mapped" > unmapped_flagstat.txt

for file in ${1}*_sorted.bam_rg.bam
do
    echo ${file} >> mapped_flagstat.txt
    samtools flagstat ${file} | awk -F "[(|%]" 'NR == 8 {print $2}' >> mapped_flagstat.txt
done

```
