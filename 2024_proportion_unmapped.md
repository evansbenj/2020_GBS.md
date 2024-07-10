# Calculate proportion of unmapped reads

This script is located here:
```
/home/ben/projects/rrg-ben/ben/2022_GBS_lotsofxennies/ben_scripts
```
and will calculate the proportion of unmapped reads of all bam files in a directory and save it in a file called "unmapped.txt"
```
#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=8gb
#SBATCH --output=samtools_unmapped.%J.out
#SBATCH --error=samtools_unmapped.%J.err
#SBATCH --account=rrg-ben

# run by passing the path to the sorted bam files like this
# sbatch ./2021_samtools_index_bamfiles.sh pathtorefgenome
# sbatch ./2021_samtools_index_bamfiles.sh /home/ben/projects/rrg-ben/ben/2020_GBS_muel_fish_allo_cliv_laev/raw_data/cutaddapted_by_species_across_three_plates/clivii/

module load StdEnv/2020 samtools/1.12

echo "Proportion_unmapped" > unmapped.txt

for file in ${1}*.bam
do
    echo ${file} >> unmapped.txt
    samtools view -c ${file} >> unmapped.txt
    samtools view -c -F 0X04 ${file} >> unmapped.txt
done
```
