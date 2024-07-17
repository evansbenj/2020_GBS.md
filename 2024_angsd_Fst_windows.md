# Fst windows with angsd

First make a text file for each sex with the names of all the bam files in it (e.g. fem_bamz.txt, mal_bamz.txt)

Now run this script on each file:
```
/home/ben/projects/rrg-ben/ben/2022_Liberia/ben_scripts/2024_angsd_fst_step1.sh
```
like this:
```
sbatch /home/ben/projects/rrg-ben/ben/2022_Liberia/ben_scripts/2024_angsd_fst_step1.sh /home/ben/projects/rrg-ben/ben/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa fem_bamz Chr3L.bed
```
Now load angsd
```
module load StdEnv/2023 angsd/0.940
```
