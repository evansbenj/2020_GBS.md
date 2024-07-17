# Fst windows with angsd

To avoid errors I had to make new bam files that had only the chromosome of interest.

```
/home/ben/projects/rrg-ben/ben/2022_GBS_lotsofxennies/ben_scripts/2024_select_one_chr_from_multiple_bams.sh
```

```
#!/bin/sh
#SBATCH --job-name=samtools
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --mem=2gb
#SBATCH --output=samtools.%J.out
#SBATCH --error=samtools.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this
# sbatch 2024_select_one_chr_from_multiple_bams.sh chr

module load StdEnv/2023  gcc/12.3 samtools/1.20

for file in *_sorted_rg.bam
do
    samtools view -b ${file} ${1} > ${file}_${1}.bam
    samtools index ${file}_${1}.bam
done
```

Now make a text file for each sex with the names of all the bam files in it (e.g. fem_bamz_Chr3L.txt, mal_bamz_Chr3L.txt)

Now run this script on each file:
```
/home/ben/projects/rrg-ben/ben/2022_Liberia/ben_scripts/2024_angsd_fst_step1.sh
```
like this:
```
sbatch /home/ben/projects/rrg-ben/ben/2022_Liberia/ben_scripts/2024_angsd_fst_step1.sh /home/ben/projects/rrg-ben/ben/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa fem_bamz_Chr3L
```
Now load angsd
```
module load StdEnv/2023 angsd/0.940
```
now step2:
```
realSFS fem_bamz_Chr3L_.saf.idx mal_bamz_Chr3L_.saf.idx >fem_mal_Chr3L.ml
```
now step3:
```
realSFS fst index fem_bamz_Chr3L.saf.idx mal_bamz_Chr3L.saf.idx -sfs fem_mal_Chr3L.ml -fstout fem_mal_Chr3L__
```
now step4:
```
realSFS fst stats2 fem_mal_Chr3L__.fst.idx -win 50000 -step 10000 > slidingwindow
```
