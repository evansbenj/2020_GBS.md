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
now step2:
```
realSFS fem_bamz_.saf.idx mal_bamz_.saf.idx >fem.mal.ml
```
now step3:
```
realSFS fst index pop1.saf.idx pop2.saf.idx -sfs pop1.pop2.ml -fstout here
```
now step4:
```
realSFS fst stats2 here.fst.idx -win 50000 -step 10000 >slidingwindow
```
