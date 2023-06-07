# AngSD for association tests

Path
```
/home/ben/projects/rrg-ben/ben/2022_GBS_lotsofxennies/individual_gvcfs_by_species/2022_allofraseri/mapped_to_XLv10_concatscaf
```

Example sbatch script:
```
#!/bin/sh
#SBATCH --job-name=angsd_allo
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=256gb
#SBATCH --output=angsd_allo.%J.out
#SBATCH --error=angsd_allo.%J.err
#SBATCH --account=def-ben


module load StdEnv/2020 angsd/0.939

angsd -yBin bin_sex.ybin -doAsso 1 -GL 1 -out out_additive_F1 -doMajorMinor 1 -doMaf 1 -SNP_pval 1e-6 -minInd 4 -ba
m bam.filelist -P 5 -doCounts 1 -setMinDepthInd 2 -setMaxDepthInd 100 -Pvalue 1
```
Filter to remove non-significant sites for plotting. 
Without the "-Pvalue 1' flag, the 6th column is a chisq value with df=1, so lets save only significant values (higher than 7):
```
zcat tempty.lrt0.gz | awk '$6 < 7 { next } { print }'> sig.only
zcat out_additive_F1.lrt0.gz | awk '$7 > 0.001 { next } { print }'> 2022_pyg_P_gt_0.001_only.txt
```
