# XL mtDNA
I am going to estimate a mtDNA phylogeny using iqtree

I have the executable here on graham (version 1.6.12):
```
/home/ben/projects/rrg-ben/ben/2020_GBS_muel_fish_allo_cliv_laev/bin/iqtree-1.6.12-Linux/bin/iqtree
```

# Input nexus file
```
/home/ben/projects/rrg-ben/ben/2020_GBS_muel_fish_allo_cliv_laev/XL_mtDNA/Nucleotide_alignment_X.laevis_nodups_alllocalitiesknown.nex
```
```
#!/bin/sh
#SBATCH --job-name=iqtree
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=8:00:00
#SBATCH --mem=2gb
#SBATCH --output=iqtree.%J.out
#SBATCH --error=iqtree.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this
# sbatch ./2022_iqtree.sh path_and_filename
/home/ben/projects/rrg-ben/ben/2020_GBS_muel_fish_allo_cliv_laev/bin/iqtree-1.6.12-Linux/bin/iqtree -s ${1} -m TEST -nt 1 -pr
e ${1}_
```
