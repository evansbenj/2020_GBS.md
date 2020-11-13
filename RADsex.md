#RADsex

I'm first going to try this approach on the raw reads: https://github.com/SexGenomicsToolkit/radsex

Issues to keep in mind are that the method really is aimed at single read data.  I can probably get around this either by concatenating forward and reverse reads just for this analysis or (and this is preferable) making the table for everything in the first step and then doing the second (distrib) step using different popmap files (one for the forward reads and one for the reverse).

RADsex is installed here:
```
/home/ben/projects/rrg-ben/ben/2020_radsex/bin
```

First step is to make a table of reads.  This uses up some memory:
```
#!/bin/sh
#SBATCH --job-name=radsex
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --mem=32gb
#SBATCH --output=radsex.%J.out
#SBATCH --error=radsex.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this
# sbatch ./2020_fastqc.sh ../raw_data/plate1
/home/ben/projects/rrg-ben/ben/2020_radsex/bin/radsex process --input-dir /home/ben/projects/rrg-ben/ben/2020_GBS_muel_fish_allo_c
liv_laev/raw_data/cutaddapted_by_species_across_three_plates/allofraseri --output-file /home/ben/projects/rrg-ben/ben/2020_GBS_mue
l_fish_allo_cliv_laev/raw_data/cutaddapted_by_species_across_three_plates/2020_allo_markers_table.tsv --threads 16 --min-depth 1
```

Once this is done I can summarize it using a file that specifies whether the individuals are male or female.

```
/home/ben/projects/rrg-ben/ben/2020_radsex/bin/radsex distrib --markers-table /home/ben/projects/rrg-ben/ben/2020_GBS_muel_fish_allo_cliv_laev/raw_data/cutaddapted_by_species_across_three_plates/2020_allo_markers_table.tsv --output-file /home/ben/projects/rrg-ben/ben/2020_GBS_muel_fish_allo_cliv_laev/raw_data/cutaddapted_by_species_across_three_plates/2020_allo_distribution.tsv --popmap /home/ben/projects/rrg-ben/ben/2020_GBS_muel_fish_allo_cliv_laev/raw_data/cutaddapted_by_species_across_three_plates_allofraseri_sex --min-depth 5 --groups M,F
```

And this can be done for clivii, muelleri, fischbergi, and parafraseri. 
