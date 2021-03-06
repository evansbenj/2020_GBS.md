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

Once this is done I can summarize it using a file that specifies whether the individuals are male or female. I turned off Bonferrini correction (-C) and made the pvalue more stringent (-S 0.0001). This worked for fischbergi, which has a massive sex linked region.

```
/home/ben/projects/rrg-ben/ben/2020_radsex/bin/radsex distrib --markers-table /home/ben/projects/rrg-ben/ben/2020_GBS_muel_fish_allo_cliv_laev/raw_data/cutaddapted_by_species_across_three_plates/2020_fisc_catR1R2_markers_table.tsv --output-file /home/ben/projects/rrg-ben/ben/2020_GBS_muel_fish_allo_cliv_laev/raw_data/cutaddapted_by_species_across_three_plates/2020_fisc_distribution.tsv --popmap /home/ben/projects/rrg-ben/ben/2020_GBS_muel_fish_allo_cliv_laev/raw_data/cutaddapted_by_species_across_three_plates/fischbergi_sex_R1R2cat --min-depth 5 --groups M,F -C -S 0.0001
```
Or output a fasta file with significant reads:

```
/home/ben/projects/rrg-ben/ben/2020_radsex/bin/radsex signif --markers-table /home/ben/projects/rrg-ben/ben/2020_GBS_muel_fish_allo_cliv_laev/raw_data/cutaddapted_by_species_across_three_plates/2020_fisc_catR1R2_markers_table.tsv --output-file /home/ben/projects/rrg-ben/ben/2020_GBS_muel_fish_allo_cliv_laev/raw_data/cutaddapted_by_species_across_three_plates/2020_fisc_significant_markers.fasta --popmap /home/ben/projects/rrg-ben/ben/2020_GBS_muel_fish_allo_cliv_laev/raw_data/cutaddapted_by_species_across_three_plates/fischbergi_sex_R1R2cat --min-depth 5 --groups M,F --output-fasta -C -S 0.0001
```

And this can be done for clivii, muelleri, fischbergi, and parafraseri. 

I can print the output like this:
```R
#https://github.com/SexGenomicsToolkit/radsex


setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/submitted/XB_sex_determining_gene/RADsex")

library(sgtr)
library(ggplot2)

distrib_plot <- radsex_distrib("clivii_short_distribution.tsv",
                               groups = c("M", "F"),
                               group_labels = c("Males", "Females"),
                               title = "Distribution of markers",
                               significance_color = "green",
                               bins = c(0,1,5,10,20,35,100,1000),
                               colors = c("white", "blue"))

ggsave("Rplot_clivii_short.pdf", width = 5, height = 4, dpi = 200)
```


It is also possible to map the results and print this by chr, including the probability of association (P), and the sex bias (positive for male biased, Y-linked and negative for female biased, W-linked:
```
#!/bin/sh
#SBATCH --job-name=radsex_map
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=36:00:00
#SBATCH --mem=32gb
#SBATCH --output=radsex_map.%J.out
#SBATCH --error=radsex_map.%J.err
#SBATCH --account=def-ben

/home/ben/projects/rrg-ben/ben/2020_radsex/bin/radsex map --markers-file /home/ben/projects/rrg-ben/ben/2020_GBS_muel_fish_a
llo_cliv_laev/raw_data/cutaddapted_by_species_across_three_plates/2020_para_catR1R2_markers_table.tsv --output-file /home/be
n/projects/rrg-ben/ben/2020_GBS_muel_fish_allo_cliv_laev/raw_data/cutaddapted_by_species_across_three_plates/2020_para_catR1
R2_alignment_results.tsv --popmap /home/ben/projects/rrg-ben/ben/2020_GBS_muel_fish_allo_cliv_laev/raw_data/cutaddapted_by_s
pecies_across_three_plates/parafraseri_sex_R1R2cat --genome-file ~/projects/rrg-ben/ben/2020_XL_v9.2_refgenome/XENLA_9.2_gen
ome.fa.gz --min-quality 20 --min-frequency 0.1 --min-depth 1 --groups M,F

```

Get rid of scaffolds like this:
```
sed '/Sca/d' /home/ben/projects/rrg-ben/ben/2020_radsex/bin/2017_cliv_markers_alignment_all.tsv > /home/ben/projects/rrg-ben/ben/2020_radsex/bin/2017_cliv_markers_alignment_all_chronly.tsv
```

The "radsex_map_circos" in the sgtr package is broken, so I just printed the results by chr using ggplot:
```
setwd("/Users/Shared/Previously Relocated Items/Security/projects/2020_fisch_clivii_muel_allo_RADseq/RADsex")

library(ggplot2)

data_SL<-read.table("./2020_allo_catR1R2_alignment_results_chronly.tsv",header=T)
data_SL_sorted <- data_SL[order(data_SL["Contig"],data_SL["Position"]),]
head(data_SL_sorted)
data_SL_sorted$Position=as.numeric(data_SL_sorted$Position)

sp <- ggplot(data_SL_sorted, aes(x=Position, y=-log(P))) + 
  geom_point(size=1) + 
  facet_wrap( ~ Contig, ncol=2)
sp

sp <- ggplot(data_SL_sorted, aes(x=Position, y=Bias)) + 
  # Bias is positive if male specific
  # neggative if female specific
  geom_point(size=1) + 
  facet_wrap( ~ Contig, ncol=2)
sp
```
