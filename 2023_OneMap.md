# Sex-specific recombination landscapes

Best to work with filtered files
directory:
```
/home/ben/projects/rrg-ben/ben/2022_GBS_lotsofxennies/individual_gvcfs_by_species/2022_allofraseri/2018_2022_genotyped_after_filtering
```
for allofraseri I had to change the name of two individuals that previously had the wrong sex:
```
sed -i 's/allo_fem1_Z23696_GATG_cuttrim_sorted.bam/allo_mal1_Z23696_GATG_cuttrim_sorted.bam/g' DB_newnew_chr9_10S_genotyped.vcf_filtered.vcf.gz_filtered_removed.vcf

sed -i 's/allo_fem2_Z23719_AACCGAGA_cuttrim_sorted.bam/allo_mal2_Z23719_AACCGAGA_cuttrim_sorted.bam/g' DB_newnew_chr9_10S_genotyped.vcf_filtered.vcf.gz_filtered_removed.vcf


```

Now filter individuals that are from one family:

```
module load bcftools
bcftools view -S allo_family_one.txt infile > outfile
bcftools view -S allo_family_one.txt DB_newnew_chr9_10S_genotyped.vcf_filtered.vcf.gz_filtered_removed.vcf > allo_family_one_chr9_10S_filtered.vcf
```

where `allo_family_one.txt` is:
```
allo_Cam_female_4_F_AGGAT_TAATA_cuttrim_sorted.bam
allo_Cam_male_1_M_GGTGT_GTCAA_cuttrim_sorted.bam
allo_fem1_CAS265119_CGCTT_cuttrim_sorted.bam
allo_fem1_CAS265165_TCAG_cuttrim_sorted.bam
allo_fem1_Z23697_AGCG_cuttrim_sorted.bam
allo_fem1_Z23698_GTAA_cuttrim_sorted.bam
allo_fem1_Z23699_CGAT_cuttrim_sorted.bam
allo_fem1_Z23705_TGCA_cuttrim_sorted.bam
allo_fem1_Z23707_CTCG_cuttrim_sorted.bam
allo_fem1_Z23708_TTCCTGGA_cuttrim_sorted.bam
allo_fem1_Z23709_TATCGGGA_cuttrim_sorted.bam
allo_fem1_Z23710_GTGAGGGT_cuttrim_sorted.bam
allo_fem1_Z23711_GGATTGGT_cuttrim_sorted.bam
allo_fem1_Z23712_GCTGTGGA_cuttrim_sorted.bam
allo_fem1_Z23720_CGCCTTAT_cuttrim_sorted.bam
allo_fem1_Z23724_TGCAAGGA_cuttrim_sorted.bam
allo_fem1_Z23725_TAGGCCAT_cuttrim_sorted.bam
allo_fem1_Z23726_TAGCATGG_cuttrim_sorted.bam
allo_fem1_Z23727_ACGACTAG_cuttrim_sorted.bam
allo_fem1_Z23728_TGCTGGA_cuttrim_sorted.bam
allo_fem1_Z23730_TCGAAGA_cuttrim_sorted.bam
allo_fem1_Z23731_TAGCGGA_cuttrim_sorted.bam
allo_fem1_Z23732_GCGGAAT_cuttrim_sorted.bam
allo_fem1_Z23733_CTACGGA_cuttrim_sorted.bam
allo_fem1_Z23734_CGGTAGA_cuttrim_sorted.bam
allo_fem1_Z23736_CATAAGT_cuttrim_sorted.bam
allo_fem1_Z23737_ATTGGAT_cuttrim_sorted.bam
allo_fem1_Z23740_AATATGG_cuttrim_sorted.bam
allo_mal1_CAS265071_TCACG_cuttrim_sorted.bam
allo_mal1_CAS265072_CTAGG_cuttrim_sorted.bam
allo_mal1_CAS265118_TGCGA_cuttrim_sorted.bam
allo_mal1_Z23696_GATG_cuttrim_sorted.bam
allo_mal1_Z23722_TCTCAGTG_cuttrim_sorted.bam
allo_mal1_Z23723_TGGTACGT_cuttrim_sorted.bam
allo_mal1_Z23729_TCTGTGA_cuttrim_sorted.bam
allo_mal1_Z23735_CGCTGAT_cuttrim_sorted.bam
allo_mal1_Z23738_ATTAATT_cuttrim_sorted.bam
allo_mal1_Z23739_ACGTGTT_cuttrim_sorted.bam
allo_mal1_Z23741_AACGCCT_cuttrim_sorted.bam
allo_mal1_Z23742_GTCGATT_cuttrim_sorted.bam
allo_mal1_Z23743_GGACCTA_cuttrim_sorted.bam
allo_mal1_Z23749_GAACTTG_cuttrim_sorted.bam
```
Then I compressed and indexed this file:
```
module load tabix
bgzip -c allo_family_one_chr9_10S_filtered.vcf > allo_family_one_chr9_10S_filtered.vcf.gz
tabix -p vcf allo_family_one_chr9_10S_filtered.vcf.gz
```

# OneMap

```R
#!/usr/bin/env Rscript 

# for computecanada first type this before starting R
# module load  StdEnv/2020 r/4.3.1
# export R_LIBS=~/.local/R/$EBVERSIONR/

## Working directory
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#   BiocManager::install("GenomicRanges")


# directions on passing variables to Rscript via unix
# https://stackoverflow.com/questions/56777529/how-to-pass-bash-variable-into-r-script

# Rscript "2023_OneMap.R" --args inputfile_C="allo_family_one_chr9_10S_filtered.vcf.gz" mom_C="allo_Cam_female_4_F_AGGAT_TAATA_cuttrim_sorted.bam" dad_C="allo_Cam_male_1_M_GGTGT_GTCAA_cuttrim_sorted.bam" prefix_C="allo_1_"

library(onemap)
library(ggplot2)
library(tidyverse)

rm(list=ls()) # removes all variables
#setwd("/Users/Shared/Previously Relocated Items/Security/projects/2022_GBS_lotsof_Xennies/2023_OneMap")


# Process the variables that were passed in from unix
arguments <- commandArgs(trailingOnly = TRUE) 
args <- strsplit(arguments, "=", fixed = TRUE)

for (e in args) {
  argname <- e[1]
  if (! is.na(e[2])) {
    argval <- e[2]
    ## regular expression to delete initial \" and trailing \"
    argval <- gsub("(^\\\"|\\\"$)", "", argval)
  }
  else {
    # If arg specified without value, assume it is bool type and TRUE
    argval <- TRUE
  }
  
  # Infer type from last character of argname, cast val
  type <- substring(argname, nchar(argname), nchar(argname))
  if (type == "I") {
    argval <- as.integer(argval)
  }
  if (type == "N") {
    argval <- as.numeric(argval)
  }
  if (type == "L") {
    argval <- as.logical(argval)
  }
  assign(argname, argval)
  cat("Assigned", argname, "=", argval, "\n")
}

# import data from a vcf file
# my_dat <- onemap_read_vcfR(vcf = "allo_family_one_chr9_10S_filtered.vcf.gz",
my_dat <- onemap_read_vcfR(vcf = inputfile_C,
                                    cross = c("outcross"),
                                    #parent1 = c("allo_Cam_female_4_F_AGGAT_TAATA_cuttrim_sorted.bam"), 
                                    parent1 = mom_C,
                                    #parent2 = c("allo_Cam_male_1_M_GGTGT_GTCAA_cuttrim_sorted.bam"), 
                                    parent2 = dad_C,
                                    only_biallelic = TRUE,
                                    verbose = TRUE); my_dat

# filter sites with lots of missing data
data_filtered <- filter_missing(my_dat, threshold = 0.25);data_filtered

# identify markers with or without segregation distortion
segreg_test <- test_segregation(data_filtered)
dist <- select_segreg(segreg_test, distorted = TRUE, numbers = TRUE) #to show the markers numbers with segregation distortion
dist
no_dist <- select_segreg(segreg_test, distorted = FALSE, numbers = TRUE) #to show the markers numbers without segregation distortion
no_dist

# Estimate the recombination fraction between all pairs of markers, 
# using two-point tests.
# rm_mks = T gets rid of ones with weird segregation
# due to excess of missing data or segregation deviation
# but this doesn't seem to really work, so I did it again in the next step
twopts <- rf_2pts(data_filtered, rm_mks = T); twopts # an object of class rf_2pts

# make a variable with the markers that do not have segregation distortion
mark_no_dist <- make_seq(twopts, c(no_dist))

# look at markers
marker_type(mark_no_dist)
p <- marker_type(mark_no_dist)
p %>% group_by(Type) %>% count()

# get indexes of double hets or maternal specific hets:
maternal_SNPs <- as.integer(p[(p$Type == 'B3.7')|(p$Type == 'D1.10'), ]$Marker)
# get indexes of double hets or paternal specific hets:
paternal_SNPs <- as.integer(p[(p$Type == 'B3.7')|(p$Type == 'D2.15'), ]$Marker)
# get indexes of mat and pat sites
# matpat_SNPs <- as.integer(p[(p$Type == 'B3.7')|(p$Type == 'D1.10')|(p$Type == 'D2.15'), ]$Marker)

# test to have no informative SNPs for mother (only pat het sites):
# maternal_SNPs <- as.integer(p[(p$Type == 'D2.15'), ]$Marker) # still shows recombination
# now try only maternal het sites
# maternal_SNPs <- as.integer(p[(p$Type == 'D1.10'), ]$Marker)

# Marker types are explained here:
# https://statgen-esalq.github.io/tutorials/onemap/Outcrossing_Populations.html

# 1 B3.7    321  ab×ab (double het)
# 2 D1.10   483  ab×aa (first parent het; maternal)
# 3 D2.15   384  aa×ab (second parent het; paternal)

# make an object with only maternal SNPs
maternal_SNPs_no_dist <- make_seq(twopts, maternal_SNPs)

# make an object with only paternal SNPs
paternal_SNPs_no_dist <- make_seq(twopts, paternal_SNPs)

# make an object with only matpat SNPs
# matpat_SNPs_no_dist <- make_seq(twopts, matpat_SNPs)

maternal_map <- onemap::map(maternal_SNPs_no_dist)
paternal_map <- onemap::map(paternal_SNPs_no_dist)
# matpat_map <- onemap::map(matpat_SNPs_no_dist)

maternal_parents_haplot <- parents_haplotypes(maternal_map)
paternal_parents_haplot <- parents_haplotypes(paternal_map)
# matpat_haplot <- parents_haplotypes(matpat_map)

write.table(maternal_parents_haplot, paste(prefix_C,"mat_parents_haplot.txt",sep="_"))
write.table(paternal_parents_haplot, paste(prefix_C,"pat_parents_haplot.txt",sep="_"))
# write.table(maternal_haplot, "mat_haplot.txt")
# write.table(paternal_haplot, "pat_haplot.txt")
# write.table(matpat_haplot, "matpat_haplot.txt")

# draw_map(maternal_map)

# Export haplotypes
mat_progeny_haplot <- progeny_haplotypes(maternal_map, 
                                     most_likely = TRUE, 
                                     ind = c(1:40), 
                                     group_names = "chr9_10S")

pat_progeny_haplot <- progeny_haplotypes(paternal_map, 
                                         most_likely = TRUE, 
                                         ind = c(1:40), 
                                         group_names = "chr9_10S")

#matpat_progeny_haplot <- progeny_haplotypes(matpat_map, 
#                                         most_likely = TRUE, 
#                                         ind = c(1:40), 
#                                         group_names = "chr9_10S")

# write.table(mat_progeny_haplot, "mat_progeny_haplot.txt")
# write.table(pat_progeny_haplot, "pat_progeny_haplot.txt")
# write.table(matpat_progeny_haplot, "matpat_progeny_haplot.txt")


mat_progeny_haplot_wide <- mat_progeny_haplot %>%
  as_tibble()%>%
  select(-parents, -parents.homologs) %>%
  mutate(prob = round(prob, 4)) %>%
  pivot_wider(names_from = allele, values_from = prob) %>%
  arrange(grp, pos)

mat_progeny_haplot_wide <- mat_progeny_haplot_wide[order(mat_progeny_haplot_wide$ind, mat_progeny_haplot_wide$pos), ]

pat_progeny_haplot_wide <- pat_progeny_haplot %>%
  as_tibble()%>%
  select(-parents, -parents.homologs) %>%
  mutate(prob = round(prob, 4)) %>%
  pivot_wider(names_from = allele, values_from = prob) %>%
  arrange(grp, pos)

pat_progeny_haplot_wide <- pat_progeny_haplot_wide[order(pat_progeny_haplot_wide$ind, pat_progeny_haplot_wide$pos), ]

#write.table(mat_progeny_haplot_wide, "mat_progeny_haplot_wide.txt")
#write.table(pat_progeny_haplot_wide, "pat_progeny_haplot_wide.txt")

write.table(mat_progeny_haplot_wide, paste(prefix_C,"mat_progeny_haplot_wide.txt",sep="_"))
write.table(pat_progeny_haplot_wide, paste(prefix_C,"pat_progeny_haplot_wide.txt",sep="_"))


#mathap <- plot(mat_progeny_haplot, position = "stack")
#mathap_split <- plot(mat_progeny_haplot, position = "split")
#pathap <- plot(pat_progeny_haplot, position = "stack")
#matpathap <- plot(matpat_progeny_haplot, position = "stack")

#ggsave(file="mathap.pdf", mathap, width=10, height=30)
#ggsave(file="mathap_split.pdf", mathap_split, width=10, height=40)
#ggsave(file="pathap.pdf", pathap, width=10, height=30)
#ggsave(file="matpathap.pdf", matpathap, width=10, height=30)






# https://search.r-project.org/CRAN/refmans/onemap/html/progeny_haplotypes_counts.html
# Generate graphic with the number of break points for each individual 
# considering the most likely genotypes estimated by the HMM.
# e <- progeny_haplotypes_counts(progeny_haplotypes(maternal_map, 
#                                                  most_likely = TRUE, 
#                                                  ind = c(1:40), 
#                                                  group_names = "chr9_10S"))



```

Run like this:
```
sbatch ../../../ben_scripts/2023_Rscript_OneMap.sh allo_family_one_chr1S_filtered.vcf allo_Cam_female_4_F_AGGAT_TAATA_cuttrim_sorted.bam allo_Cam_male_1_M_GGTGT_GTCAA_cuttrim_sorted.bam allo_1_chr1S
```
where the sbatch file is:
```sh
#!/bin/sh
#SBATCH --job-name=Rscript
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=6:00:00
#SBATCH --mem=2gb
#SBATCH --output=Rscript.%J.out
#SBATCH --error=Rscript.%J.err
#SBATCH --account=def-ben

module load  StdEnv/2020 r/4.3.1
export R_LIBS=~/.local/R/$EBVERSIONR/

# execute in the directory that has the 2023_OneMap.R script and also the inputfile

Rscript "2023_OneMap.R" --args inputfile_C=${1} mom_C=${2} dad_C=${3} prefix_C=${4}
```

# Draft processing script for OneMap output
```R
setwd('/Users/Shared/Previously Relocated Items/Security/projects/2022_GBS_lotsof_Xennies/2023_OneMap')

# read in the OneMap data 
chr7S_mat <- read.table("./allo_1_chr7S_mat_progeny_haplot_wide.txt", header = T)
chr7S_pat <- read.table("./allo_1_chr7S_pat_progeny_haplot_wide.txt", header = T)
chr9_10S_mat <- read.table("./allo_1_9_10S_mat_progeny_haplot_wide.txt", header = T)
chr9_10S_pat <- read.table("./allo_1_chr9_10S_pat_progeny_haplot_wide.txt", header = T)

# make a column with the genomic coordinates
chr7S_mat$coord <- as.numeric(str_split_i(chr7S_mat$marker, "_",2))
chr7S_pat$coord <- as.numeric(str_split_i(chr7S_pat$marker, "_",2))
chr9_10S_mat$coord <- as.numeric(str_split_i(chr9_10S_mat$marker, "_",3)) # need to use a 3 because Chr9_10 has an underscore
chr9_10S_pat$coord <- as.numeric(str_split_i(chr9_10S_pat$marker, "_",3)) # need to use a 3 because Chr9_10 has an underscore

# make a column with the chr
chr7S_mat$CHR <- str_split_i(chr7S_mat$marker, "_",1)
chr7S_pat$CHR <- str_split_i(chr7S_pat$marker, "_",1)
chr9_10S_mat$CHR<- paste(str_split_i(chr9_10S_mat$marker, "_",1),"_",str_split_i(chr9_10S_mat$marker, "_",2), sep = "") # need to use a 3 because Chr9_10L has an underscore
chr9_10S_pat$CHR <- paste(str_split_i(chr9_10S_pat$marker, "_",1),"_",str_split_i(chr9_10S_pat$marker, "_",2), sep = "") # need to use a 3 because Chr9_10S has an underscore

#chr7S_mat_recomb <- data.frame()

chr7S_recomb <- data.frame(Positions=integer(),
                 Parent=character(),
                 stringsAsFactors=F)

buffer <- 500000
switch=0



# Figure out where recombination occurred during oogenesis
for(i in 1:(nrow(chr7S_mat)-1)) {       # for-loop over rows
  if(chr7S_mat[i,"P1_H1"] != chr7S_mat[i+1,"P1_H1"]){ # recombination occurred here
    # populate a vector with the locations of maternal recombination events
    print(paste(chr7S_mat[i,"P1_H1"]," ",chr7S_mat[i+1,"P1_H1"]," ",mean(chr7S_mat[i,"coord"],chr7S_mat[i+1,"coord"]),sep=""))
    # test whether this marker conflicts with the next marker and if so whether it it nearby
    if(switch == 0){
      chr7S_recomb[(nrow(chr7S_recomb) + 1),"Positions"] <- mean(chr7S_mat[i,"coord"],chr7S_mat[i+1,"coord"])
      chr7S_recomb[nrow(chr7S_recomb),"Parent"] <- "mat"
      switch=1
    }
    if((switch != 0)
       &&
       ((chr7S_mat[i+1,"P1_H1"] == chr7S_mat[i+2,"P1_H1"]) # this means that at least two consecutive markers support a recombination event
       &&
       (as.numeric(chr7S_mat[i,"coord"]) - as.numeric(chr7S_recomb[nrow(chr7S_recomb),"Positions"]) > buffer)
       )) # this means that the previous recombination event was at least $buffer before this one
      { 
      chr7S_recomb[(nrow(chr7S_recomb) + 1),"Positions"] <- mean(chr7S_mat[i,"coord"],chr7S_mat[i+1,"coord"])
      chr7S_recomb[nrow(chr7S_recomb),"Parent"] <- "mat" 
    }  
  }
}
#View(chr7S_recomb)

switch=0
# Figure out where recombination occurred during oogenesis
for(i in 1:(nrow(chr7S_pat)-1)) {       # for-loop over rows
  if(chr7S_pat[i,"P1_H1"] != chr7S_pat[i+1,"P1_H1"]){ # recombination occurred here
    # populate a vector with the locations of paternal recombination events
    print(paste(chr7S_pat[i,"P1_H1"]," ",chr7S_pat[i+1,"P1_H1"]," ",mean(chr7S_pat[i,"coord"],chr7S_pat[i+1,"coord"]),sep=""))
    # test whether this marker conflicts with the next marker and if so whether it it nearby
    if(switch == 0){
      chr7S_recomb[(nrow(chr7S_recomb) + 1),"Positions"] <- mean(chr7S_pat[i,"coord"],chr7S_pat[i+1,"coord"])
      chr7S_recomb[nrow(chr7S_recomb),"Parent"] <- "pat"
      switch=1
    }
    if((switch != 0)
       &&
       ((chr7S_pat[i+1,"P1_H1"] == chr7S_pat[i+2,"P1_H1"]) # this means that at least two consecutive markers support a recombination event
        &&
        (as.numeric(chr7S_pat[i,"coord"]) - as.numeric(chr7S_recomb[nrow(chr7S_recomb),"Positions"]) > buffer)
       )) # this means that the previous recombination event was at least $buffer before this one
    { 
      chr7S_recomb[(nrow(chr7S_recomb) + 1),"Positions"] <- mean(chr7S_pat[i,"coord"],chr7S_pat[i+1,"coord"])
      chr7S_recomb[nrow(chr7S_recomb),"Parent"] <- "pat" 
    }  
  }
}




ggplot(chr7S_recomb, aes(x = Positions)) +
  geom_density(aes(color = Parent))+
  geom_vline(xintercept=49000000)



recombination <- function(r,matpat){
  switch=0
  # Figure out where recombination occurred during oogenesis
  for(i in 1:(nrow(r)-1)) {       # for-loop over rows
    if(r[i,"P1_H1"] != r[i+1,"P1_H1"]){ # recombination occurred here
      # populate a vector with the locations of paternal recombination events
      print(paste(r[i,"P1_H1"]," ",r[i+1,"P1_H1"]," ",mean(r[i,"coord"],r[i+1,"coord"]),sep=""))
      # test whether this marker conflicts with the next marker and if so whether it it nearby
      if(switch == 0){
        chr7S_recomb[(nrow(chr7S_recomb) + 1),"Positions"] <- mean(r[i,"coord"],r[i+1,"coord"])
        chr7S_recomb[nrow(chr7S_recomb),"Parent"] <- matpat
        switch=1
      }
      if((switch != 0)
         &&
         ((r[i+1,"P1_H1"] == r[i+2,"P1_H1"]) # this means that at least two consecutive markers support a recombination event
          &&
          (as.numeric(r[i,"coord"]) - as.numeric(chr7S_recomb[nrow(chr7S_recomb),"Positions"]) > buffer)
         )) # this means that the previous recombination event was at least $buffer before this one
      { 
        chr7S_recomb[(nrow(chr7S_recomb) + 1),"Positions"] <- mean(r[i,"coord"],r[i+1,"coord"])
        chr7S_recomb[nrow(chr7S_recomb),"Parent"] <- matpat 
      }  
    }
  }
}
recombination(chr7S_pat,"pat")
```
# Concatenate all files but save the header
```sh
awk 'FNR==1 && NR!=1 { while (/^<header>/) getline; } 1 {print}' 2023*.txt >all.txt
```

