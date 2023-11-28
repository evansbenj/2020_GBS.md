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
If the OneMap analysis failed I (very) lightly thinned the vcf file and tried again like this:
```
vcftools --vcf GE_Chr1_removed.vcf --thin 20 --recode --out GE_Chr1_removed_thinned
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


library(onemap)
library(ggplot2)
library(tidyverse)

rm(list=ls()) # removes all variables
# setwd("/Users/Shared/Previously Relocated Items/Security/projects/2022_GBS_lotsof_Xennies/2023_OneMap/2023_trop_GOOD/GE_family")

# load the data ----
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
                           #parent1 = c("./3897mom_trim_sorted"), 
                           parent1 = mom_C,
                           #parent2 = c("./3896dad_trim_sorted"), 
                           parent2 = dad_C,
                           only_biallelic = TRUE,
                           verbose = TRUE); my_dat


# filter sites with lots of missing data ----
# higher threshold means more data is kept
data_filtered <- filter_missing(my_dat, threshold = 0.8);data_filtered


# do not try to find redundant markers; this changes the order of markers ----
# exact = FALSE means missing data will not be considered
# as recommended here: https://statgen-esalq.github.io/tutorials/onemap/Outcrossing_Populations.html#find-redundant-markers
#bins <- find_bins(data_filtered, exact = FALSE);bins

# Create a new onemap object without redundant markers ----
# this object does not have markers in order of coordinates
#new_onemap_object_no_redundant <- create_data_bins(my_dat, bins); new_onemap_object_no_redundant # This is an object of class 'onemap'

# plot the data
#plot(new_onemap_object_no_redundant)

# If you want to create a new input file with the dataset you are working on 
# after using these functions, you can use the function write_onemap_raw.
# write_onemap_raw(new_onemap_object_no_redundant, file.name = "Chr2_noredundant_onemapobj.raw")


# identify markers with or without segregation distortion ----
segreg_test <- test_segregation(data_filtered)
dist <- select_segreg(segreg_test, distorted = TRUE, numbers = TRUE) #to show the markers numbers with segregation distortion
dist
no_dist <- select_segreg(segreg_test, distorted = FALSE, numbers = TRUE) #to show the markers numbers without segregation distortion
no_dist

# Estimate the recombination fraction between all pairs of markers using two-point tests ----
# rm_mks = T gets rid of ones with weird segregation
# due to excess of missing data or segregation deviation
# but this doesn't seem to really work, so I did it again in the next step
twopts <- rf_2pts(data_filtered, rm_mks = T); twopts # an object of class rf_2pts

# make a variable with the markers that do not have segregation distortion ----
mark_no_dist <- make_seq(twopts, c(no_dist)) # this is an object of class 'sequence'


# look at markers
#marker_type(mark_no_dist)
p <- marker_type(mark_no_dist)
#p %>% count(p$Type)

# get indexes of mat and pat sites ----
# the positions should be sorted now
#matpat_SNPs <- as.integer(p_sorted[(p_sorted$Type == 'B3.7')|(p_sorted$Type == 'D1.10')|(p_sorted$Type == 'D2.15'), ]$Marker)
matpat_SNPs <- as.integer(p[(p$Type == 'B3.7')|(p$Type == 'D1.10')|(p$Type == 'D2.15'), ]$Marker)


# it seems that if the B3.7 type is not included, that no LG can be formed from 
# only 'D1.10' or only 'D2.15'
# but when I include B3.7 sites, the largest mat and pat LGs appears 
# to only have B3.7 type markers...!
#mat_SNPs <- as.integer(p[(p$Type == 'B3.7')|(p$Type == 'D1.10'), ]$Marker)
#pat_SNPs <- as.integer(p[(p$Type == 'B3.7')|(p$Type == 'D2.15'), ]$Marker)

# I'm not including the informative doublehet sites because of the large LGs had only this
# type of site in them and prevented maternal and paternal LGs from then being inferred

# Marker types are explained here:
# https://statgen-esalq.github.io/tutorials/onemap/Outcrossing_Populations.html

# 1 B3.7    321  ab×ab (double het)
# 2 D1.10   483  ab×aa (first parent het; maternal)
# 3 D2.15   384  aa×ab (second parent het; paternal)

# make an object with only matpat SNPs ----
matpat_SNPs_no_dist <- make_seq(twopts, unique(matpat_SNPs)) # these are objects of class 'sequence'
#mat_SNPs_no_dist <- make_seq(twopts, mat_SNPs) # these are objects of class 'sequence'
#pat_SNPs_no_dist <- make_seq(twopts, pat_SNPs) # these are objects of class 'sequence'

# look at markers
#marker_type(matpat_SNPs_no_dist)
p <- marker_type(matpat_SNPs_no_dist)
p %>% count(p$Type)


# Don't use the function suggest_lod to calculate a suggested LOD score 
# considering that multiple tests are being performed. If you do, you will get very
# samll linkage groups that often don't contain maternal and paternal markers
# LOD_sug_matpat <- suggest_lod(matpat_SNPs_no_dist)

# Group the two-point tests and set the maximum recombination fraction to 0.40:
# Furman et al. 2020 used LOD=5; here I use LOD=3 to be more inclusive with markers
matpat_LGs <- group(matpat_SNPs_no_dist, LOD=3, max.rf = 0.4)

# figure out which linkage group has the most markers
matpat_df <- as.data.frame(table(matpat_LGs$groups))

# get rid of first row, which is the unlinked markers
matpat_df <- matpat_df[-1, ]
matpat_biggest_LG_group <- as.vector(matpat_df$Var1[matpat_df$Freq == max(matpat_df$Freq)]);matpat_biggest_LG_group
print(paste("The biggest matpat LG is",matpat_biggest_LG_group,sep=" "))



# set_map_fun(type = "kosambi") # this is the default
# extract the largest linkage group for mapping ----
# the order of the markers in this object does not match the coordinates
matpat_biggest_LG <- make_seq(matpat_LGs,as.integer(matpat_biggest_LG_group))
rf_graph_table(matpat_biggest_LG)
matpat_LG_1_graph <- rf_graph_table(matpat_biggest_LG);matpat_LG_1_graph
markers_to_keep <- matpat_biggest_LG$seq.num

# marker_type(matpat_biggest_LG)
p <- marker_type(matpat_biggest_LG)
p %>% count(p$Type)

# reload the data ----
# now reload the data and save only the markers in the largest LG group we have identified
# then we can hopefully force coordinates without trying to estimate the LG


# filter sites with lots of missing data ----
data_filtered2 <- filter_missing(my_dat, threshold = 0.8);data_filtered2

# no not find redundant markers ----
# exact = FALSE means missing data will not be considered
# as recommended here: https://statgen-esalq.github.io/tutorials/onemap/Outcrossing_Populations.html#find-redundant-markers
#bins2 <- find_bins(data_filtered2, exact = FALSE)
#bins2

# Create a new onemap object without redundant markers ----
# this step messes with the marker order
#new_onemap_object_no_redundant2 <- create_data_bins(my_dat2, bins2)
#new_onemap_object_no_redundant2 # This is an object of class 'onemap'

# Estimate the recombination fraction between all pairs of markers using two-point tests ----
# rm_mks = T gets rid of ones with weird segregation
twopts2 <- rf_2pts(data_filtered2, rm_mks = T); twopts2 # an object of class rf_2pts

# make a sequence of markers from the chromosome of interest
# this has all markers in it, including the ones that don't form a LG


new_biggest_LG <- make_seq(twopts2, prefix_C) # new_biggest_LG is class 'sequence'

# Make a vector with a list of markers to remove
markers_to_remove <- setdiff(new_biggest_LG$seq.num, markers_to_keep)

# Now drop all markers except markers_to_keep, which are the ones that form a LG
matpat_biggest_LG_coord <- drop_marker(matpat_biggest_LG, markers_to_remove)

# Check that everything looks good
rf_graph_table(matpat_biggest_LG_coord) 

#marker_type(matpat_biggest_LG_coord)
p <- marker_type(matpat_biggest_LG_coord)
p %>% count(p$Type)

# print the recombination map
ggsave(file=paste(prefix_C,"_matpat_graph.pdf",sep="_"), matpat_LG_1_graph, width=38, height=8)


# At this point no parameters have been estimated ----
# Make two vectors that have the markers we want to remove for the mat and pat maps
not_mat_D1.10 <- marker_type(matpat_biggest_LG_coord)[marker_type(matpat_biggest_LG_coord)$Type != "D1.10", ]$Marker
not_pat_D2.15 <- marker_type(matpat_biggest_LG_coord)[marker_type(matpat_biggest_LG_coord)$Type != "D2.15", ]$Marker

# drop non-mat or non-pat markers and estimate parameters
matonly_biggest_LG <- drop_marker(matpat_biggest_LG_coord, not_mat_D1.10)
patonly_biggest_LG <- drop_marker(matpat_biggest_LG_coord, not_pat_D2.15)

# how many markers are there in the mat and pat LGs?
length(matonly_biggest_LG$seq.num)
length(patonly_biggest_LG$seq.num)

# Check that everything looks good
#rf_graph_table(matonly_biggest_LG) 
#rf_graph_table(patonly_biggest_LG) 


# this step now estimates the multipoint log-likelihood, 
# linkage phases and recombination frequencies for 
# a sequence of markers in a given order.
print("Making the mat map ")
mat_map <- onemap::map(matonly_biggest_LG)
print("Making the pat map ")
pat_map <- onemap::map(patonly_biggest_LG)



#mat_LG_graph <- rf_graph_table(mat_map)
#mat_LG_graph_LOD <- rf_graph_table(mat_map, graph.LOD = TRUE)
#ggsave(file=paste(prefix_C,"_mat_graph.pdf",sep="_"), mat_LG_1_graph, width=8, height=8)
#ggsave(file=paste(prefix_C,"_mat_graph_LOD.pdf",sep="_"), mat_LG_1_graph_LOD, width=8, height=8)

#pat_LG_graph <- rf_graph_table(pat_map)
#pat_LG_graph_LOD <- rf_graph_table(pat_map, graph.LOD = TRUE)
#ggsave(file=paste(prefix_C,"_pat_graph.pdf",sep="_"), pat_LG_1_graph, width=8, height=8)
#ggsave(file=paste(prefix_C,"_pat_graph_LOD.pdf",sep="_"), pat_LG_1_graph_LOD, width=8, height=8)

maternal_parents_haplot <- parents_haplotypes(mat_map)
paternal_parents_haplot <- parents_haplotypes(pat_map)
# matpat_haplot <- parents_haplotypes(matpat_map)

# Add column for matpat
maternal_parents_haplot$matpat <- "mat"
paternal_parents_haplot$matpat <- "pat"

#prefix_C<-"temp"
write.table(maternal_parents_haplot, paste(prefix_C,"_mat_parents_haplot_LG.txt",sep="_"), row.names=F)
write.table(paternal_parents_haplot, paste(prefix_C,"_pat_parents_haplot_LG.txt",sep="_"), row.names=F)

# draw_map(maternal_map)

# Export haplotypes
mat_progeny_haplot <- progeny_haplotypes(mat_map, 
                                     most_likely = TRUE, 
                                     ind = c(1:data_filtered$n.ind), 
                                     group_names = prefix_C)

pat_progeny_haplot <- progeny_haplotypes(pat_map, 
                                         most_likely = TRUE, 
                                         ind = c(1:data_filtered$n.ind), 
                                         group_names = prefix_C)


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

write.table(mat_progeny_haplot_wide, paste(prefix_C,"_mat_progeny_haplot_wide_LG.txt",sep="_"),row.names = F)
write.table(pat_progeny_haplot_wide, paste(prefix_C,"_pat_progeny_haplot_wide_LG.txt",sep="_"),row.names = F)

# how many markers are there in the mat and pat LGs?
length(matonly_biggest_LG$seq.num)
length(patonly_biggest_LG$seq.num)


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
# Concatenate all files but save the header
```sh
head -1 allo_1_chr1S_mat_progeny_haplot_wide.txt > all_mat.txt; awk 'FNR>1{print}' *mat*wide_LG.txt >> all_mat.txt
```

# Processing script for OneMap output
```R
library(tidyverse)

# concatenate "wide" haplotypes from all individuals for maternal and paternal recombination events like this:
#  head -1 one_of_the_files_haplot_wide.txt > all_mat.txt; awk 'FNR>1{print}' *mat*wide*.txt >> all_mat.txt
#  head -1 one_of_the_files_pat_progeny_haplot_wide.txt > all_pat.txt; awk 'FNR>1{print}' *pat*wide*.txt >> all_pat.txt
# if there are rownames, you can add a bogus header for them to prevent a redundant rowname error like this:
# perl -pi -e 's/^/"bogus" / if $.==1' all_mat.txt
options(scipen=999)
#setwd("/Users/Shared/Previously Relocated Items/Security/projects/2022_GBS_lotsof_Xennies/2023_OneMap/2023_trop_GOOD/GE_family")
setwd('/Users/Shared/Previously Relocated Items/Security/projects/2022_GBS_lotsof_Xennies/2023_OneMap/2023_trop_GOOD/GW_family')
setwd('/Users/Shared/Previously Relocated Items/Security/projects/2022_GBS_lotsof_Xennies/2023_OneMap/2023_laevis_good')
all_mat <- read.table("./all_mat.txt", header = T)
all_pat <- read.table("./all_pat.txt", header = T)

# ad matpat column
all_mat$matpat <- "mat"
all_pat$matpat <- "pat"

all<- rbind(all_mat,all_pat)

# make a column with chrs and one with coordinates
all[c('Chr', 'Coord')] <- str_split_fixed(all$marker, '_', 2)

# the two lines below are only for XL
all[all$Chr == 'Chr9',]$Coord <- as.numeric(str_split_i(all[all$Chr == 'Chr9',]$marker,"_",3)) # need to use a 3 because Chr9_10L has an underscore
all[all$Chr == 'Chr9',]$Chr <- paste(str_split_i(all[all$Chr == 'Chr9',]$marker, "_",1),str_split_i(all[all$Chr == 'Chr9',]$marker, "_",2), sep = "_") # need to use a 3 because Chr9_10L has an underscore

# check the chrs
unique(all$Chr)

all$Coord <- as.numeric(all$Coord)

# now make a column with chromosome proportions for XL
all$Proportions <- all$Coord/233740090 # this is the length of XL Chr1L
# now update the ones that are not Chr1L
all[all$Chr == 'Chr1S',]$Proportions <- all[all$Chr == 'Chr1S',]$Coord/202412970 # this is the length of XL Chr1S
all[all$Chr == 'Chr2L',]$Proportions <- all[all$Chr == 'Chr2L',]$Coord/191000146 # this is the length of XL Chr2L
all[all$Chr == 'Chr2S',]$Proportions <- all[all$Chr == 'Chr2S',]$Coord/169306100 # this is the length of XL Chr2S
all[all$Chr == 'Chr3L',]$Proportions <- all[all$Chr == 'Chr3L',]$Coord/161426101 # this is the length of XL Chr3L
all[all$Chr == 'Chr3S',]$Proportions <- all[all$Chr == 'Chr3S',]$Coord/131962816 # this is the length of XL Chr3S
all[all$Chr == 'Chr4L',]$Proportions <- all[all$Chr == 'Chr4L',]$Coord/155250554 # this is the length of XL Chr4L
all[all$Chr == 'Chr4S',]$Proportions <- all[all$Chr == 'Chr4S',]$Coord/132731174 # this is the length of XL Chr4S
all[all$Chr == 'Chr5L',]$Proportions <- all[all$Chr == 'Chr5L',]$Coord/171415384 # this is the length of XL Chr5L
all[all$Chr == 'Chr5S',]$Proportions <- all[all$Chr == 'Chr5S',]$Coord/143394103 # this is the length of XL Chr5S
all[all$Chr == 'Chr6L',]$Proportions <- all[all$Chr == 'Chr6L',]$Coord/164223595 # this is the length of XL Chr6L
all[all$Chr == 'Chr6S',]$Proportions <- all[all$Chr == 'Chr6S',]$Coord/137316286 # this is the length of XL Chr6S
all[all$Chr == 'Chr7L',]$Proportions <- all[all$Chr == 'Chr7L',]$Coord/139837618 # this is the length of XL Chr7L
all[all$Chr == 'Chr7S',]$Proportions <- all[all$Chr == 'Chr7S',]$Coord/113060389 # this is the length of XL Chr7S
all[all$Chr == 'Chr8L',]$Proportions <- all[all$Chr == 'Chr8L',]$Coord/135449133 # this is the length of XL Chr8L
all[all$Chr == 'Chr8S',]$Proportions <- all[all$Chr == 'Chr8S',]$Coord/103977862 # this is the length of XL Chr8S
all[all$Chr == 'Chr9_10L',]$Proportions <- all[all$Chr == 'Chr9_10L',]$Coord/137811819 # this is the length of XL Chr9_10L
all[all$Chr == 'Chr9_10S',]$Proportions <- all[all$Chr == 'Chr9_10S',]$Coord/117266291 # this is the length of XL Chr9_10S

# or for XT
all$Proportions <- all$Coord/217471165 # this is the length of XT Chr1
# now update the ones that are not Chr1
all[all$Chr == 'Chr2',]$Proportions <- all[all$Chr == 'Chr2',]$Coord/181034960 # this is the length of XT Chr2
all[all$Chr == 'Chr3',]$Proportions <- all[all$Chr == 'Chr3',]$Coord/153873356 # this is the length of XT Chr3
all[all$Chr == 'Chr4',]$Proportions <- all[all$Chr == 'Chr4',]$Coord/153961318 # this is the length of XT Chr4
all[all$Chr == 'Chr5',]$Proportions <- all[all$Chr == 'Chr5',]$Coord/164033574 # this is the length of XT Chr5
all[all$Chr == 'Chr6',]$Proportions <- all[all$Chr == 'Chr6',]$Coord/154486311 # this is the length of XT Chr6
all[all$Chr == 'Chr7',]$Proportions <- all[all$Chr == 'Chr7',]$Coord/133565929 # this is the length of XT Chr7
all[all$Chr == 'Chr8',]$Proportions <- all[all$Chr == 'Chr8',]$Coord/147241509 # this is the length of XT Chr8
all[all$Chr == 'Chr9',]$Proportions <- all[all$Chr == 'Chr9',]$Coord/91218943 # this is the length of XT Chr9
all[all$Chr == 'Chr10',]$Proportions <- all[all$Chr == 'Chr10',]$Coord/52432565 # this is the length of XT Chr10





max(all$Proportions) # should < 1
min(all$Proportions) # should be > 0
# for Xl
# Chr1L	233740090
# Chr1S	202412970
# Chr2L	191000146
# Chr2S	169306100
# Chr3L	161426101
# Chr3S	131962816
# Chr4L	155250554
# Chr4S	132731174
# Chr5L	171415384
# Chr5S	143394103
# Chr6L	164223595
# Chr6S	137316286
# Chr7L	139837618
# Chr7S	113060389
# Chr8L	135449133
# Chr8S	103977862
# Chr9_10L	137811819
# Chr9_10S	117266291



# now make a column with chromosome Proportions_relative_to_centromeres for XL
# this just makes a column with a dummy value
all$Proportions_relative_to_centromeres <- all$Coord/233740090 # this is the length of XL Chr1L
# now update for each chr for XL
all[(all$Chr == 'Chr1L')&(all$Coord < 97110544),]$Proportions_relative_to_centromeres <- all[(all$Chr == 'Chr1L')&(all$Coord < 97110544),]$Coord/97110544 # this is the first part of XL Chr1L
all[(all$Chr == 'Chr1L')&(all$Coord >= 97110544),]$Proportions_relative_to_centromeres <- (233740090 - all[(all$Chr == 'Chr1L')&(all$Coord >= 97110544),]$Coord)/(233740090-97110544) # this is the first part of XL Chr1L
all[(all$Chr == 'Chr1S')&(all$Coord < 78965477),]$Proportions_relative_to_centromeres <- all[(all$Chr == 'Chr1S')&(all$Coord < 78965477),]$Coord/78965477 # this is the first part of XL Chr1S
all[(all$Chr == 'Chr1S')&(all$Coord >= 78965477),]$Proportions_relative_to_centromeres <- (202412970 - all[(all$Chr == 'Chr1S')&(all$Coord >= 78965477),]$Coord)/(202412970-78965477) # this is the first part of XL Chr1L

all[(all$Chr == 'Chr2L')&(all$Coord < 70167544),]$Proportions_relative_to_centromeres <- all[(all$Chr == 'Chr2L')&(all$Coord < 70167544),]$Coord/70167544 # this is the first part of XL Chr2L
all[(all$Chr == 'Chr2L')&(all$Coord >= 70167544),]$Proportions_relative_to_centromeres <- (191000146 - all[(all$Chr == 'Chr2L')&(all$Coord >= 70167544),]$Coord)/(191000146-70167544) # this is the first part of XL Chr2L
all[(all$Chr == 'Chr2S')&(all$Coord < 54769165),]$Proportions_relative_to_centromeres <- all[(all$Chr == 'Chr2S')&(all$Coord < 54769165),]$Coord/54769165 # this is the first part of XL Chr2S
all[(all$Chr == 'Chr2S')&(all$Coord >= 54769165),]$Proportions_relative_to_centromeres <- (169306100 - all[(all$Chr == 'Chr2S')&(all$Coord >= 54769165),]$Coord)/(169306100-54769165) # this is the first part of XL Chr2L

all[(all$Chr == 'Chr3L')&(all$Coord < 19199725),]$Proportions_relative_to_centromeres <- all[(all$Chr == 'Chr3L')&(all$Coord < 19199725),]$Coord/19199725 # this is the first part of XL Chr3L
all[(all$Chr == 'Chr3L')&(all$Coord >= 19199725),]$Proportions_relative_to_centromeres <- (161426101 - all[(all$Chr == 'Chr3L')&(all$Coord >= 19199725),]$Coord)/(161426101-19199725) # this is the first part of XL Chr3L
all[(all$Chr == 'Chr3S')&(all$Coord < 19205290),]$Proportions_relative_to_centromeres <- all[(all$Chr == 'Chr3S')&(all$Coord < 19205290),]$Coord/19205290 # this is the first part of XL Chr3S
all[(all$Chr == 'Chr3S')&(all$Coord >= 19205290),]$Proportions_relative_to_centromeres <- (131962816 - all[(all$Chr == 'Chr3S')&(all$Coord >= 19205290),]$Coord)/(131962816-19205290) # this is the first part of XL Chr3L

all[(all$Chr == 'Chr4L')&(all$Coord < 36455216),]$Proportions_relative_to_centromeres <- all[(all$Chr == 'Chr4L')&(all$Coord < 36455216),]$Coord/36455216 # this is the first part of XL Chr4L
all[(all$Chr == 'Chr4L')&(all$Coord >= 36455216),]$Proportions_relative_to_centromeres <- (155250554 - all[(all$Chr == 'Chr4L')&(all$Coord >= 36455216),]$Coord)/(155250554-36455216) # this is the first part of XL Chr4L
all[(all$Chr == 'Chr4S')&(all$Coord < 29375848),]$Proportions_relative_to_centromeres <- all[(all$Chr == 'Chr4S')&(all$Coord < 29375848),]$Coord/29375848 # this is the first part of XL Chr4S
all[(all$Chr == 'Chr4S')&(all$Coord >= 29375848),]$Proportions_relative_to_centromeres <- (132731174 - all[(all$Chr == 'Chr4S')&(all$Coord >= 29375848),]$Coord)/(132731174-29375848) # this is the first part of XL Chr4L

all[(all$Chr == 'Chr5L')&(all$Coord < 63286813),]$Proportions_relative_to_centromeres <- all[(all$Chr == 'Chr5L')&(all$Coord < 63286813),]$Coord/63286813 # this is the first part of XL Chr5L
all[(all$Chr == 'Chr5L')&(all$Coord >= 63286813),]$Proportions_relative_to_centromeres <- (171415384 - all[(all$Chr == 'Chr5L')&(all$Coord >= 63286813),]$Coord)/(171415384-63286813) # this is the first part of XL Chr5L
all[(all$Chr == 'Chr5S')&(all$Coord < 54131023),]$Proportions_relative_to_centromeres <- all[(all$Chr == 'Chr5S')&(all$Coord < 54131023),]$Coord/54131023 # this is the first part of XL Chr5S
all[(all$Chr == 'Chr5S')&(all$Coord >= 54131023),]$Proportions_relative_to_centromeres <- (143394103 - all[(all$Chr == 'Chr5S')&(all$Coord >= 54131023),]$Coord)/(143394103-54131023) # this is the first part of XL Chr5L

all[(all$Chr == 'Chr6L')&(all$Coord < 78606137),]$Proportions_relative_to_centromeres <- all[(all$Chr == 'Chr6L')&(all$Coord < 78606137),]$Coord/78606137 # this is the first part of XL Chr6L
all[(all$Chr == 'Chr6L')&(all$Coord >= 78606137),]$Proportions_relative_to_centromeres <- (164223595 - all[(all$Chr == 'Chr6L')&(all$Coord >= 78606137),]$Coord)/(164223595-78606137) # this is the first part of XL Chr6L
all[(all$Chr == 'Chr6S')&(all$Coord < 55279263),]$Proportions_relative_to_centromeres <- all[(all$Chr == 'Chr6S')&(all$Coord < 55279263),]$Coord/55279263 # this is the first part of XL Chr6S
all[(all$Chr == 'Chr6S')&(all$Coord >= 55279263),]$Proportions_relative_to_centromeres <- (137316286 - all[(all$Chr == 'Chr6S')&(all$Coord >= 55279263),]$Coord)/(137316286-55279263) # this is the first part of XL Chr6L

all[(all$Chr == 'Chr7L')&(all$Coord < 58257035),]$Proportions_relative_to_centromeres <- all[(all$Chr == 'Chr7L')&(all$Coord < 58257035),]$Coord/58257035 # this is the first part of XL Chr7L
all[(all$Chr == 'Chr7L')&(all$Coord >= 58257035),]$Proportions_relative_to_centromeres <- (139837618 - all[(all$Chr == 'Chr7L')&(all$Coord >= 58257035),]$Coord)/(139837618-58257035) # this is the first part of XL Chr7L
all[(all$Chr == 'Chr7S')&(all$Coord < 47593023),]$Proportions_relative_to_centromeres <- all[(all$Chr == 'Chr7S')&(all$Coord < 47593023),]$Coord/47593023 # this is the first part of XL Chr7S
all[(all$Chr == 'Chr7S')&(all$Coord >= 47593023),]$Proportions_relative_to_centromeres <- (113060389 - all[(all$Chr == 'Chr7S')&(all$Coord >= 47593023),]$Coord)/(113060389-47593023) # this is the first part of XL Chr7L

all[(all$Chr == 'Chr8L')&(all$Coord < 21985291),]$Proportions_relative_to_centromeres <- all[(all$Chr == 'Chr8L')&(all$Coord < 21985291),]$Coord/21985291 # this is the first part of XL Chr8L
all[(all$Chr == 'Chr8L')&(all$Coord >= 21985291),]$Proportions_relative_to_centromeres <- (135449133 - all[(all$Chr == 'Chr8L')&(all$Coord >= 21985291),]$Coord)/(135449133-21985291) # this is the first part of XL Chr8L
all[(all$Chr == 'Chr8S')&(all$Coord < 49209017),]$Proportions_relative_to_centromeres <- all[(all$Chr == 'Chr8S')&(all$Coord < 49209017),]$Coord/49209017 # this is the first part of XL Chr8S
all[(all$Chr == 'Chr8S')&(all$Coord >= 49209017),]$Proportions_relative_to_centromeres <- (103977862 - all[(all$Chr == 'Chr8S')&(all$Coord >= 49209017),]$Coord)/(103977862-49209017) # this is the first part of XL Chr8L


all[(all$Chr == 'Chr9_10L')&(all$Coord < 23877437),]$Proportions_relative_to_centromeres <- all[(all$Chr == 'Chr9_10L')&(all$Coord < 23877437),]$Coord/23877437 # this is the first part of XL Chr9_10L
all[(all$Chr == 'Chr9_10L')&(all$Coord >= 23877437),]$Proportions_relative_to_centromeres <- (137811819 - all[(all$Chr == 'Chr9_10L')&(all$Coord >= 23877437),]$Coord)/(137811819-23877437) # this is the first part of XL Chr9_10L
all[(all$Chr == 'Chr9_10S')&(all$Coord < 25122470),]$Proportions_relative_to_centromeres <- all[(all$Chr == 'Chr9_10S')&(all$Coord < 25122470),]$Coord/25122470 # this is the first part of XL Chr9_10S
all[(all$Chr == 'Chr9_10S')&(all$Coord >= 25122470),]$Proportions_relative_to_centromeres <- (117266291 - all[(all$Chr == 'Chr9_10S')&(all$Coord >= 25122470),]$Coord)/(117266291-25122470) # this is the first part of XL Chr9_10L

# now make a column with chromosome Proportions_relative_to_centromeres for XT
# for now I have just used the midpoint - so this is not correct and should not be used for XT

# this just makes a column with a dummy value
all$Proportions_relative_to_centromeres <- all$Coord/217471165 # this is the length of XT Chr1
# now update for each chr for XT
all[(all$Chr == 'Chr1')&(all$Coord < 108735582.5),]$Proportions_relative_to_centromeres <- all[(all$Chr == 'Chr1')&(all$Coord < 108735582.5),]$Coord/108735582.5 # this is the first part of XT Chr1
all[(all$Chr == 'Chr1')&(all$Coord >= 108735582.5),]$Proportions_relative_to_centromeres <- (217471165 - all[(all$Chr == 'Chr1')&(all$Coord >= 108735582.5),]$Coord)/(217471165-108735582.5) # this is the first part of XT Chr1

all[(all$Chr == 'Chr2')&(all$Coord < 90517480),]$Proportions_relative_to_centromeres <- all[(all$Chr == 'Chr2')&(all$Coord < 90517480),]$Coord/90517480 # this is the first part of XT Chr2
all[(all$Chr == 'Chr2')&(all$Coord >= 90517480),]$Proportions_relative_to_centromeres <- (181034960 - all[(all$Chr == 'Chr2')&(all$Coord >= 90517480),]$Coord)/(181034960-90517480) # this is the first part of XT Chr2

all[(all$Chr == 'Chr3')&(all$Coord < 76936678),]$Proportions_relative_to_centromeres <- all[(all$Chr == 'Chr3')&(all$Coord < 76936678),]$Coord/76936678 # this is the first part of XT Chr3
all[(all$Chr == 'Chr3')&(all$Coord >= 76936678),]$Proportions_relative_to_centromeres <- (153873356 - all[(all$Chr == 'Chr3')&(all$Coord >= 76936678),]$Coord)/(153873356-76936678) # this is the first part of XT Chr3

all[(all$Chr == 'Chr4')&(all$Coord < 76980659),]$Proportions_relative_to_centromeres <- all[(all$Chr == 'Chr4')&(all$Coord < 76980659),]$Coord/76980659 # this is the first part of XT Chr4
all[(all$Chr == 'Chr4')&(all$Coord >= 76980659),]$Proportions_relative_to_centromeres <- (153961318 - all[(all$Chr == 'Chr4')&(all$Coord >= 76980659),]$Coord)/(153961318-76980659) # this is the first part of XT Chr4

all[(all$Chr == 'Chr5')&(all$Coord < 82016787),]$Proportions_relative_to_centromeres <- all[(all$Chr == 'Chr5')&(all$Coord < 82016787),]$Coord/82016787 # this is the first part of XT Chr5
all[(all$Chr == 'Chr5')&(all$Coord >= 82016787),]$Proportions_relative_to_centromeres <- (164033574 - all[(all$Chr == 'Chr5')&(all$Coord >= 82016787),]$Coord)/(164033574-82016787) # this is the first part of XT Chr5

all[(all$Chr == 'Chr6')&(all$Coord < 77243155.5),]$Proportions_relative_to_centromeres <- all[(all$Chr == 'Chr6')&(all$Coord < 77243155.5),]$Coord/77243155.5 # this is the first part of XT Chr6
all[(all$Chr == 'Chr6')&(all$Coord >= 77243155.5),]$Proportions_relative_to_centromeres <- (154486311 - all[(all$Chr == 'Chr6')&(all$Coord >= 77243155.5),]$Coord)/(154486311-77243155.5) # this is the first part of XT Chr6

all[(all$Chr == 'Chr7')&(all$Coord < 66782964.5),]$Proportions_relative_to_centromeres <- all[(all$Chr == 'Chr7')&(all$Coord < 66782964.5),]$Coord/66782964.5 # this is the first part of XT Chr7
all[(all$Chr == 'Chr7')&(all$Coord >= 66782964.5),]$Proportions_relative_to_centromeres <- (133565929 - all[(all$Chr == 'Chr7')&(all$Coord >= 66782964.5),]$Coord)/(133565929-66782964.5) # this is the first part of XT Chr7

all[(all$Chr == 'Chr8')&(all$Coord < 73620754.5),]$Proportions_relative_to_centromeres <- all[(all$Chr == 'Chr8')&(all$Coord < 73620754.5),]$Coord/73620754.5 # this is the first part of XT Chr8
all[(all$Chr == 'Chr8')&(all$Coord >= 73620754.5),]$Proportions_relative_to_centromeres <- (147241509 - all[(all$Chr == 'Chr8')&(all$Coord >= 73620754.5),]$Coord)/(147241509-73620754.5) # this is the first part of XT Chr8

all[(all$Chr == 'Chr9')&(all$Coord < 45609471.5),]$Proportions_relative_to_centromeres <- all[(all$Chr == 'Chr9')&(all$Coord < 45609471.5),]$Coord/45609471.5 # this is the first part of XT Chr9
all[(all$Chr == 'Chr9')&(all$Coord >= 45609471.5),]$Proportions_relative_to_centromeres <- (91218943 - all[(all$Chr == 'Chr9')&(all$Coord >= 45609471.5),]$Coord)/(91218943-45609471.5) # this is the first part of XT Chr9

all[(all$Chr == 'Chr10')&(all$Coord < 26216282.5),]$Proportions_relative_to_centromeres <- all[(all$Chr == 'Chr10')&(all$Coord < 26216282.5),]$Coord/26216282.5 # this is the first part of XT Chr10
all[(all$Chr == 'Chr10')&(all$Coord >= 26216282.5),]$Proportions_relative_to_centromeres <- (52432565 - all[(all$Chr == 'Chr10')&(all$Coord >= 26216282.5),]$Coord)/(52432565-26216282.5) # this is the first part of XT Chr10


max(all$Proportions_relative_to_centromeres) # should < 1
min(all$Proportions_relative_to_centromeres) # should > 0

# From the Smith et al paper	GSE153058_xla_v10.2_cen 1.bed	https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153058	
# Chr   Start     Stop      Average_of_centromere_coordinates	Chr_Length
# Chr1L	96582536	97638551	97110544	233740090
# Chr1S	78922130	79008823	78965477	202412970
# Chr2L	70076389	70258699	70167544	191000146
# Chr2S	54347978	55190352	54769165	169306100
# Chr3L	19140049	19259400	19199725	161426101
# Chr3S	18944075	19466505	19205290	131962816
# Chr4L	36375516	36534915	36455216	155250554
# Chr4S	29328671	29423025	29375848	132731174
# Chr5L	63210646	63362980	63286813	171415384
# Chr5S	54099924	54162122	54131023	143394103
# Chr6L	78555074	78657199	78606137	164223595
# Chr6S	55253614	55304911	55279263	137316286
# Chr7L	58158940	58355130	58257035	139837618
# Chr7S	47582934	47603112	47593023	113060389
# Chr8L	21900593	22069988	21985291	135449133
# Chr8S	49180152	49237882	49209017	103977862
# Chr9_10L	23812947	23941927	23877437	137811819
# Chr9_10S	25065099	25179840	25122470	117266291
    
all_recomb <- data.frame(Chr=character(),
                         Positions=integer(),
                         Proportions=integer(),
                         Proportions_relative_to_centromeres=integer(),
                         Parent=character(),
                         stringsAsFactors=F)
buffer <- 10000000
previous_chr <- "temp"
previous_position <- 0
previous_individual <- "temp"

# make a df to track the number of mat and pat recombinations per individual per chromosome
# this will hopefully be useful to sanity check results
# first make a matrix with the correct dimensions
new_matrix <- matrix(rep(0,length(unique(all$ind))*length(unique(all$Chr))), ncol = length(unique(all$ind)));new_matrix
mat_recomb_per_chr <- as.data.frame(new_matrix) 
pat_recomb_per_chr <- as.data.frame(new_matrix) 
colnames(mat_recomb_per_chr) <- unique(all$ind)
row.names(mat_recomb_per_chr) <- unique(all$Chr);mat_recomb_per_chr
colnames(pat_recomb_per_chr) <- unique(all$ind)
row.names(pat_recomb_per_chr) <- unique(all$Chr);pat_recomb_per_chr

# Figure out where recombination occurred during oogenesis
for(i in 1:(nrow(all)-2)) {       # for-loop over rows; need to only go up to the second to last to allow
                                  # a check for at least two sites supporting recombination
  if((all[i,]$grp != previous_chr)|(all[i,]$ind != previous_individual)){
    previous_chr <- all[i,]$grp
    previous_position <- 0
    previous_individual <- all[i,]$ind
  }
  if((all[i,]$matpat == "mat")&
   #  (is.na(all[i,"P1_H1"]) == F)& #commented because NA should only be present in bad files
   #   (is.na(all[i+1,"P1_H1"]) == F)& # so I want this to fail if an NA is encountered
     (all[i,]$grp == all[i+1,]$grp)&
     (all[i,]$grp == all[i+2,]$grp)){ # this is for maternal SNPS on the same chr
    if((all[i,"P1_H1"] != all[i+1,"P1_H1"])&
       (as.numeric(all[i,]$pos) != 0)&
       (all[i,]$ind == previous_individual)){ # recombination occurred here in the mother
      # potentially populate a vector with the locations of maternal recombination events if two downstream sites support this
      # print(paste(all[i,"P1_H1"]," ",all[i+1,"P1_H1"]," ",all[i+1,"Coord"],sep=""))
      # test whether this marker conflicts with the next two markers and if so whether it is nearby
      if((all[i-1,"P1_H1"] == all[i,"P1_H1"]) & # these five conditions mean that at least two previous and two consecutive markers support a recombination event
         (all[i+1,"P1_H1"] == all[i+2,"P1_H1"]) &
         (all[i-1,"ind"] == all[i,"ind"]) & 
         (all[i+1,"ind"] == all[i+2,"ind"]) &
         (all[i,"ind"] == all[i+1,"ind"]) &
          (((as.numeric(all[i,"Coord"])) - as.numeric(previous_position)) > buffer)
         ) # this means that the previous recombination event was at least $buffer before this one
      { 
       # print(paste(all[i,"ind"]," ",all[i,"grp"]," ",all[i-1,"P1_H1"]," ",all[i,"P1_H1"]," ",all[i+1,"P1_H1"]," ",all[i+2,"P1_H1"]," ",
       #             all[i-1,"Coord"]," ",all[i,"Coord"]," ",all[i+1,"Coord"]," ",all[i+2,"Coord"]," ",mean(c(all[i,"Coord"],all[i+1,"Coord"]))," ",previous_position,sep=""))
        all_recomb[(nrow(all_recomb) + 1),"Positions"] <- mean(c(all[i,"Coord"],all[i+1,"Coord"]))
        all_recomb[nrow(all_recomb),"Proportions"] <- mean(c(all[i,"Proportions"],all[i+1,"Proportions"]))
        all_recomb[nrow(all_recomb),"Proportions_relative_to_centromeres"] <- mean(c(all[i,"Proportions_relative_to_centromeres"],all[i+1,"Proportions_relative_to_centromeres"]))
        all_recomb[nrow(all_recomb),"Parent"] <- all[i,"matpat"] 
        all_recomb[nrow(all_recomb),"Chr"] <- all[i,"Chr"] 
        previous_chr <- all[i,]$grp
        previous_position <- mean(c(all[i,"Coord"],all[i+1,"Coord"]))
        previous_individual <- all[i,"ind"]
        mat_recomb_per_chr[all[i,]$Chr,all[i,"ind"]] <-  mat_recomb_per_chr[all[i,]$Chr,all[i,"ind"]]+1
      }  
    }
  } # end of check for maternal recomb
  else if((all[i,]$matpat == "pat")&
        #  (is.na(all[i,"P2_H1"]) == F)& #commented because NA should only be present in bad files
        #  (is.na(all[i+1,"P2_H1"]) == F)& # so I want this to fail if an NA is encountered
          (all[i,]$grp == all[i+1,]$grp)&
          (all[i,]$grp == all[i+2,]$grp)){ # this is for paternal SNPS on the same chr
    print(paste(i," ",all[i+3,"ind"]," ",all[i+3,"grp"]," ",all[i+3,]$pos," ",all[i-1,"P2_H1"]," ",all[i,"P2_H1"]," ",all[i+1,"P2_H1"]," ",
                all[i+2,"P2_H1"]," ", all[i-1,"Coord"]," ",all[i,"Coord"]," ",all[i+1,"Coord"]," ",
                all[i+2,"Coord"]," ", mean(c(all[i,"Coord"],all[i+1,"Coord"]))," ",previous_position,sep=""))
    if((all[i,"P2_H1"] != all[i+1,"P2_H1"])&
       (as.numeric(all[i,]$pos) != 0)&
       (all[i,]$ind == previous_individual)){ # recombination occurred here in the father
      # populate a vector with the locations of maternal recombination events
      # print(paste(all[i,"P2_H1"]," ",all[i+1,"P2_H1"]," ",mean(c(all[i,"Coord"],all[i+1,"Coord"])),sep=""))
      # test whether this marker conflicts with the next marker and if so whether it it nearby
      if((all[i-1,"P2_H1"] == all[i,"P2_H1"]) & # these five conditions mean that at least two previous and two consecutive markers support a recombination event
         (all[i+1,"P2_H1"] == all[i+2,"P2_H1"]) &
         (all[i-1,"ind"] == all[i,"ind"]) & 
         (all[i+1,"ind"] == all[i+2,"ind"]) &
         (all[i,"ind"] == all[i+1,"ind"]) &
         (((as.numeric(all[i,"Coord"])) - as.numeric(previous_position)) > buffer)
         ) # this means that the previous recombination event was at least $buffer before this one
        { 
        #print(paste(all[i,"ind"]," ",all[i,"grp"]," ",all[i-1,"P2_H1"]," ",all[i,"P2_H1"]," ",all[i+1,"P2_H1"]," ",all[i+2,"P2_H1"]," ",
        #            all[i-1,"Coord"]," ",all[i,"Coord"]," ",all[i+1,"Coord"]," ",all[i+2,"Coord"]," ",previous_position,sep=""))
        all_recomb[(nrow(all_recomb) + 1),"Positions"] <- mean(c(all[i,"Coord"],all[i+1,"Coord"]))
        all_recomb[nrow(all_recomb),"Proportions"] <- mean(c(all[i,"Proportions"],all[i+1,"Proportions"]))
        all_recomb[nrow(all_recomb),"Proportions_relative_to_centromeres"] <- mean(c(all[i,"Proportions_relative_to_centromeres"],all[i+1,"Proportions_relative_to_centromeres"]))
        all_recomb[nrow(all_recomb),"Parent"] <- all[i,"matpat"] 
        all_recomb[nrow(all_recomb),"Chr"] <- all[i,"Chr"] 
        previous_chr <- all[i,]$grp
        previous_position <- mean(c(all[i,"Coord"],all[i+1,"Coord"]))
        previous_individual <- all[i,"ind"]
        pat_recomb_per_chr[all[i,]$Chr,all[i,"ind"]] <-  pat_recomb_per_chr[all[i,]$Chr,all[i,"ind"]]+1
      }  
    }
  } # end of check for paternal recomb
}


# This is the number of maternal recombination events per chr for each indiv:
View(mat_recomb_per_chr)
# this is a vector of the number of recombination events per individual across the whole genome
colSums(as.data.frame(mat_recomb_per_chr))
# This is the number of paternal recombination events per chr for each indiv:
View(pat_recomb_per_chr)
colSums(as.data.frame(pat_recomb_per_chr))

# number of recombinations per chr ----
library(reshape2)
mat_recomb_per_chr$Chrs <- rownames(mat_recomb_per_chr)
mat_recomb_per_chr_long <- gather(mat_recomb_per_chr, individual, count, unique(all$ind))

pat_recomb_per_chr$Chrs <- rownames(pat_recomb_per_chr)
pat_recomb_per_chr_long <- gather(pat_recomb_per_chr, individual, count, unique(all$ind))



dim(all_recomb)
head(all_recomb)

# Number of maternal recombination events across all individuals and all chrs
n_mat_recomb <- nrow(all_recomb[all_recomb$Parent == "mat",]);n_mat_recomb
# Number of paternal recombination events across all individuals and all chrs
n_pat_recomb <- nrow(all_recomb[all_recomb$Parent == "pat",]);n_pat_recomb


# number of maternal and paternal recombination event s ----
print(paste("Number mat recombination events across all individuals and all chrs ",n_mat_recomb,sep=" "))
print(paste("Number pat recombination events across all individuals and all chrs ",n_pat_recomb,sep=" "))

print(paste("Average number mat recombination events per individual across all chrs ",n_mat_recomb/length(unique(all$ind)),sep=" "))
print(paste("Average number pat recombination events per individual across all chrs ",n_pat_recomb/length(unique(all$ind)),sep=" "))


# Length of maternal LG
mat_bp_length <- 0
pat_bp_length <- 0
mat_cM_length <- 0
pat_cM_length <- 0
mat_bp_length_per_individual <- 0
pat_bp_length_per_individual <- 0
mat_cM_length_per_individual <- 0
pat_cM_length_per_individual <- 0
# make a df that has the max and min coordinates of each LG for each chromosome
# do this only for individual 1 because the others should be the same
LG_min_max <- data.frame(CHR = unique(all$Chr), 
                         min = rep(NA,length(unique(all$Chr))), 
                         max = rep(NA,length(unique(all$Chr))),
                         min_prop = rep(NA,length(unique(all$Chr))), 
                         max_prop = rep(NA,length(unique(all$Chr))),
                         chr_length = rep(NA,length(unique(all$Chr))))
for(j in unique(all$ind)) { # cycle through individuals
  print(j)
  for(i in unique(all$Chr)) { # cycle through chromosomes
    print(i)
    # add the max cM for each Chr for each individual; pos is in cM; Coord is bp positions on chr
    mat_cM_length <- mat_cM_length + max(all[(all$ind == j)&(all$Chr == i)&(all$matpat == "mat"),]$pos);mat_cM_length
    pat_cM_length <- pat_cM_length + max(all[(all$ind == j)&(all$Chr == i)&(all$matpat == "pat"),]$pos);pat_cM_length
    # now add the difference in min/max bp coordinates for each Chr for each individual
    mat_bp_length <- mat_bp_length + (max(all[(all$ind == j)&(all$Chr == i)&(all$matpat == "mat"),]$Coord)-min(all[(all$ind == j)&(all$Chr == i)&(all$matpat == "mat"),]$Coord))
    pat_bp_length <- pat_bp_length + (max(all[(all$ind == j)&(all$Chr == i)&(all$matpat == "pat"),]$Coord)-min(all[(all$ind == j)&(all$Chr == i)&(all$matpat == "pat"),]$Coord))
    if(j == unique(all$ind)[1]){ # get values for the first individual
      mat_bp_length_per_individual <- mat_bp_length_per_individual + (max(all[(all$ind == j)&(all$Chr == i)&(all$matpat == "mat"),]$Coord)-min(all[(all$ind == j)&(all$Chr == i)&(all$matpat == "mat"),]$Coord))
      pat_bp_length_per_individual <- pat_bp_length_per_individual + (max(all[(all$ind == j)&(all$Chr == i)&(all$matpat == "pat"),]$Coord)-min(all[(all$ind == j)&(all$Chr == i)&(all$matpat == "pat"),]$Coord))
      mat_cM_length_per_individual <- mat_cM_length_per_individual + max(all[(all$ind == j)&(all$Chr == i)&(all$matpat == "mat"),]$pos);mat_cM_length
      pat_cM_length_per_individual <- pat_cM_length_per_individual + max(all[(all$ind == j)&(all$Chr == i)&(all$matpat == "pat"),]$pos);pat_cM_length
      #print(paste(i,j,mat_cM_length_per_individual,sep=" "))
      # getting coordinates range for mat recombination here
      LG_min_max[(LG_min_max$CHR == i),"min"] <- min(all[(all$ind == j)&(all$Chr == i)&(all$matpat == "mat"),]$Coord)
      LG_min_max[(LG_min_max$CHR == i),"max"] <- max(all[(all$ind == j)&(all$Chr == i)&(all$matpat == "mat"),]$Coord)
    }
  }
}  

#print(paste("Inaccurate mat_cM_length summed over all individuals",mat_cM_length,sep=" "))
#print(paste("Inaccurate pat_cM_length summed over all individuals",pat_cM_length,sep=" "))
#print(paste("mat_bp_length summed over all individuals",mat_bp_length,sep=" "))
#print(paste("pat_bp_length summed over all individuals",pat_bp_length,sep=" "))

# the first two metrics are inaccurate because they include recombination events
# that were not well supported by the buffer and strings of consistent genotype criteria
#print(paste("Inaccurate mat_cM_length_for_one_individual in cM ",mat_cM_length_per_individual,sep=" "))
#print(paste("Inaccurate pat_cM_length_for_one_individual in cM ",pat_cM_length_per_individual,sep=" "))
print(paste("mat_bp_length_for_one_individual in bp ",mat_bp_length_per_individual/1000000000,' Gb',sep=" "))
print(paste("pat_bp_length_for_one_individual in bp ",pat_bp_length_per_individual/1000000000,' Gb',sep=" "))
# a centimorgand is length of chromosome in which an average of 0.01 crossover occurs per generation.

# according to wikipedia: https://en.wikipedia.org/wiki/Centimorgan
# the distance in cM is equal to 50*ln(1 / (1-2P) ); where p is the probability of recombination
# but this doesn't work for long distances where P>0.5 because the denominator is negative

# according to ChatGPT (sketchy!), cM = (# recombination events / total gametes) * 100

P_mat <- ((n_mat_recomb/length(unique(all$ind))))*100;P_mat
P_pat <- ((n_pat_recomb/length(unique(all$ind))))*100;P_pat

# It can be calculated as the length in bp/(# crossover events/100)
print(paste("Accurate mat_cM_length_for_one_individual in cM ",
            P_mat,
            sep=" "))
print(paste("Accurate pat_cM_length_for_one_individual in cM ",
            P_pat,
            sep=" "))
# human genome is ~3300 centimorgans, so these numbers for xennies are either way off or much bigger.
# maybe more recombination is needed to maintain disomy in a polyploid
# ******** check calculations *********
# for XL the "accurate" pat cM is lower than the "inaccurate" pat cM. It should be higher!
# ******** ******** ******** ******** ******** ******** ******** ******** ******** ********



# now calculate proportion boundaries of each LG for XL
for(i in LG_min_max$CHR) {
  if(i == "Chr1L"){
    LG_min_max[which(LG_min_max$CHR == i),]$chr_length <- 233740090
  } 
  else if(i == "Chr1S"){
    LG_min_max[which(LG_min_max$CHR == i),]$chr_length <- 202412970
  } 
  else if(i == "Chr2L"){
    LG_min_max[which(LG_min_max$CHR == i),]$chr_length <- 191000146
  } 
  else if(i == "Chr2S"){
    LG_min_max[which(LG_min_max$CHR == i),]$chr_length <- 169306100
  } 
  else if(i == "Chr3L"){
    LG_min_max[which(LG_min_max$CHR == i),]$chr_length <- 161426101
  } 
  else if(i == "Chr3S"){
    LG_min_max[which(LG_min_max$CHR == i),]$chr_length <- 131962816
  } 
  else if(i == "Chr4L"){
    LG_min_max[which(LG_min_max$CHR == i),]$chr_length <- 155250554
  } 
  else if(i == "Chr4S"){
    LG_min_max[which(LG_min_max$CHR == i),]$chr_length <- 132731174
  } 
  else if(i == "Chr5L"){
    LG_min_max[which(LG_min_max$CHR == i),]$chr_length <- 171415384
  } 
  else if(i == "Chr5S"){
    LG_min_max[which(LG_min_max$CHR == i),]$chr_length <- 143394103
  } 
  else if(i == "Chr6L"){
    LG_min_max[which(LG_min_max$CHR == i),]$chr_length <- 164223595
  } 
  else if(i == "Chr6S"){
    LG_min_max[which(LG_min_max$CHR == i),]$chr_length <- 137316286
  } 
  else if(i == "Chr7L"){
    LG_min_max[which(LG_min_max$CHR == i),]$chr_length <- 139837618
  } 
  else if(i == "Chr7S"){
    LG_min_max[which(LG_min_max$CHR == i),]$chr_length <- 113060389
  } 
  else if(i == "Chr8L"){
    LG_min_max[which(LG_min_max$CHR == i),]$chr_length <- 135449133
  } 
  else if(i == "Chr8S"){
    LG_min_max[which(LG_min_max$CHR == i),]$chr_length <- 103977862
  } 
  else if(i == "Chr9_10L"){
    LG_min_max[which(LG_min_max$CHR == i),]$chr_length <- 137811819
  } 
  else if(i == "Chr9_10S"){
    LG_min_max[which(LG_min_max$CHR == i),]$chr_length <- 117266291
  } 
}


# now calculate proportion boundaries of each LG for XT
for(i in LG_min_max$CHR) {
  if(i == "Chr1"){
    LG_min_max[which(LG_min_max$CHR == i),]$chr_length <- 217471165
  } 
  else if(i == "Chr2"){
    LG_min_max[which(LG_min_max$CHR == i),]$chr_length <- 181034960
  } 
  else if(i == "Chr3"){
    LG_min_max[which(LG_min_max$CHR == i),]$chr_length <- 153873356
  } 
  else if(i == "Chr4"){
    LG_min_max[which(LG_min_max$CHR == i),]$chr_length <- 153961318
  } 
  else if(i == "Chr5"){
    LG_min_max[which(LG_min_max$CHR == i),]$chr_length <- 164033574
  } 
  else if(i == "Chr6"){
    LG_min_max[which(LG_min_max$CHR == i),]$chr_length <- 154486311
  } 
  else if(i == "Chr7"){
    LG_min_max[which(LG_min_max$CHR == i),]$chr_length <- 133565929
  } 
  else if(i == "Chr8"){
    LG_min_max[which(LG_min_max$CHR == i),]$chr_length <- 147241509
  } 
  else if(i == "Chr9"){
    LG_min_max[which(LG_min_max$CHR == i),]$chr_length <- 91218943
  } 
  else if(i == "Chr10"){
    LG_min_max[which(LG_min_max$CHR == i),]$chr_length <- 52432565
  } 
}


# proportion of each chromosome covered by the matpat LG ----
LG_min_max$min_prop <- LG_min_max$min/LG_min_max$chr_length
LG_min_max$max_prop <- LG_min_max$max/LG_min_max$chr_length

LG_min_max$total_prop <- LG_min_max$max_prop - LG_min_max$min_prop;LG_min_max

LG_min_max$CHR <- factor(LG_min_max$CHR, 
                                       levels = c('Chr1L','Chr2L','Chr3L','Chr4L','Chr5L',
                                                  'Chr6L','Chr7L','Chr8L','Chr9_10L',
                                                  'Chr1S','Chr2S','Chr3S','Chr4S','Chr5S',
                                                  'Chr6S','Chr7S','Chr8S','Chr9_10S',
                                                  'Chr1','Chr2','Chr3','Chr4','Chr5',
                                                  'Chr6','Chr7','Chr8','Chr9','Chr10'),
                                       ordered = T)

# Plot the proportions of each chr that are covered by the LG
proportions <- ggplot(LG_min_max, aes(x=CHR))+
  geom_linerange(aes(ymin=min_prop,ymax=max_prop),linetype=2,color="blue")+
  geom_point(aes(y=min_prop),size=3,color="red")+
  geom_point(aes(y=max_prop),size=3,color="red")+
  xlab("Chromosome") + ylab("Proportions covered by largest linkage group") +
  theme_bw();proportions

ggsave(file="LG_proportion_of_CHRs.pdf", w=10, h=6, proportions)



# Plot proportion on chromosomes irrespective of centromere
densities <- ggplot(all_recomb, aes(x = Proportions)) +
  geom_density(aes(color = Parent))+
  #geom_vline(xintercept=49000000)+ 
  #facet_wrap( ~ Chr, ncol=2, scales = "free_x")+
  xlab("Relative chromosome position") + ylab("Density") +
  xlim(0,1) +
  theme_classic(); densities

ggsave(file="Recombination_density_relative_to_chr_tips.pdf", w=10, h=6, densities)



# Plot proportion on chromosomes relative of centromere (0 is a tip, 100 is a centromere)
ggplot(all_recomb, aes(x = Proportions_relative_to_centromeres)) +
  geom_density(aes(color = Parent), adjust=0.5)+
  #geom_vline(xintercept=49000000)+ 
  #facet_wrap( ~ Chr, ncol=2, scales = "free_x")+
  xlab("Proportion to Centromere (Tip = 0, Centromere = 1)") + ylab("Density") +
  xlim(0,1) +
  theme_classic()


# Plot proportion on chromosomes relative of centromere (0 is a tip, 100 is a centromere)
# stacked density
ggplot(all_recomb, aes(x = Proportions_relative_to_centromeres)) +
  geom_density(aes(color = Parent), position="stack")+
  xlab("Proportion to Centromere (Tip = 0, Centromere = 1)") + ylab("Density") +
  xlim(0,1) +
  theme_classic()

# all together
#ggplot(all_recomb, aes(x = Proportions_relative_to_centromeres)) +
#  geom_density()+
#  xlab("Proportion to Centromere (Tip = 0, Centromere = 1)") + ylab("Density") +
#  xlim(0,1) +
#  theme_classic()



#ggplot(all_recomb, aes(x=Proportions_relative_to_centromeres, group=Parent, fill=Parent)) +
#  geom_density(adjust=1.5, position="fill") +
#  xlab("Proportion to Centromere (Tip = 0, Centromere = 1)") + ylab("Density") +
#  theme_classic()


# Histogram with density plot
ggplot(all_recomb, aes(x=Proportions_relative_to_centromeres, color=Parent, fill=Parent)) + 
  geom_histogram(aes(y=..density..), alpha=0.5, 
                 position="identity")+
  geom_density(alpha=.2, adjust=0.5) +
  theme_classic()

# Histogram 
overlay_hist <- ggplot(all_recomb,aes(x=Proportions_relative_to_centromeres)) + 
  geom_histogram(data=subset(all_recomb,Parent == 'mat'),fill = "red", alpha = 0.2, binwidth = .025) +
  geom_histogram(data=subset(all_recomb,Parent == 'pat'),fill = "blue", alpha = 0.2, binwidth = .025) +
  xlab("Proportion to Centromere (Tip = 0, Centromere = 1)") + ylab("Count") +
  theme_classic();overlay_hist


ggsave(file="Recombination_events_relative_to_centromere.pdf", w=10, h=6, overlay_hist)

# Histogram 
overlay_hist <- ggplot(all_recomb,aes(x=Proportions)) + 
  geom_histogram(data=subset(all_recomb,Parent == 'mat'),fill = "red", alpha = 0.2, binwidth = .025) +
  geom_histogram(data=subset(all_recomb,Parent == 'pat'),fill = "blue", alpha = 0.2, binwidth = .025) +
  xlab("Proportions (Tip = 0, Other tip = 1)") + ylab("Count") +
  theme_classic();overlay_hist


ggsave(file="Recombination_events_on_scaled_chrs.pdf", w=10, h=6, overlay_hist)

# Now make a plot that shows how recombination events accumulate
# over the length of a chromosome
# first subset the sex chromosome
mysex_chr <- "Chr2L"
sex_chr_mat_recomb <- all_recomb[((all_recomb$Parent == "mat")&(all_recomb$Chr == mysex_chr)),];sex_chr_mat_recomb
sex_chr_pat_recomb <- all_recomb[((all_recomb$Parent == "pat")&(all_recomb$Chr == mysex_chr)),];sex_chr_pat_recomb
not_sex_chr_mat_recomb <- all_recomb[((all_recomb$Parent == "mat")&(all_recomb$Chr != mysex_chr)),];sex_chr_mat_recomb
not_sex_chr_pat_recomb <- all_recomb[((all_recomb$Parent == "pat")&(all_recomb$Chr != mysex_chr)),];sex_chr_pat_recomb

# now sort by chr and position
sex_chr_mat_recomb_sorted <- sex_chr_mat_recomb[order(sex_chr_mat_recomb$Chr, 
                                                      sex_chr_mat_recomb$Positions),];sex_chr_mat_recomb_sorted
sex_chr_pat_recomb_sorted <- sex_chr_pat_recomb[order(sex_chr_pat_recomb$Chr, 
                                                      sex_chr_pat_recomb$Positions),];sex_chr_pat_recomb_sorted
not_sex_chr_mat_recomb_sorted <- not_sex_chr_mat_recomb[order(not_sex_chr_mat_recomb$Positions),];not_sex_chr_mat_recomb_sorted
not_sex_chr_pat_recomb_sorted <- not_sex_chr_pat_recomb[order(not_sex_chr_pat_recomb$Positions),];not_sex_chr_pat_recomb_sorted

# add a column for the y-axis
sex_chr_mat_recomb_sorted$n_recombination <- 1:nrow(sex_chr_mat_recomb_sorted)
sex_chr_pat_recomb_sorted$n_recombination <- 1:nrow(sex_chr_pat_recomb_sorted)
not_sex_chr_mat_recomb_sorted$n_recombination <- 1:nrow(not_sex_chr_mat_recomb_sorted)
not_sex_chr_pat_recomb_sorted$n_recombination <- 1:nrow(not_sex_chr_pat_recomb_sorted)
# make column for color
sex_chr_mat_recomb_sorted$sexchr_or_not <- "red"
sex_chr_pat_recomb_sorted$sexchr_or_not <- "red"
not_sex_chr_mat_recomb_sorted$sexchr_or_not <- "black"
not_sex_chr_pat_recomb_sorted$sexchr_or_not <- "black"
# make column for matpat
sex_chr_mat_recomb_sorted$matpat <- "mat"
sex_chr_pat_recomb_sorted$matpat <- "pat"
not_sex_chr_mat_recomb_sorted$matpat <- "mat"
not_sex_chr_pat_recomb_sorted$matpat <- "pat"

# bind it all together
recomb_df <- rbind(sex_chr_mat_recomb_sorted,not_sex_chr_mat_recomb_sorted,
                   sex_chr_pat_recomb_sorted,not_sex_chr_pat_recomb_sorted)

# make the chr into a factor
recomb_df$Chr <- factor(recomb_df$Chr, 
                         levels = c('Chr1L','Chr2L','Chr3L','Chr4L','Chr5L',
                                    'Chr6L','Chr7L','Chr8L','Chr9_10L',
                                    'Chr1S','Chr2S','Chr3S','Chr4S','Chr5S',
                                    'Chr6S','Chr7S','Chr8S','Chr9_10S',
                                    'Chr1','Chr2','Chr3','Chr4','Chr5',
                                    'Chr6','Chr7','Chr8','Chr9','Chr10'),
                         ordered = T)

# plot
MatPat_recombination_accumulation_perchr<-ggplot(recomb_df, aes(x=Proportions, y=n_recombination, color=sexchr_or_not)) + 
  #geom_line() +
  geom_point(size=1, alpha = 0.7) +
  # get rid of gray background
  theme_bw() + theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
  labs(x = "Position", y= "Number of recombination events") +
  #scale_colour_hue(l=50) + # Use a slightly darker palette than normal
  theme(axis.text=element_text(size=5), axis.title=element_text(size=10)) +
  facet_grid(matpat ~ Chr) +
  # get rid of legend title
  # theme(legend.title = element_blank()) +
  # remove  the legend
  theme(legend.position="none") +
  # increase font size
  #theme(legend.text=element_text(size=14)) +
  # fix x-axis labels
  scale_x_continuous(breaks=seq(0,1,by=0.5)) +
  # color the stuff the way I want
  scale_color_manual(breaks = c("black", "red"),values=c("black", "red"),labels=c("Autosomes", "Sex Chr")) +
  theme(panel.background=element_rect(fill="transparent",colour="transparent"),
        plot.background=element_rect(fill="transparent",colour="transparent"),
        legend.key = element_rect(fill = "transparent", colour = "transparent"));MatPat_recombination_accumulation_perchr

ggsave(file="MatPat_recombination_accumulation_perchr.pdf", w=10, h=4, MatPat_recombination_accumulation_perchr)


# Histogram of recombination events per chr for mat ----
mat_recomb_per_chr_plot <- ggplot(mat_recomb_per_chr_long) +
    geom_histogram(aes(x=count), position="dodge", binwidth = 1) +
    scale_x_continuous(breaks=seq(0, 20, 2)) +
    scale_y_continuous(breaks=seq(0, 20, 5)) +
    facet_wrap(~Chrs, nrow=9,ncol=2) +
    labs(x = "Recombinations per chromosome") +
    theme_classic(); mat_recomb_per_chr_plot  
  
ggsave(file="Mat_recombination_perchr.pdf", w=4, h=10, mat_recomb_per_chr_plot)

# Histogram of recombination events per chr for pat ----
pat_recomb_per_chr_plot <- ggplot(pat_recomb_per_chr_long) +
  geom_histogram(aes(x=count), position="dodge", binwidth = 1) +
  scale_x_continuous(breaks=seq(0, 20, 2)) +
  scale_y_continuous(breaks=seq(0, 20, 5)) +
  facet_wrap(~Chrs, nrow=9,ncol=2) +
  labs(x = "Recombinations per chromosome") +
  theme_classic(); pat_recomb_per_chr_plot

ggsave(file="Pat_recombination_perchr.pdf", w=4, h=6, pat_recomb_per_chr_plot)


```


