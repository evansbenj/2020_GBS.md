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

# Rscript "2023_OneMap.R" --args inputfile_C="allo_family_one_chr7L_filtered.vcf" mom_C="allo_Cam_female_4_F_AGGAT_TAATA_cuttrim_sorted.bam" dad_C="allo_Cam_male_1_M_GGTGT_GTCAA_cuttrim_sorted.bam" prefix_C="allo_1_chr7L"

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

 my_dat <- onemap_read_vcfR(vcf = "allo_family_one_chr7L_filtered.vcf",
                           cross = c("outcross"),
                           parent1 = c("allo_Cam_female_4_F_AGGAT_TAATA_cuttrim_sorted.bam"), 
                           parent2 = c("allo_Cam_male_1_M_GGTGT_GTCAA_cuttrim_sorted.bam"), 
                           only_biallelic = TRUE,
                           verbose = TRUE); my_dat

 my_dat <- onemap_read_vcfR(vcf = "DB__Chr5S_out.vcf_filtered.vcf.gz_selected.vcf",
                           cross = c("outcross"),
                           parent1 = c("./3897mom_trim_sorted"), 
                           parent2 = c("./3896dad_trim_sorted"), 
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

# get indexes of mat and pat sites
matpat_SNPs <- as.integer(p[(p$Type == 'B3.7')|(p$Type == 'D1.10')|(p$Type == 'D2.15'), ]$Marker)

# Marker types are explained here:
# https://statgen-esalq.github.io/tutorials/onemap/Outcrossing_Populations.html

# 1 B3.7    321  ab×ab (double het)
# 2 D1.10   483  ab×aa (first parent het; maternal)
# 3 D2.15   384  aa×ab (second parent het; paternal)

# make an object with only matpat SNPs
matpat_SNPs_no_dist <- make_seq(twopts, matpat_SNPs)

# Use the function suggest_lod to calculate a suggested LOD score considering that multiple tests are being performed.
LOD_sug_matpat <- suggest_lod(matpat_SNPs_no_dist)

# And apply this suggested value to the two-point tests and set the maximum recombination fraction to 0.40:
matpat_LGs <- group(matpat_SNPs_no_dist, LOD=LOD_sug_matpat, max.rf = 0.4)

# figure out which linkage group has the most markers
matpat_df <- as.data.frame(table(matpat_LGs$groups))

# get rid of first row, which is the unlinked markers
matpat_df <- matpat_df[-1, ]

matpat_biggest_LG_group <- as.vector(matpat_df$Var1[matpat_df$Freq == max(matpat_df$Freq)]);matpat_biggest_LG_group

matpat_biggest_LG <- make_seq(matpat_LGs,as.integer(matpat_biggest_LG_group))

matpat_LG_1_graph <- rf_graph_table(matpat_biggest_LG)
ggsave(file=paste(prefix_C,"_matpat_graph_LG.pdf",sep="_"), matpat_LG_1_graph, width=8, height=8)

# make a combined map
matpat_map <- onemap::map(matpat_biggest_LG)
marker_type(matpat_map)

q <- marker_type(matpat_map)
q %>% group_by(Type) %>% count()
  
print(matpat_map, detailed = T)
  
matpat_biggest_LG_matonly <- marker_type(matpat_biggest_LG)[(marker_type(matpat_biggest_LG)$Type == 'D1.10'), ]$Marker
matpat_biggest_LG_patonly <- marker_type(matpat_biggest_LG)[(marker_type(matpat_biggest_LG)$Type == 'D2.15'), ]$Marker


matpat_biggest_LG_matonly_no_dist <- make_seq(twopts, c(matpat_biggest_LG_matonly))
matpat_biggest_LG_patonly_no_dist <- make_seq(twopts, c(matpat_biggest_LG_patonly))

# look at markers
marker_type(matpat_biggest_LG_matonly_no_dist)
s <- marker_type(matpat_biggest_LG_matonly_no_dist)
s %>% group_by(Type) %>% count()

# look at markers
marker_type(matpat_biggest_LG_patonly_no_dist)
r <- marker_type(matpat_biggest_LG_patonly_no_dist)
r %>% group_by(Type) %>% count()


# Use the function suggest_lod to calculate a suggested LOD score considering that multiple tests are being performed.
LOD_sug_mat <- suggest_lod(matpat_biggest_LG_matonly_no_dist)
LOD_sug_pat <- suggest_lod(matpat_biggest_LG_patonly_no_dist)

# And apply this suggested value to the two-point tests and set the maximum recombination fraction to 0.40:
# mat_LGs_ord <- group(matpat_biggest_LG_matonly_no_dist, LOD=LOD_sug_mat, max.rf = 0.4)
# pat_LGs_ord <- group(matpat_biggest_LG_patonly_no_dist, LOD=LOD_sug_pat, max.rf = 0.4)

mat_biggest_LGs_ordd <- order_seq(matpat_biggest_LG_matonly_no_dist, n.init = 7, THRES = LOD_sug_mat)
pat_biggest_LGs_ordd <- order_seq(matpat_biggest_LG_patonly_no_dist, n.init = 7, THRES = LOD_sug_pat)

mat_LGs_ord_force <- make_seq(mat_biggest_LGs_ordd, "force")
pat_LGs_ord_force <- make_seq(pat_biggest_LGs_ordd, "force")

mat_LG_1_graph <- rf_graph_table(mat_LGs_ord_force)
mat_LG_1_graph_LOD <- rf_graph_table(mat_LGs_ord_force, graph.LOD = TRUE)
ggsave(file=paste(prefix_C,"_mat_graph_LG_longest_force.pdf",sep="_"), mat_LG_1_graph, width=8, height=8)
ggsave(file=paste(prefix_C,"_mat_graph_LG_longest_force_LOD.pdf",sep="_"), mat_LG_1_graph_LOD, width=8, height=8)

pat_LG_1_graph <- rf_graph_table(pat_LGs_ord_force)
pat_LG_1_graph_LOD <- rf_graph_table(pat_LGs_ord_force, graph.LOD = TRUE)
ggsave(file=paste(prefix_C,"_pat_graph_LG_longest_force.pdf",sep="_"), pat_LG_1_graph, width=8, height=8)
ggsave(file=paste(prefix_C,"_pat_graph_LG_longest_force_LOD.pdf",sep="_"), pat_LG_1_graph_LOD, width=8, height=8)


maternal_map <- onemap::map(mat_LGs_ord_force)
paternal_map <- onemap::map(pat_LGs_ord_force)






# get indexes of double hets or maternal specific hets:
#maternal_SNPs <- as.integer(p[(p$Type == 'D1.10'), ]$Marker)
# get indexes of double hets or paternal specific hets:
#paternal_SNPs <- as.integer(p[(p$Type == 'D2.15'), ]$Marker)


# make an object with only maternal SNPs
# maternal_SNPs_no_dist <- make_seq(twopts, maternal_SNPs)
# make an object with only paternal SNPs
# paternal_SNPs_no_dist <- make_seq(twopts, paternal_SNPs)


# Use the function suggest_lod to calculate a suggested LOD score considering that multiple tests are being performed.
# LOD_sug_mat <- suggest_lod(maternal_SNPs_no_dist)
# LOD_sug_pat <- suggest_lod(paternal_SNPs_no_dist)

# And apply this suggested value to the two-point tests and set the maximum recombination fraction to 0.40:
# mat_LGs <- group(maternal_SNPs_no_dist, LOD=LOD_sug_mat, max.rf = 0.4)
# pat_LGs <- group(paternal_SNPs_no_dist, LOD=LOD_sug_pat, max.rf = 0.4)

# figure out which linkage group has the most markers
# mat_df <- as.data.frame(table(mat_LGs$groups))
# pat_df <- as.data.frame(table(pat_LGs$groups))

# get rid of first row, which is the unlinked markers
# mat_df <- mat_df[-1, ]
# pat_df <- pat_df[-1, ]

# mat_biggest_LG_group <- as.vector(mat_df$Var1[mat_df$Freq == max(mat_df$Freq)]);mat_biggest_LG_group
# pat_biggest_LG_group <- as.vector(pat_df$Var1[pat_df$Freq == max(pat_df$Freq)]);pat_biggest_LG_group

# mat_biggest_LG <- make_seq(mat_LGs,as.integer(mat_biggest_LG_group))
# pat_biggest_LG <- make_seq(pat_LGs,as.integer(pat_biggest_LG_group))

# mat_LG_1_graph <- rf_graph_table(mat_biggest_LG)
# pat_LG_1_graph <- rf_graph_table(pat_biggest_LG)

# maternal_map <- onemap::map(mat_biggest_LG)
# paternal_map <- onemap::map(pat_biggest_LG)


maternal_parents_haplot <- parents_haplotypes(maternal_map)
paternal_parents_haplot <- parents_haplotypes(paternal_map)
# matpat_haplot <- parents_haplotypes(matpat_map)

# Add column for matpat
maternal_parents_haplot$matpat <- "mat"
paternal_parents_haplot$matpat <- "pat"

#prefix_C<-"temp"
write.table(maternal_parents_haplot, paste(prefix_C,"mat_parents_haplot_LG.txt",sep="_"), row.names=F)
write.table(paternal_parents_haplot, paste(prefix_C,"pat_parents_haplot_LG.txt",sep="_"), row.names=F)
# write.table(maternal_haplot, "mat_haplot.txt")
# write.table(paternal_haplot, "pat_haplot.txt")
# write.table(matpat_haplot, "matpat_haplot.txt")

# draw_map(maternal_map)

# Export haplotypes
mat_progeny_haplot <- progeny_haplotypes(maternal_map, 
                                     most_likely = TRUE, 
                                     ind = c(1:data_filtered$n.ind), 
                                     group_names = prefix_C)

pat_progeny_haplot <- progeny_haplotypes(paternal_map, 
                                         most_likely = TRUE, 
                                         ind = c(1:data_filtered$n.ind), 
                                         group_names = prefix_C)

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

write.table(mat_progeny_haplot_wide, paste(prefix_C,"mat_progeny_haplot_wide_LG.txt",sep="_"),row.names = F)
write.table(pat_progeny_haplot_wide, paste(prefix_C,"pat_progeny_haplot_wide_LG.txt",sep="_"),row.names = F)

# ggsave(file=paste(prefix_C,"_mat_graph_LG.pdf",sep="_"), mat_LG_1_graph, width=8, height=8)
# ggsave(file=paste(prefix_C,"_pat_graph_LG.pdf",sep="_"), pat_LG_1_graph, width=8, height=8)

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
# Concatenate all files but save the header
```sh
head -1 allo_1_chr1S_mat_progeny_haplot_wide.txt > all_mat.txt; awk 'FNR>1{print}' *mat*wide.txt >> all_mat.txt
```

# Processing script for OneMap output
```R
# concatenate "wide" haplotypes from all individuals for maternal and paternal recombination events like this:
#  head -1 one_of_the_files_haplot_wide.txt > all_mat.txt; awk 'FNR>1{print}' *mat*wide.txt >> all_mat.txt
#  head -1 one_of_the_files_pat_progeny_haplot_wide.txt > all_pat.txt; awk 'FNR>1{print}' *pat*wide.txt >> all_pat.txt
# if there are rownames, you can add a bogus header for them to prevent a redundant rowname error like this:
# perl -pi -e 's/^/"bogus" / if $.==1' all_mat.txt

setwd('/Users/Shared/Previously Relocated Items/Security/projects/2022_GBS_lotsof_Xennies/2023_OneMap/2023_muelleri')
library(tidyverse)
# read in the OneMap data 
all_mat <- read.table("./all_mat.txt", header = T)
all_pat <- read.table("./all_pat.txt", header = T)

# ad matpat column
all_mat$matpat <- "mat"
all_pat$matpat <- "pat"

all<- rbind(all_mat,all_pat)

# make a column with chrs and one with coordinates
# this probably needs to be different for trop
all[c('Chr', 'Coord')] <- str_split_fixed(all$marker, '_', 2)
all[all$Chr == 'Chr9',]$Coord <- as.numeric(str_split_i(all[all$Chr == 'Chr9',]$marker,"_",3)) # need to use a 3 because Chr9_10L has an underscore
all[all$Chr == 'Chr9',]$Chr <- paste(str_split_i(all[all$Chr == 'Chr9',]$marker, "_",1),str_split_i(all[all$Chr == 'Chr9',]$marker, "_",2), sep = "_") # need to use a 3 because Chr9_10L has an underscore
unique(all$Chr)

all$Coord <- as.numeric(all$Coord)

# now make a column with chromosome proportions
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
max(all$Proportions) # should < 1
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


    
all_recomb <- data.frame(Chr=character(),
                         Positions=integer(),
                         Proportions=integer(),
                         Parent=character(),
                         stringsAsFactors=F)
buffer <- 5000000
# Figure out where recombination occurred during oogenesis
for(i in 1:(nrow(all)-2)) {       # for-loop over rows; need to only go up to the second to last to allow
                                  # a check for at least two sites supporting recombination
  if(as.numeric(all[i,]$pos) == 0){
    switch=0
  }
  if(all[i,"P1_H1"] != all[i+1,"P1_H1"]){ # recombination occurred here
    # populate a vector with the locations of maternal recombination events
    print(paste(all[i,"P1_H1"]," ",all[i+1,"P1_H1"]," ",mean(all[i,"Coord"],all[i+1,"Coord"]),sep=""))
    # test whether this marker conflicts with the next marker and if so whether it it nearby
    if(switch == 0){
      all_recomb[(nrow(all_recomb) + 1),"Positions"] <- mean(all[i,"Coord"],all[i+1,"Coord"])
      all_recomb[nrow(all_recomb),"Proportions"] <- mean(all[i,"Proportions"],all[i+1,"Proportions"])
      all_recomb[nrow(all_recomb),"Parent"] <- all[i,"matpat"] 
      all_recomb[nrow(all_recomb),"Chr"] <- all[i,"Chr"] 
      switch=1
    }
    if((switch != 0)
       &&
       ((all[i+1,"P1_H1"] == all[i+2,"P1_H1"]) # this means that at least two consecutive markers support a recombination event
       &&
       (as.numeric(all[i,"Coord"]) - as.numeric(all_recomb[nrow(all_recomb),"Positions"]) > buffer)
       )) # this means that the previous recombination event was at least $buffer before this one
      { 
      all_recomb[(nrow(all_recomb) + 1),"Positions"] <- mean(all[i,"Coord"],all[i+1,"Coord"])
      all_recomb[nrow(all_recomb),"Proportions"] <- mean(all[i,"Proportions"],all[i+1,"Proportions"])
      all_recomb[nrow(all_recomb),"Parent"] <- all[i,"matpat"] 
      all_recomb[nrow(all_recomb),"Chr"] <- all[i,"Chr"] 
    }  
  }
}
dim(all_recomb)
head(all_recomb)

# Number of maternal recombination events
n_mat_recomb <- nrow(all_recomb[all_recomb$Parent == "mat",]);n_mat_recomb
# Number of maternal recombination events
n_pat_recomb <- nrow(all_recomb[all_recomb$Parent == "pat",]);n_pat_recomb

# Length of maternal LG
mat_bp_length <- 0
pat_bp_length <- 0
mat_cM_length <- 0
pat_cM_length <- 0
for(i in unique(all$Chr)) {
  print(i)
  # add the max cM for each Chr
  mat_cM_length <- mat_cM_length + max(all[(all$Chr == i)&(all$matpat == "mat"),]$pos);mat_cM_length
  pat_cM_length <- pat_cM_length + max(all[(all$Chr == i)&(all$matpat == "pat"),]$pos);pat_cM_length
  # now add the difference in min/max bp coordinates for each char
  mat_bp_length <- mat_bp_length + (max(all[(all$Chr == i)&(all$matpat == "mat"),]$Coord)-min(all[(all$Chr == i)&(all$matpat == "mat"),]$Coord))
  pat_bp_length <- pat_bp_length + (max(all[(all$Chr == i)&(all$matpat == "pat"),]$Coord)-min(all[(all$Chr == i)&(all$matpat == "pat"),]$Coord))
}

print(paste("mat_cM_length",mat_cM_length,sep=" "))
print(paste("pat_cM_length",pat_cM_length,sep=" "))
print(paste("mat_bp_length",mat_bp_length,sep=" "))
print(paste("pat_bp_length",pat_bp_length,sep=" "))



ggplot(all_recomb, aes(x = Proportions)) +
  geom_density(aes(color = Parent))+
  #geom_vline(xintercept=49000000)+ 
  #facet_wrap( ~ Chr, ncol=2, scales = "free_x")+
  xlab("Chromosome proportion") + ylab("Density") +
  xlim(0,1) +
  theme_bw()



```


