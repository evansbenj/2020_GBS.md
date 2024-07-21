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
realSFS fst index fem_bamz_Chr3L_.saf.idx mal_bamz_Chr3L_.saf.idx -sfs fem_mal_Chr3L.ml -fstout fem_mal_Chr3L__
```
now step4:
```
realSFS fst stats2 fem_mal_Chr3L__.fst.idx -win 500000 -step 500000 > fisch_fem_mal_fst_500000_windows.txt
```

# Plotting
```R
setwd("/Users/Shared/Previously Relocated Items/Security/projects/submitted/2022_GBS_lotsof_Xennies/2024_fst_windows")
library(ggplot2)
library(readr)
library(dplyr)
library(ggh4x)
dir <- "/Users/Shared/Previously Relocated Items/Security/projects/submitted/2022_GBS_lotsof_Xennies/2024_fst_windows"
list.files(dir)

# load the data 

allo_nonkin <- read.table("allo_nonkin_fem_mal_fst_5000000_windows.txt", header = F, skip = 1)
colnames(allo_nonkin) <- c("region","chr","midPos","NSites","Fst")
allo_nonkin$kin_nonkin <- "nonkin"
allo_nonkin$species <- "allo"

allo_kin <- read.table("allo_kin_fem_mal_fst_5000000_windows.txt", header = F, skip = 1)
colnames(allo_kin) <- c("region","chr","midPos","NSites","Fst")
allo_kin$kin_nonkin <- "kin"
allo_kin$species <- "allo"

bor_wildwest <- read.table("bor_wildwest_fem_mal_fst_5000000_windows.txt", header = F, skip = 1)
colnames(bor_wildwest) <- c("region","chr","midPos","NSites","Fst")
bor_wildwest$kin_nonkin <- "nonkin"
bor_wildwest$species <- "bore"

bor_wildeast <- read.table("bor_wildeast_fem_mal_fst_5000000_windows.txt", header = F, skip = 1)
colnames(bor_wildeast) <- c("region","chr","midPos","NSites","Fst")
bor_wildeast$kin_nonkin <- "nonkin"
bor_wildeast$species <- "bore"

bor_kin <- read.table("borfam_fem_mal_fst_5000000_windows.txt", header = F, skip = 1)
colnames(bor_kin) <- c("region","chr","midPos","NSites","Fst")
bor_kin$kin_nonkin <- "kin"
bor_kin$species <- "bore"

boum_nonkin <- read.table("boum_nonkin_fem_mal_fst_5000000_windows.txt", header = F, skip = 1)
colnames(boum_nonkin) <- c("region","chr","midPos","NSites","Fst")
boum_nonkin$kin_nonkin <- "nonkin"
boum_nonkin$species <- "boum"

boum_kin <- read.table("boum_kin_fem_mal_fst_5000000_windows.txt", header = F, skip = 1)
colnames(boum_kin) <- c("region","chr","midPos","NSites","Fst")
boum_kin$kin_nonkin <- "kin"
boum_kin$species <- "boum"

muel_nonkin <- read.table("muel_fem_mal_nonkinonly_fst_5000000_windows.txt", header = F, skip = 1)
colnames(muel_nonkin) <- c("region","chr","midPos","NSites","Fst")
muel_nonkin$kin_nonkin <- "nonkin"
muel_nonkin$species <- "muel"

muel_kin <- read.table("muel_kin_fem_mal_fst_5000000_windows.txt", header = F, skip = 1)
colnames(muel_kin) <- c("region","chr","midPos","NSites","Fst")
muel_kin$kin_nonkin <- "kin"
muel_kin$species <- "muel"

fisc_nonkin <- read.table("fisch_fem_mal_onlynonkin_fst_5000000_windows.txt", header = F, skip = 1)
colnames(fisc_nonkin) <- c("region","chr","midPos","NSites","Fst")
fisc_nonkin$kin_nonkin <- "nonkin"
fisc_nonkin$species <- "fisc"

fisc_kin <- read.table("fisch_kin_fem_mal_fst_5000000_windows.txt", header = F, skip = 1)
colnames(fisc_kin) <- c("region","chr","midPos","NSites","Fst")
fisc_kin$kin_nonkin <- "kin"
fisc_kin$species <- "fisc"

pygm_nonkin <- read.table("pygm_fem_mal_fst_5000000_windows_nonkin.txt", header = F, skip = 1)
colnames(pygm_nonkin) <- c("region","chr","midPos","NSites","Fst")
pygm_nonkin$kin_nonkin <- "nonkin"
pygm_nonkin$species <- "pygm"

pygm_kin <- read.table("pygm_kin_fem_mal_fst_5000000_windows.txt", header = F, skip = 1)
colnames(pygm_kin) <- c("region","chr","midPos","NSites","Fst")
pygm_kin$kin_nonkin <- "kin"
pygm_kin$species <- "pygm"

laev_nonkin <- read.table("laev_wild_fem_mal_fst_5000000_windows.txt", header = F, skip = 1)
colnames(laev_nonkin) <- c("region","chr","midPos","NSites","Fst")
laev_nonkin$kin_nonkin <- "nonkin"
laev_nonkin$species <- "laev"

laev_kin <- read.table("laev_lab_fem_mal_fst_5000000_windows.txt", header = F, skip = 1)
colnames(laev_kin) <- c("region","chr","midPos","NSites","Fst")
laev_kin$kin_nonkin <- "kin"
laev_kin$species <- "laev"

lend_nonkin <- read.table("lend_nonkin_fem_mal_fst_5000000_windows.txt", header = F, skip = 1)
colnames(lend_nonkin) <- c("region","chr","midPos","NSites","Fst")
lend_nonkin$kin_nonkin <- "nonkin"
lend_nonkin$species <- "lend"


my_df <- rbind(allo_kin, allo_nonkin,boum_kin, boum_nonkin,
               bor_kin,bor_wildeast,fisc_kin,fisc_nonkin,laev_kin,laev_nonkin,
               lend_nonkin,
               muel_kin,muel_nonkin,pygm_kin,pygm_nonkin)

head(my_df)
my_df$midPos <- as.numeric(my_df$midPos)
my_df$Fst <- as.numeric(my_df$Fst)

my_df$color <- "black"
my_df$color[my_df$Fst > 0.25 ] <- "orange"
my_df$color[my_df$Fst > 0.5 ] <- "red"

# make rectangles to highlight SL regions for each species
rect<-data.frame(xmin = c(8,8,117,117,1,1,41,41,177,177,0,15,111,111,117,117), 
                 xmax = c(17,17,135,135,54.1, 54.1,105,105,190,190,0,16.2,147,147,135,135), 
                 ymin = c(rep(-0.1,16)), ymax = c(rep(1,16)), 
                 alpha = c(rep(0.1,16)),
                 fill = c(rep("blue",16)))

pdf("./Xenopus_Fst.pdf",w=6, h=8.0, version="1.4", bg="transparent")
  p<-ggplot(my_df, aes(x=midPos/1000000, y=Fst, color = color)) + 
    # add points
    geom_point(size=2, alpha = 0.7) +
    # color the stuff the way I want
    scale_color_manual(name=expression(paste(italic(F[ST]))),
                       values = c("black" = "lightgray", "orange" = "orange","red" = "red"),
                       labels = c("<0.5", "0.5-0.8", ">0.8"))+
    #ylim(-0.1,1) +
    xlab("Position (Mb)") + ylab(expression(paste(italic(F[ST])))) +
    # highlight the region of recombination suppression
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              alpha = 0.1, fill = "blue",
              data = transform(rect, species = c("allo","allo","boum","boum",
                                                 "bore","bore","fisc","fisc",
                                                 "laev","laev","lend","lend",
                                                 "muel","muel",
                                                 "pygm","pygm"),
                               kin_nonkin = c("kin","nonkin","kin","nonkin",
                                              "kin","nonkin","kin","nonkin",
                                              "kin","nonkin","kin","nonkin",
                                              "kin","nonkin",
                                              "kin","nonkin")),
              inherit.aes = FALSE) +
    #annotate("rect", xmin = 0, xmax = 54.1, ymin = -0.1, ymax = 1,
    #         alpha = .1,fill = "blue") +
    ggh4x::facet_grid2(species ~ kin_nonkin, scales = "free_y", independent = "y") +
    #facet_wrap(kin_nonkin ~ species, scales = "free_y") +
    # get rid of gray background
    theme_bw() +
    theme(legend.position="none") +
    theme(strip.background = element_blank(), strip.placement = "outside")
    p
dev.off()

# subgenus Silurana
trop_nonkin <- read.table("trop_nonkin_fem_mal_fst_5000000_windows.txt", header = F, skip = 1)
colnames(trop_nonkin) <- c("region","chr","midPos","NSites","Fst")
trop_nonkin$kin_nonkin <- "nonkin"
trop_nonkin$species <- "trop"

trop_kin <- read.table("trop_kin_fem_mal_fst_5000000_windows.txt", header = F, skip = 1)
colnames(trop_kin) <- c("region","chr","midPos","NSites","Fst")
trop_kin$kin_nonkin <- "kin"
trop_kin$species <- "trop"

mell_kin <- read.table("mell_kin_fem_mal_fst_5000000_windows.txt", header = F, skip = 1)
colnames(mell_kin) <- c("region","chr","midPos","NSites","Fst")
mell_kin$kin_nonkin <- "kin"
mell_kin$species <- "mell"


my_df <- rbind(trop_kin, trop_nonkin,mell_kin)

head(my_df)
my_df$midPos <- as.numeric(my_df$midPos)
my_df$Fst <- as.numeric(my_df$Fst)

my_df$color <- "black"
my_df$color[my_df$Fst > 0.25 ] <- "orange"
my_df$color[my_df$Fst > 0.5 ] <- "red"

# make rectangles to highlight SL regions for each species
rect<-data.frame(xmin = c(0,0,0,0), 
                 xmax = c(11,0,11,11), 
                 ymin = c(rep(-0.1,4)), ymax = c(rep(1,4)), 
                 alpha = c(rep(0.1,4)),
                 fill = c(rep("blue",4)))

pdf("./Silurana_Fst.pdf",w=6, h=2.5, version="1.4", bg="transparent")
p<-ggplot(my_df, aes(x=midPos/1000000, y=Fst, color = color)) + 
  # add points
  geom_point(size=2, alpha = 0.7) +
  # color the stuff the way I want
  scale_color_manual(name=expression(paste(italic(F[ST]))),
                     values = c("black" = "lightgray", "orange" = "orange","red" = "red"),
                     labels = c("<0.5", "0.5-0.8", ">0.8"))+
  #ylim(-0.1,1) +
  xlab("Position (Mb)") + ylab(expression(paste(italic(F[ST])))) +
  # highlight the region of recombination suppression
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            alpha = 0.1, fill = "blue",
            data = transform(rect, species = c("mell","mell",
                                               "trop","trop"),
                             kin_nonkin = c("kin","nonkin","kin","nonkin")),
            inherit.aes = FALSE) +
  #annotate("rect", xmin = 0, xmax = 54.1, ymin = -0.1, ymax = 1,
  #         alpha = .1,fill = "blue") +
  ggh4x::facet_grid2(species ~ kin_nonkin, scales = "free_y", independent = "y") +
  #facet_wrap(kin_nonkin ~ species, scales = "free_y") +
  # get rid of gray background
  theme_bw() +
  theme(legend.position="none") +
  theme(strip.background = element_blank(), strip.placement = "outside")
p
dev.off()
```
