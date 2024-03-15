# Gene density

First harvest the start position of all annotated "genes" like this:
```
zcat < XENTR_10.0_Xenbase_longest.gff3.gz | grep 'gene        ' | cut -f1,4,9 > XENTR_10.0_Xenbase_longest.gff3_gene_starts.txt zcat < XENLA_10.1_Xenbase_longest.gff3.gz | grep 'gene        ' | cut -f1,4,9 > XENLA_10.1_Xenbase_longest.gff3_gene_starts.txt 
```
This resulted in 28864 for XENTR_10.0 and 44457 for XENLA_10.1. Perfect!

An interesting question is whether sex-linked regions are found in areas that have high gene density. To explore this I plotted the density of start sites of annotated genes in XL and the locations of sex-linked regions. 

Here is the R code for the XL plot:
```
setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/2023_angsd/mapped_to_XL_v10_only")
library(ggplot2)
library(readr)
library(dplyr)
dir <- "/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/2023_angsd/mapped_to_XL_v10_only"
list.files(dir)

# load the data 
my_gene_startsites <- read.table("chr_and_gene_start_longest.txt", header = F)
colnames(my_gene_startsites) <- c("CHR","POS")
my_gene_startsites_chronly$MB <- as.numeric(my_gene_startsites_chronly$POS)/1000000
# subset to only include large chrs
my_gene_startsites_chronly <- my_gene_startsites[(my_gene_startsites$CHR == "Chr1L")|
                                                   (my_gene_startsites$CHR == "Chr2L")|
                                                   (my_gene_startsites$CHR == "Chr3L")|
                                                   (my_gene_startsites$CHR == "Chr4L")|
                                                   (my_gene_startsites$CHR == "Chr5L")|
                                                   (my_gene_startsites$CHR == "Chr6L")|
                                                   (my_gene_startsites$CHR == "Chr7L")|
                                                   (my_gene_startsites$CHR == "Chr8L")|
                                                   (my_gene_startsites$CHR == "Chr9_10L")|
                                                   (my_gene_startsites$CHR == "Chr1S")|
                                                   (my_gene_startsites$CHR == "Chr2S")|
                                                   (my_gene_startsites$CHR == "Chr3S")|
                                                   (my_gene_startsites$CHR == "Chr4S")|
                                                   (my_gene_startsites$CHR == "Chr5S")|
                                                   (my_gene_startsites$CHR == "Chr6S")|
                                                   (my_gene_startsites$CHR == "Chr7S")|
                                                   (my_gene_startsites$CHR == "Chr8S")|
                                                   (my_gene_startsites$CHR == "Chr9_10S"),]
  
  factor(my_gene_startsites_chronly$CHR,
                    levels = c("Chr1L","Chr1S","Chr2L","Chr2S",
                               "Chr3L","Chr3S","Chr4L","Chr4S",
                               "Chr5L","Chr5S","Chr6L","Chr6S",
                               "Chr7L","Chr7S","Chr8L","Chr8S",
                               "Chr9_10L","Chr9_10S"), ordered = T)

library(ggplot2)
library(scales)
library(tidyverse)
  
# A data frame with labels for each facet
  
# allo
allo_labels <- data.frame(CHR = c("Chr1L","Chr1S","Chr2L","Chr2S",
                                    "Chr3L","Chr3S","Chr4L","Chr4S",
                                    "Chr5L","Chr5S","Chr6L","Chr6S",
                                    "Chr7L","Chr7S","Chr8L","Chr8S",
                                    "Chr9_10L","Chr9_10S"), 
                            label = c("","","","",
                                      "","","","",
                                      "","","","",
                                      "X. allofraseri","","","",
                                      "",""))   
# bore
bore_labels <- data.frame(CHR = c("Chr1L","Chr1S","Chr2L","Chr2S",
                                    "Chr3L","Chr3S","Chr4L","Chr4S",
                                    "Chr5L","Chr5S","Chr6L","Chr6S",
                                    "Chr7L","Chr7S","Chr8L","Chr8S",
                                    "Chr9_10L","Chr9_10S"), 
                            label = c("","","","",
                                      "","","","",
                                      "","","","",
                                      "","","X. borealis","",
                                      "",""))    
# fisc
fisc_labels <- data.frame(CHR = c("Chr1L","Chr1S","Chr2L","Chr2S",
                                   "Chr3L","Chr3S","Chr4L","Chr4S",
                                   "Chr5L","Chr5S","Chr6L","Chr6S",
                                   "Chr7L","Chr7S","Chr8L","Chr8S",
                                   "Chr9_10L","Chr9_10S"), 
                           label = c("","","","",
                                     "X. fischbergi","","","",
                                     "","","","",
                                     "","","","",
                                     "",""))

# laev
laev_labels <- data.frame(CHR = c("Chr1L","Chr1S","Chr2L","Chr2S",
                                  "Chr3L","Chr3S","Chr4L","Chr4S",
                                  "Chr5L","Chr5S","Chr6L","Chr6S",
                                  "Chr7L","Chr7S","Chr8L","Chr8S",
                                  "Chr9_10L","Chr9_10S"), 
                          label = c("","","X. laevis","",
                                    "","","","",
                                    "","","","",
                                    "","","","",
                                    "",""))
# lend
lend_labels <- data.frame(CHR = c("Chr1L","Chr1S","Chr2L","Chr2S",
                                  "Chr3L","Chr3S","Chr4L","Chr4S",
                                  "Chr5L","Chr5S","Chr6L","Chr6S",
                                  "Chr7L","Chr7S","Chr8L","Chr8S",
                                  "Chr9_10L","Chr9_10S"), 
                          label = c("","","","",
                                    "X. lenduensis","","","",
                                    "","","","",
                                    "","","","",
                                    "",""))

# muel
muel_labels <- data.frame(CHR = c("Chr1L","Chr1S","Chr2L","Chr2S",
                                  "Chr3L","Chr3S","Chr4L","Chr4S",
                                  "Chr5L","Chr5S","Chr6L","Chr6S",
                                  "Chr7L","Chr7S","Chr8L","Chr8S",
                                  "Chr9_10L","Chr9_10S"), 
                          label = c("","","","",
                                    "","","X. muelleri","",
                                    "","","","",
                                    "","","","",
                                    "",""))  
# pygm
pygm_labels <- data.frame(CHR = c("Chr1L","Chr1S","Chr2L","Chr2S",
                                  "Chr3L","Chr3S","Chr4L","Chr4S",
                                  "Chr5L","Chr5S","Chr6L","Chr6S",
                                  "Chr7L","Chr7S","Chr8L","Chr8S",
                                  "Chr9_10L","Chr9_10S"), 
                          label = c("","","","",
                                    "","","","",
                                    "","","","",
                                    "","","X. pygmaeus","",
                                    "",""))    


jpeg("./SI_FigX_gene_density_SLregions.jpg",w=7, h=8.0, units ="in", bg="transparent", res = 200)
  ggplot(my_gene_startsites_chronly, aes(x=MB)) + 
    #ggplot(my_data, aes(FW_H, group=group, col=group)) + 
    geom_density(alpha = 0.5, adjust=0.1) +
    facet_wrap(~CHR, ncol=2, scales = "free_x") +
    xlab("Position (Mb)") + ylab("Density") +
    scale_y_continuous(breaks=seq(0,0.05, 0.05)) +
    # allofraseri
    geom_area(data = ~subset(., CHR == "Chr7L"), aes(x = stage(MB, 
      after_stat = oob_censor(x, c(17.2, 19.8)))),
      stat = "density", adjust=0.1, col="red", fill="red") +
    geom_text(x = 20, y = 0.04, aes(label = label), fontface = "italic", col = "red",data = allo_labels) +
    # fischbergi
    geom_area(data = ~subset(., CHR == "Chr3L"), aes(x = stage(MB, 
              after_stat = oob_censor(x, c(41, 105)))),
              stat = "density", adjust=0.1, col="red", fill="red") +
    geom_text(x = 75, y = 0.04, aes(label = label), fontface = "italic", col = "red",data = fisc_labels) +
    # muelleri
    geom_area(data = ~subset(., CHR == "Chr4L"), aes(x = stage(MB, 
              after_stat = oob_censor(x, c(111, 147)))),
              stat = "density", adjust=0.1, col="red", fill="red") +
    geom_text(x = 125, y = 0.04, aes(label = label), fontface = "italic", col = "red",data = muel_labels) +
    # pygmaeus
    geom_area(data = ~subset(., CHR == "Chr8L"), aes(x = stage(MB, 
              after_stat = oob_censor(x, c(130, 133)))),
              stat = "density", adjust=0.1, col="red", fill="red") +
    geom_text(x = 120, y = 0.04, aes(label = label), fontface = "italic", col = "red",data = pygm_labels) +
    # lenduensis
    geom_area(data = ~subset(., CHR == "Chr3L"), aes(x = stage(MB, 
              after_stat = oob_censor(x, c(15, 16.2)))),
              stat = "density", adjust=0.1, col="red", fill="red") +
    geom_text(x = 20, y = 0.04, aes(label = label), fontface = "italic", col = "red",data = lend_labels) +
    # borealis east
    geom_area(data = ~subset(., CHR == "Chr8L"), aes(x = stage(MB, 
              after_stat = oob_censor(x, c(0, 54)))),
              stat = "density", adjust=0.1, col="red", fill="red") +
    geom_text(x = 25, y = 0.04, aes(label = label), fontface = "italic", col = "red",data = bore_labels) +
    # laevis
    geom_area(data = ~subset(., CHR == "Chr2L"), aes(x = stage(MB, 
              after_stat = oob_censor(x, c(181.6, 182.7)))),
              stat = "density", adjust=0.1, col="red", fill="red") +
    geom_text(x = 180, y = 0.04, aes(label = label), fontface = "italic", col = "red", data = laev_labels) +
    theme_classic()
dev.off()  



```
