# Gene density

First harvest the start position of all annotated "genes" like this:
```
zcat < XENTR_10.0_Xenbase_longest.gff3.gz | grep 'gene        ' | cut -f1,4 > XENTR_10.0_Xenbase_longest.gff3_gene_starts.txt zcat < XENLA_10.1_Xenbase_longest.gff3.gz | grep 'gene        ' | cut -f1,4 > XENLA_10.1_Xenbase_longest.gff3_gene_starts.txt 
```
This resulted in 28773 for XENTR_10.0 and 44438 for XENLA_10.1 on chromosomes (not including mtDNA and scaffolds in either genome). For XL, this was 26020 from the L subgenome and 18418 from the S subgenome. Perfect!

An interesting question is whether sex-linked regions are found in areas that have high gene density. To explore this I plotted a faceted histogram of start sites of annotated genes in XL and the locations of sex-linked regions. This is better than a density plot because the chromosomes are comparable across facets with a histogram.

Here is the R code for the XL plot:
```
library(ggplot2)
library(readr)
library(dplyr)
library(scales)
library(tidyverse)

setwd("/Users/Shared/Previously Relocated Items/Security/projects/2022_GBS_lotsof_Xennies/gene_density")

# load the data for XL
my_gene_startsites <- read.table("XENLA_10.1_Xenbase_longest.gff3_gene_starts.txt", header = F, sep = "\t")
unique(my_gene_startsites$V1)
colnames(my_gene_startsites) <- c("CHR","POS")
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

my_gene_startsites_chronly$MB <- as.numeric(my_gene_startsites_chronly$POS)/1000000

my_gene_startsites_chronly$CHR <- factor(my_gene_startsites_chronly$CHR,
                    levels = c("Chr1L","Chr1S","Chr2L","Chr2S",
                               "Chr3L","Chr3S","Chr4L","Chr4S",
                               "Chr5L","Chr5S","Chr6L","Chr6S",
                               "Chr7L","Chr7S","Chr8L","Chr8S",
                               "Chr9_10L","Chr9_10S"), ordered = T)

  
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
                                      "",""),
                            COLOR = c("red","red","red","red",
                                    "red","red","red","red",
                                    "red","red","red","red",
                                    "red","red","red","red",
                                    "red","red"))   
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
                                      "",""),
                          COLOR = c("red","red","red","red",
                                    "red","red","red","red",
                                    "red","red","red","red",
                                    "red","red","red","red",
                                    "red","red"))    
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
                                     "",""),
                          COLOR = c("red","red","red","red",
                                    "red","red","red","red",
                                    "red","red","red","red",
                                    "red","red","red","red",
                                    "red","red"))

# laev
laev_labels <- data.frame(CHR = c("Chr1L","Chr1S","Chr2L","Chr2S",
                                  "Chr3L","Chr3S","Chr4L","Chr4S",
                                  "Chr5L","Chr5S","Chr6L","Chr6S",
                                  "Chr7L","Chr7S","Chr8L","Chr8S",
                                  "Chr9_10L","Chr9_10S"), 
                          label = c("","","X. laevis, X. gilli","",
                                    "","","","",
                                    "","","","",
                                    "","","","",
                                    "",""),
                          COLOR = c("red","red","red","red",
                                    "red","red","red","red",
                                    "red","red","red","red",
                                    "red","red","red","red",
                                    "red","red"))
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
                                    "",""),
                          COLOR = c("red","red","red","red",
                                    "red","red","red","red",
                                    "red","red","red","red",
                                    "red","red","red","red",
                                    "red","red"))

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
                                    "",""),
                          COLOR = c("red","red","red","red",
                                    "red","red","red","red",
                                    "red","red","red","red",
                                    "red","red","red","red",
                                    "red","red"))  
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
                                    "",""),
                          COLOR = c("red","red","red","red",
                                    "red","red","red","red",
                                    "red","red","red","red",
                                    "red","red","red","red",
                                    "red","red"))    


dummy_data <- data.frame (CHR = c("Chr1L","Chr1S","Chr2L","Chr2S",
                                     "Chr3L","Chr3S","Chr4L","Chr4S",
                                     "Chr5L","Chr5S","Chr6L","Chr6S",
                                     "Chr7L","Chr7S","Chr8L","Chr8S",
                                     "Chr9_10L","Chr9_10S"),
                          x = rep(234,18),
                          y = rep(0,18),
                          COLOR = c("white","white","white","white",
                                    "white","white","white","white",
                                    "white","white","white","white",
                                    "white","white","white","white",
                                    "white","white"))


# now make a color variable
my_gene_startsites_chronly$COLOR <- "gray"

for(i in 1:dim(my_gene_startsites_chronly)[1]) {             
  if( # beginif
      # laevis
      (my_gene_startsites_chronly[i,]$CHR == "Chr2L")&&
      (my_gene_startsites_chronly[i,]$MB > 181.6)&&
      (my_gene_startsites_chronly[i,]$MB < 182.7)
      ||
      # borealis
      (my_gene_startsites_chronly[i,]$CHR == "Chr8L")&&
      (my_gene_startsites_chronly[i,]$MB > 0)&&
      (my_gene_startsites_chronly[i,]$MB < 54)
      ||
      # allofraseri
      (my_gene_startsites_chronly[i,]$CHR == "Chr7L")&&
      (my_gene_startsites_chronly[i,]$MB > 18.1)&&
      (my_gene_startsites_chronly[i,]$MB < 19.8)
      ||
      # fischbergi
      (my_gene_startsites_chronly[i,]$CHR == "Chr3L")&&
      (my_gene_startsites_chronly[i,]$MB > 41)&&
      (my_gene_startsites_chronly[i,]$MB < 105)
      ||
      # mulleri
      (my_gene_startsites_chronly[i,]$CHR == "Chr4L")&&
      (my_gene_startsites_chronly[i,]$MB > 111)&&
      (my_gene_startsites_chronly[i,]$MB < 147)
      ||
      # pygmaeus
      (my_gene_startsites_chronly[i,]$CHR == "Chr8L")&&
      (my_gene_startsites_chronly[i,]$MB > 130)&&
      (my_gene_startsites_chronly[i,]$MB < 133)
      ||
      # lenduensis
      (my_gene_startsites_chronly[i,]$CHR == "Chr3L")&&
      (my_gene_startsites_chronly[i,]$MB > 15)&&
      (my_gene_startsites_chronly[i,]$MB < 16.2)
    ) # endif
    {
                my_gene_startsites_chronly[i,]$COLOR <- "red"
    }
}

# count genes in each region
# allofraseri
count <- 0
for(i in 1:dim(my_gene_startsites_chronly)[1]) {             
  if( # beginif
    (my_gene_startsites_chronly[i,]$CHR == "Chr7L")&&
    (my_gene_startsites_chronly[i,]$MB > 18.1)&&
    (my_gene_startsites_chronly[i,]$MB < 19.8)
  ){count <- count +1}
}  
count
# borealis
count <- 0
for(i in 1:dim(my_gene_startsites_chronly)[1]) {             
  if( # beginif
    (my_gene_startsites_chronly[i,]$CHR == "Chr8L")&&
    (my_gene_startsites_chronly[i,]$MB > 0)&&
    (my_gene_startsites_chronly[i,]$MB < 54)
  ){count <- count +1}
}  
count

# fischbergi
count <- 0
for(i in 1:dim(my_gene_startsites_chronly)[1]) {             
  if( # beginif
    (my_gene_startsites_chronly[i,]$CHR == "Chr3L")&&
    (my_gene_startsites_chronly[i,]$MB > 41)&&
    (my_gene_startsites_chronly[i,]$MB < 105)
  ){count <- count +1}
}  
count

# laevis
count <- 0
for(i in 1:dim(my_gene_startsites_chronly)[1]) {             
  if( # beginif
    (my_gene_startsites_chronly[i,]$CHR == "Chr2L")&&
    (my_gene_startsites_chronly[i,]$MB > 181.6)&&
    (my_gene_startsites_chronly[i,]$MB < 182.7)
  ){count <- count +1}
}  
count

# lenduensis
count <- 0
for(i in 1:dim(my_gene_startsites_chronly)[1]) {             
  if( # beginif
    (my_gene_startsites_chronly[i,]$CHR == "Chr3L")&&
    (my_gene_startsites_chronly[i,]$MB > 15)&&
    (my_gene_startsites_chronly[i,]$MB < 16.2)
  ){count <- count +1}
}  
count

# muelleri
count <- 0
for(i in 1:dim(my_gene_startsites_chronly)[1]) {             
  if( # beginif
    (my_gene_startsites_chronly[i,]$CHR == "Chr4L")&&
    (my_gene_startsites_chronly[i,]$MB > 111)&&
    (my_gene_startsites_chronly[i,]$MB < 147)
  ){count <- count +1}
}  
count

# pygmaeus
count <- 0
for(i in 1:dim(my_gene_startsites_chronly)[1]) {             
  if( # beginif
    (my_gene_startsites_chronly[i,]$CHR == "Chr8L")&&
    (my_gene_startsites_chronly[i,]$MB > 130)&&
    (my_gene_startsites_chronly[i,]$MB < 133)
  ){count <- count +1}
}  
count

#jpeg("./SI_FigX_gene_density_SLregions.jpg",w=7, h=8.0, units ="in", bg="transparent", res = 400)
pdf("./SI_FigX_gene_density_SLregions.pdf",w=7, h=8.0, bg="transparent")
  ggplot(my_gene_startsites_chronly, aes(x = MB, fill = COLOR, col = COLOR)) + 
    #ggplot(my_data, aes(FW_H, group=group, col=group)) + 
    #stat_density(alpha = 0.5, adjust=0.1, position = "identity") +
    facet_wrap(~CHR, ncol=2, scales = "free_x") +
    geom_histogram(binwidth=1, pad=TRUE) +
    scale_color_manual(values = c("gray" , "red", "white"))+
    xlab("Position (Mb)") + ylab("log10(Count)") +
    # the expand = c(0,0) removes lines at zero; looks nicer!
    # scale_y_continuous(breaks=seq(0,350, 100), expand = c(0,0), limits = c(0,125)) +
    scale_y_continuous(breaks=c(25,100,500), expand = c(0,0), trans = "log10") + #limits = c(0,200)) +
    # allofraseri
    geom_text(x = 40, y = 2.2, aes(label = label), fontface = "italic", col = "red",data = allo_labels) +
    # fischbergi
    geom_text(x = 130, y = 2.2, aes(label = label), fontface = "italic", col = "red",data = fisc_labels) +
    # muelleri
    geom_text(x = 140, y = 2.2, aes(label = label), fontface = "italic", col = "red",data = muel_labels) +
    # pygmaeus
    geom_text(x = 120, y = 2.2, aes(label = label), fontface = "italic", col = "red",data = pygm_labels) +
    # lenduensis
    geom_text(x = 40, y = 2.2, aes(label = label), fontface = "italic", col = "red",data = lend_labels) +
    # borealis east
    geom_text(x = 35, y = 2.2, aes(label = label), fontface = "italic", col = "red",data = bore_labels) +
    # laevis
    geom_text(x = 130, y = 2.2, aes(label = label), fontface = "italic", col = "red", data = laev_labels) +
    # now overlay points to make the x-axis the way I want it
    geom_point(data=dummy_data, aes(x,y), alpha = 0) +
    theme_classic()+
    theme(strip.background = element_blank(), legend.position="none")
dev.off()  

# no labels
pdf("./SI_FigX_gene_density_SLregions_nolabeles.pdf",w=7, h=8.0, bg="transparent")
ggplot(my_gene_startsites_chronly, aes(x = MB, fill = COLOR, col = COLOR)) + 
  #ggplot(my_data, aes(FW_H, group=group, col=group)) + 
  #stat_density(alpha = 0.5, adjust=0.1, position = "identity") +
  facet_wrap(~CHR, ncol=2, scales = "free_x") +
  geom_histogram(binwidth=1, pad=TRUE) +
  scale_color_manual(values = c("gray" , "red", "white"))+
  xlab("Position (Mb)") + ylab("log10(Count)") +
  # the expand = c(0,0) removes lines at zero; looks nicer!
  # scale_y_continuous(breaks=seq(0,350, 100), expand = c(0,0), limits = c(0,125)) +
  scale_y_continuous(breaks=c(25,100,500), expand = c(0,0), trans = "log10") + #limits = c(0,200)) +
  # allofraseri
  geom_point(data=dummy_data, aes(x,y), alpha = 0) +
  theme_classic()+
  theme(strip.background = element_blank(), legend.position="none")
dev.off()  



# Now do XT (only one region)
# load the data for XL
my_gene_startsites_XT <- read.table("XENTR_10.0_Xenbase_longest.gff3_gene_starts.txt", header = F, sep = "\t")
unique(my_gene_startsites_XT$V1)
colnames(my_gene_startsites_XT) <- c("CHR","POS")
# subset to only include large chrs
my_gene_startsites_chronly_XT <- my_gene_startsites_XT[(my_gene_startsites_XT$CHR == "Chr1")|
                                                   (my_gene_startsites_XT$CHR == "Chr2")|
                                                   (my_gene_startsites_XT$CHR == "Chr3")|
                                                   (my_gene_startsites_XT$CHR == "Chr4")|
                                                   (my_gene_startsites_XT$CHR == "Chr5")|
                                                   (my_gene_startsites_XT$CHR == "Chr6")|
                                                   (my_gene_startsites_XT$CHR == "Chr7")|
                                                   (my_gene_startsites_XT$CHR == "Chr8")|
                                                   (my_gene_startsites_XT$CHR == "Chr9")|
                                                   (my_gene_startsites_XT$CHR == "Chr10"),]

my_gene_startsites_chronly_XT$MB <- as.numeric(my_gene_startsites_chronly_XT$POS)/1000000

my_gene_startsites_chronly_XT$CHR <- factor(my_gene_startsites_chronly_XT$CHR,
       levels = c("Chr1","Chr2","Chr3","Chr4",
                  "Chr5","Chr6","Chr7","Chr8",
                  "Chr9","Chr10"), ordered = T)

my_gene_startsites_chronly_XT$chr_order <- 1
my_gene_startsites_chronly_XT[(my_gene_startsites_chronly_XT$CHR == "Chr2"),]$chr_order <-2
my_gene_startsites_chronly_XT[(my_gene_startsites_chronly_XT$CHR == "Chr3"),]$chr_order <-3
my_gene_startsites_chronly_XT[(my_gene_startsites_chronly_XT$CHR == "Chr4"),]$chr_order <-4
my_gene_startsites_chronly_XT[(my_gene_startsites_chronly_XT$CHR == "Chr5"),]$chr_order <-5
my_gene_startsites_chronly_XT[(my_gene_startsites_chronly_XT$CHR == "Chr6"),]$chr_order <-6
my_gene_startsites_chronly_XT[(my_gene_startsites_chronly_XT$CHR == "Chr7"),]$chr_order <-7
my_gene_startsites_chronly_XT[(my_gene_startsites_chronly_XT$CHR == "Chr8"),]$chr_order <-8
my_gene_startsites_chronly_XT[(my_gene_startsites_chronly_XT$CHR == "Chr9"),]$chr_order <-9
my_gene_startsites_chronly_XT[(my_gene_startsites_chronly_XT$CHR == "Chr10"),]$chr_order <-10

# trop
trop_labels <- data.frame(CHR = c("Chr1","Chr2","Chr3","Chr4",
                                  "Chr5","Chr6","Chr7","Chr8",
                                  "Chr9","Chr10"), 
                          label = c("","","","",
                                    "","","X. tropicalis, X. mellotropicalis","",
                                    "",""),
                          COLOR = c("red","red","red","red",
                                    "red","red","red","red",
                                    "red","red"))   


dummy_data <- data.frame (CHR = c("Chr1","Chr2","Chr3","Chr4",
                                  "Chr5","Chr6","Chr7","Chr8",
                                  "Chr9","Chr10"),
                          x = rep(215,10),
                          y = rep(0,10),
                          COLOR = c("white","white","white","white",
                                    "white","white","white","white",
                                    "white","white"))


# now make a color variable
my_gene_startsites_chronly_XT$COLOR <- "gray"

for(i in 1:dim(my_gene_startsites_chronly_XT)[1]) {             
  if( # beginif
    # laevis
    (my_gene_startsites_chronly_XT[i,]$CHR == "Chr7")&&
    (my_gene_startsites_chronly_XT[i,]$MB > 8)&&
    (my_gene_startsites_chronly_XT[i,]$MB < 12)
  ) # endif
  {
    my_gene_startsites_chronly_XT[i,]$COLOR <- "red"
  }
}

# tropicalis
count <- 0
for(i in 1:dim(my_gene_startsites_chronly_XT)[1]) {             
  if( # beginif
    (my_gene_startsites_chronly_XT[i,]$CHR == "Chr7")&&
    (my_gene_startsites_chronly_XT[i,]$MB > 8)&&
    (my_gene_startsites_chronly_XT[i,]$MB < 12)
  ){count <- count +1}
}  
count

#jpeg("./SI_FigX_gene_density_SLregions_XT.jpg",w=7, h=5.0, units ="in", bg="transparent", res = 400)
pdf("./SI_FigX_gene_density_SLregions_XT.pdf",w=7, h=5.0, bg="transparent")
    ggplot(my_gene_startsites_chronly_XT, aes(x = MB, fill = COLOR, col = COLOR)) + 
    #ggplot(my_data, aes(FW_H, group=group, col=group)) + 
    #stat_density(alpha = 0.5, adjust=0.1, position = "identity") +
    facet_wrap(~fct_reorder(CHR,chr_order), ncol=2, scales = "free_x") +
    geom_histogram(binwidth=1) +
    scale_color_manual(values = c("gray" , "red", "white"))+
    xlab("Position (Mb)") + ylab("log10(Count)") +
    # the expand = c(0,0) removes lines at zero; looks nicer!
    scale_y_continuous(breaks=c(25,100,400), expand = c(0,0), trans = "log10") + #limits = c(0,200)) +
    # trop (this label is broken)
    # geom_text(x = 100, y = 2, aes(label = label), fontface = "italic", col = "red",data = trop_labels) +
    # now overlay points to make the x-axis the way I want it
    geom_point(data=dummy_data, aes(x,y), alpha = 0) +
    theme_classic()+
    theme(strip.background = element_blank(), legend.position="none")
dev.off()  
```
