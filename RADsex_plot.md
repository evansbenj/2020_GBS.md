# Plotting RADsex

```R
library(ggplot2)
library(ggpubr)
library(gridExtra)
options(scipen=999)

# laevis ----
setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/2023_RADsex/laev/XL_lab_fam")

my_df <- read.table("XL_family_dist.tsv", header = T)
my_df$groups <- cut(my_df$Markers,               # Add group column
                    breaks = c(-1,0, 1, 5, 10, 20, 35, 100, 500,10000000000),
                    labels = c("0","1", "1-4", "5-9", "10-20", "21-35","36-100","101-500",">500"))

#colors <- colorRampPalette(c("white","blue", "yellow", "red"))(8)
#colors <-paletteer_d("grDevices::blues9")
#colors <- colorRampPalette(c("white","snow2","blue"))(8)
#colors <-paletteer_d("RColorBrewer::Greys")
#colors <-c("#FFFFFFFF", "#F0F0F0FF", "#D9D9D9FF", "#BDBDBDFF", "#969696FF", "#737373FF", "#525252FF", "#252525FF", "#000000FF","black")
colors <-c("#FFFFFFFF","#F4FAFCFF", "#DEE8EBFF", "#C8D6DBFF", "#B2C4CAFF", "#9DB2B9FF", "#87A0A8FF", "#718E98FF", "#5B7C87FF", "#456A76FF")

#library(paletteer)
#library(ggthemes)
#colors <-paletteer_d("grDevices::blues9")
max_n <- min(max(my_df$M),max(my_df$F))+0.5;max_n

# plot
XLfam_RADsex <- ggplot(my_df, aes(M, F, fill= groups)) + 
  geom_tile()+
  scale_fill_manual(breaks = levels(my_df$groups),
                    values = colors) +
  annotate("text", x=15, y=35, label=as.character(expression(italic("X. laevis"))), size=6, parse = T) +
  geom_segment(
    x = 4, y = 17,
    xend = 2, yend = 12,
    lineend = "round", # See available arrow types in example above
    linejoin = "round",
    linewidth = 1, 
    arrow = arrow(length = unit(0.2, "cm")),
    colour = "red") + 
  coord_fixed(ratio=1, xlim=c(-1,45), ylim=c(-1,45), expand=F) +
  #coord_fixed(ratio=1, xlim=c(0,max(my_df$M)), ylim=c(0,max(my_df$F)), expand=F) +
  theme_classic()+
  theme(legend.position= c(0.8, 0.5), legend.title = element_blank()) +
  theme(plot.margin=unit(c(0,0,0,0), 'cm'))+ 
  labs(x = " ") +
  theme(axis.text.x = element_blank())
XLfam_RADsex
# add border to significant tiles
Significant <- my_df[(my_df$Signif == "True")&(my_df$Markers > 0),]

dat <- merge(my_df, Significant)

XLfam_RADsex <- XLfam_RADsex + geom_tile(data=dat, aes(M, F), fill="transparent", colour="red", linewidth=0.25)

XLfam_RADsex

# allo ----
setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/2023_RADsex/allo")

my_df <- read.table("2023_allo_distribution.tsv", header = T)
my_df$groups <- cut(my_df$Markers,               # Add group column
                    breaks = c(-1,0, 1, 5, 10, 20, 35, 100, 500,10000000000),
                    labels = c("0","1", "1-4", "5-9", "10-20", "21-35","36-100","101-500",">500"))

#colors <- colorRampPalette(c("white","blue", "yellow", "red"))(8)
#colors <-paletteer_d("grDevices::blues9")
#colors <- colorRampPalette(c("white","snow2","blue"))(8)
#colors <-paletteer_d("RColorBrewer::Greys")
#colors <-c("#FFFFFFFF", "#F0F0F0FF", "#D9D9D9FF", "#BDBDBDFF", "#969696FF", "#737373FF", "#525252FF", "#252525FF", "#000000FF","black")
colors <-c("#FFFFFFFF","#F4FAFCFF", "#DEE8EBFF", "#C8D6DBFF", "#B2C4CAFF", "#9DB2B9FF", "#87A0A8FF", "#718E98FF", "#5B7C87FF", "#456A76FF")

max_n <- min(max(my_df$M),max(my_df$F))+0.5;max_n

# plot
allo_RADsex <- ggplot(my_df, aes(M, F, fill= groups)) + 
  geom_tile()+
  scale_fill_manual(breaks = levels(my_df$groups),
                    values = colors) +
  annotate("text", x=30, y=35, label=as.character(expression(italic("X. allofraseri"))), size=6, parse = T) +
  geom_segment(
    x = 40, y = 3,
    xend = 37, yend = 3,
    lineend = "round", # See available arrow types in example above
    linejoin = "round",
    linewidth = 1, 
    arrow = arrow(length = unit(0.2, "cm")),
    colour = "red") + 
  coord_fixed(ratio=1, xlim=c(-1,45), ylim=c(-1,45), expand=F) +
  #coord_fixed(ratio=1, xlim=c(0,max(my_df$M)), ylim=c(0,max(my_df$F)), expand=F) +
  theme_classic()+
  theme(legend.position="none") + 
  labs(x = " ") +
  theme(axis.text.x = element_blank())

# add border to significant tiles
Significant <- my_df[(my_df$Signif == "True")&(my_df$Markers > 0),]

dat <- merge(my_df, Significant)

allo_RADsex <- allo_RADsex + geom_tile(data=dat, aes(M, F), fill="transparent", colour="red", linewidth=0.25) 

allo_RADsex


# fisch ----
setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/2023_RADsex/fisc")

my_df <- read.table("2023_fisc_distribution_wo_Z23798.tsv", header = T)
my_df$groups <- cut(my_df$Markers,               # Add group column
                    breaks = c(-1,0, 1, 5, 10, 20, 35, 100, 500,10000000000),
                    labels = c("0","1", "1-4", "5-9", "10-20", "21-35","36-100","101-500",">500"))

#colors <- colorRampPalette(c("white","blue", "yellow", "red"))(8)
#colors <-paletteer_d("grDevices::blues9")
#colors <- colorRampPalette(c("white","snow2","blue"))(8)
#colors <-paletteer_d("RColorBrewer::Greys")
#colors <-c("#FFFFFFFF", "#F0F0F0FF", "#D9D9D9FF", "#BDBDBDFF", "#969696FF", "#737373FF", "#525252FF", "#252525FF", "#000000FF","black")
colors <-c("#FFFFFFFF","#F4FAFCFF", "#DEE8EBFF", "#C8D6DBFF", "#B2C4CAFF", "#9DB2B9FF", "#87A0A8FF", "#718E98FF", "#5B7C87FF", "#456A76FF")

max_n <- min(max(my_df$M),max(my_df$F))+0.5;max_n

# plot
fisc_RADsex <- ggplot(my_df, aes(M, F, fill= groups)) + 
  geom_tile()+
  scale_fill_manual(breaks = levels(my_df$groups),
                    values = colors) +
  annotate("text", x=25, y=35, label=as.character(expression(italic("X. fischbergi"))), size=6, parse = T) +
  geom_segment(
    x = 4, y = 27,
    xend = 2, yend = 22,
    lineend = "round", # See available arrow types in example above
    linejoin = "round",
    linewidth = 1, 
    arrow = arrow(length = unit(0.2, "cm")),
    colour = "red") + 
  coord_fixed(ratio=1, xlim=c(-1,45), ylim=c(-1,45), expand=F) +
  #coord_fixed(ratio=1, xlim=c(0,max(my_df$M)), ylim=c(0,max(my_df$F)), expand=F) +
  theme_classic()+
  theme(legend.position="none") + 
  labs(x = " ") +
  theme(axis.text.x = element_blank())
# add border to significant tiles
Significant <- my_df[(my_df$Signif == "True")&(my_df$Markers > 0),]

dat <- merge(my_df, Significant)

fisc_RADsex <- fisc_RADsex + geom_tile(data=dat, aes(M, F), fill="transparent", colour="red", linewidth=0.25)

fisc_RADsex
#ggsave("2023_wild_bore_eastonly_Rplot.pdf", width = 5, height = 4, dpi = 200)



# bore ----
setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/2023_RADsex/bor/lab_borealis")

my_df <- read.table("2023_bore_distribution.tsv", header = T)
my_df$groups <- cut(my_df$Markers,               # Add group column
                    breaks = c(-1,0, 1, 5, 10, 20, 35, 100, 500,10000000000),
                    labels = c("0","1", "1-4", "5-9", "10-20", "21-35","36-100","101-500",">500"))

#colors <- colorRampPalette(c("white","blue", "yellow", "red"))(8)
#colors <-paletteer_d("grDevices::blues9")
#colors <- colorRampPalette(c("white","snow2","blue"))(8)
#colors <-paletteer_d("RColorBrewer::Greys")
#colors <-c("#FFFFFFFF", "#F0F0F0FF", "#D9D9D9FF", "#BDBDBDFF", "#969696FF", "#737373FF", "#525252FF", "#252525FF", "#000000FF","black")
colors <-c("#FFFFFFFF","#F4FAFCFF", "#DEE8EBFF", "#C8D6DBFF", "#B2C4CAFF", "#9DB2B9FF", "#87A0A8FF", "#718E98FF", "#5B7C87FF", "#456A76FF")

max_n <- min(max(my_df$M),max(my_df$F))+0.5;max_n

# plot
bore_RADsex <- ggplot(my_df, aes(M, F, fill= groups)) + 
  geom_tile()+
  scale_fill_manual(breaks = levels(my_df$groups),
                    values = colors) +
  annotate("text", x=25, y=35, label=as.character(expression(italic("X. borealis"))), size=6, parse = T) +
  geom_segment(
    x = 4, y = 31,
    xend = 2, yend = 26,
    lineend = "round", # See available arrow types in example above
    linejoin = "round",
    linewidth = 1, 
    arrow = arrow(length = unit(0.2, "cm")),
    colour = "red") + 
  coord_fixed(ratio=1, xlim=c(-1,45), ylim=c(-1,45), expand=F) +
  #coord_fixed(ratio=1, xlim=c(0,max(my_df$M)), ylim=c(0,max(my_df$F)), expand=F) +
  theme_classic()+
  theme(legend.position="none") 

# add border to significant tiles
Significant <- my_df[(my_df$Signif == "True")&(my_df$Markers > 0),]

dat <- merge(my_df, Significant)

bore_RADsex <- bore_RADsex + geom_tile(data=dat, aes(M, F), fill="transparent", colour="red", linewidth=0.25)

bore_RADsex


# pygm ----
setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/2023_RADsex/pygm")

my_df <- read.table("2023_pygm_distribution.tsv", header = T)
my_df$groups <- cut(my_df$Markers,               # Add group column
                    breaks = c(-1,0, 1, 5, 10, 20, 35, 100, 500,10000000000),
                    labels = c("0","1", "1-4", "5-9", "10-20", "21-35","36-100","101-500",">500"))

#colors <- colorRampPalette(c("white","blue", "yellow", "red"))(8)
#colors <-paletteer_d("grDevices::blues9")
#colors <- colorRampPalette(c("white","snow2","blue"))(8)
#colors <-paletteer_d("RColorBrewer::Greys")
#colors <-c("#FFFFFFFF", "#F0F0F0FF", "#D9D9D9FF", "#BDBDBDFF", "#969696FF", "#737373FF", "#525252FF", "#252525FF", "#000000FF","black")
colors <-c("#FFFFFFFF","#F4FAFCFF", "#DEE8EBFF", "#C8D6DBFF", "#B2C4CAFF", "#9DB2B9FF", "#87A0A8FF", "#718E98FF", "#5B7C87FF", "#456A76FF")

max_n <- min(max(my_df$M),max(my_df$F))+0.5;max_n

# plot
pygm_RADsex <- ggplot(my_df, aes(M, F, fill= groups)) + 
  geom_tile()+
  scale_fill_manual(breaks = levels(my_df$groups),
                    values = colors) +
  annotate("text", x=30, y=35, label=as.character(expression(italic("X. pygmaeus"))), size=6, parse = T) +
  geom_segment(
    x = 8, y = 43,
    xend = 4, yend = 42,
    lineend = "round", # See available arrow types in example above
    linejoin = "round",
    linewidth = 1, 
    arrow = arrow(length = unit(0.2, "cm")),
    colour = "red") + 
  coord_fixed(ratio=1, xlim=c(-1,45), ylim=c(-1,45), expand=F) +
  #coord_fixed(ratio=1, xlim=c(0,max(my_df$M)), ylim=c(0,max(my_df$F)), expand=F) +
  theme_classic()+
  theme(legend.position="none") + 
  labs(x = " ", y = " ") +
  theme(axis.text.x = element_blank(),axis.text.y =element_blank())
pygm_RADsex
# add border to significant tiles
Significant <- my_df[(my_df$Signif == "True")&(my_df$Markers > 0),]

dat <- merge(my_df, Significant)

pygm_RADsex <- pygm_RADsex + geom_tile(data=dat, aes(M, F), fill="transparent", colour="red", linewidth=0.25)

pygm_RADsex 



# lend ----
setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/2023_RADsex/lend")

my_df <- read.table("2023_lend_distribution.tsv", header = T)
my_df$groups <- cut(my_df$Markers,               # Add group column
                    breaks = c(-1,0, 1, 5, 10, 20, 35, 100, 500,10000000000),
                    labels = c("0","1", "1-4", "5-9", "10-20", "21-35","36-100","101-500",">500"))

#colors <- colorRampPalette(c("white","blue", "yellow", "red"))(8)
#colors <-paletteer_d("grDevices::blues9")
#colors <- colorRampPalette(c("white","snow2","blue"))(8)
#colors <-paletteer_d("RColorBrewer::Greys")
#colors <-c("#FFFFFFFF", "#F0F0F0FF", "#D9D9D9FF", "#BDBDBDFF", "#969696FF", "#737373FF", "#525252FF", "#252525FF", "#000000FF","black")
colors <-c("#FFFFFFFF","#F4FAFCFF", "#DEE8EBFF", "#C8D6DBFF", "#B2C4CAFF", "#9DB2B9FF", "#87A0A8FF", "#718E98FF", "#5B7C87FF", "#456A76FF")

max_n <- min(max(my_df$M),max(my_df$F))+0.5;max_n

# plot
lend_RADsex <- ggplot(my_df, aes(M, F, fill= groups)) + 
  geom_tile()+
  scale_fill_manual(breaks = levels(my_df$groups),
                    values = colors) +
  annotate("text", x=25, y=35, label=as.character(expression(italic("X. lenduensis"))), size=6, parse = T) +
  geom_segment(
    x = 4, y = 32,
    xend = 2, yend = 28,
    lineend = "round", # See available arrow types in example above
    linejoin = "round",
    linewidth = 1, 
    arrow = arrow(length = unit(0.2, "cm")),
    colour = "red") + 
  coord_fixed(ratio=1, xlim=c(-1,45), ylim=c(-1,45), expand=F) +
  #coord_fixed(ratio=1, xlim=c(0,max(my_df$M)), ylim=c(0,max(my_df$F)), expand=F) +
  theme_classic()+
  theme(legend.position="none")+
  labs(x = " ", y = " ") +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank()) 

# add border to significant tiles
Significant <- my_df[(my_df$Signif == "True")&(my_df$Markers > 0),]

dat <- merge(my_df, Significant)

lend_RADsex <- lend_RADsex + geom_tile(data=dat, aes(M, F), fill="transparent", colour="red", linewidth=0.25)

lend_RADsex


# muel ----
setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/2023_RADsex/muel")

my_df <- read.table("2023_muel_distribution.tsv", header = T)
my_df$groups <- cut(my_df$Markers,               # Add group column
                    breaks = c(-1,0, 1, 5, 10, 20, 35, 100, 500,10000000000),
                    labels = c("0","1", "1-4", "5-9", "10-20", "21-35","36-100","101-500",">500"))

#colors <- colorRampPalette(c("white","blue", "yellow", "red"))(8)
#colors <-paletteer_d("grDevices::blues9")
#colors <- colorRampPalette(c("white","snow2","blue"))(8)
#colors <-paletteer_d("RColorBrewer::Greys")
#colors <-c("#FFFFFFFF", "#F0F0F0FF", "#D9D9D9FF", "#BDBDBDFF", "#969696FF", "#737373FF", "#525252FF", "#252525FF", "#000000FF","black")
colors <-c("#FFFFFFFF","#F4FAFCFF", "#DEE8EBFF", "#C8D6DBFF", "#B2C4CAFF", "#9DB2B9FF", "#87A0A8FF", "#718E98FF", "#5B7C87FF", "#456A76FF")

max_n <- min(max(my_df$M),max(my_df$F))+0.5;max_n

# plot
muel_RADsex <- ggplot(my_df, aes(M, F, fill= groups)) + 
  geom_tile()+
  scale_fill_manual(breaks = levels(my_df$groups),
                    values = colors) +
  annotate("text", x=21, y=35, label=as.character(expression(italic("X. muelleri"))), size=6, parse = T) +
  geom_segment(
    x = 4, y =18,
    xend = 2, yend = 14,
    lineend = "round", # See available arrow types in example above
    linejoin = "round",
    linewidth = 1, 
    arrow = arrow(length = unit(0.2, "cm")),
    colour = "red") + 
  coord_fixed(ratio=1, xlim=c(-1,45), ylim=c(-1,45), expand=F) +
  #coord_fixed(ratio=1, xlim=c(0,max(my_df$M)), ylim=c(0,max(my_df$F)), expand=F) +
  theme_classic()+
  theme(legend.position="none")+ 
  labs(x = " ", y = " ") +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank())

# add border to significant tiles
Significant <- my_df[(my_df$Signif == "True")&(my_df$Markers > 0),]

dat <- merge(my_df, Significant)

muel_RADsex <- muel_RADsex + geom_tile(data=dat, aes(M, F), fill="transparent", colour="red", linewidth=0.25)

muel_RADsex


# trop_Mitros_combined ----
setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/2023_RADsex/trop/both_Mitros_combined")

my_df <- read.table("Mitros_bothfam_dist.tsv", header = T)
my_df$groups <- cut(my_df$Markers,               # Add group column
                    breaks = c(-1,0, 1, 5, 10, 20, 35, 100, 500,10000000000),
                    labels = c("0","1", "1-4", "5-9", "10-20", "21-35","36-100","101-500",">500"))

#colors <- colorRampPalette(c("white","blue", "yellow", "red"))(8)
#colors <-paletteer_d("grDevices::blues9")
#colors <- colorRampPalette(c("white","snow2","blue"))(8)
#colors <-paletteer_d("RColorBrewer::Greys")
#colors <-c("#FFFFFFFF", "#F0F0F0FF", "#D9D9D9FF", "#BDBDBDFF", "#969696FF", "#737373FF", "#525252FF", "#252525FF", "#000000FF","black")
colors <-c("#FFFFFFFF","#F4FAFCFF", "#DEE8EBFF", "#C8D6DBFF", "#B2C4CAFF", "#9DB2B9FF", "#87A0A8FF", "#718E98FF", "#5B7C87FF", "#456A76FF")

max_n <- min(max(my_df$M),max(my_df$F))+0.5;max_n

# plot
trop_RADsex <- ggplot(my_df, aes(M, F, fill= groups)) + 
  geom_tile()+
  scale_fill_manual(breaks = levels(my_df$groups),
                    values = colors) +
  annotate("text", x=5, y=38, label=as.character(expression(italic("X. tropicalis"))), size=6, parse = T) +
  geom_segment(
    x = 0, y = 43,
    xend = 0, yend = 40,
    lineend = "round", # See available arrow types in example above
    linejoin = "round",
    linewidth = 1, 
    arrow = arrow(length = unit(0.2, "cm")),
    colour = "red") + 
  coord_fixed(ratio=1, xlim=c(0,48), ylim=c(0,52), expand=F) +
  #coord_fixed(ratio=1, xlim=c(0,max(my_df$M)), ylim=c(0,max(my_df$F)), expand=F) +
  theme_classic()+
  theme(legend.position="none") 

# add border to significant tiles
Significant <- my_df[(my_df$Signif == "True")&(my_df$Markers > 0),]

dat <- merge(my_df, Significant)

muel_RADsex <- muel_RADsex + geom_tile(data=dat, aes(M, F), fill="transparent", colour="red", linewidth=0.25)

muel_RADsex

# add border to significant tiles
Significant <- my_df[(my_df$Signif == "True")&(my_df$Markers > 0),]

dat <- merge(my_df, Significant)

trop_RADsex <- trop_RADsex + geom_tile(data=dat, aes(M, F), fill="transparent", colour="red", linewidth=0.25)

trop_RADsex




# trop_GhE_combined ----
setwd("/Users/Shared/Previously Relocated Items/Security/projects/2022_GBS_lotsof_Xennies/2023_RADsex/trop/GhE")

my_df <- read.table("GhE_dist.tsv", header = T)
my_df$groups <- cut(my_df$Markers,               # Add group column
                    breaks = c(-1,0, 1, 5, 10, 20, 35, 100, 500,10000000000),
                    labels = c("0","1", "1-4", "5-9", "10-20", "21-35","36-100","101-500",">500"))

#colors <- colorRampPalette(c("white","blue", "yellow", "red"))(8)
#colors <-paletteer_d("grDevices::blues9")
#colors <- colorRampPalette(c("white","snow2","blue"))(8)
#colors <-paletteer_d("RColorBrewer::Greys")
#colors <-c("#FFFFFFFF", "#F0F0F0FF", "#D9D9D9FF", "#BDBDBDFF", "#969696FF", "#737373FF", "#525252FF", "#252525FF", "#000000FF","black")
colors <-c("#FFFFFFFF","#F4FAFCFF", "#DEE8EBFF", "#C8D6DBFF", "#B2C4CAFF", "#9DB2B9FF", "#87A0A8FF", "#718E98FF", "#5B7C87FF", "#456A76FF")

max_n <- min(max(my_df$M),max(my_df$F))+0.5;max_n

# plot
trop_GhE_RADsex <- ggplot(my_df, aes(M, F, fill= groups)) + 
  geom_tile()+
  scale_fill_manual(breaks = levels(my_df$groups),
                    values = colors) +
  annotate("text", x=5, y=38, label=as.character(expression(italic("X. tropicalis"))), size=6, parse = T) +
  geom_segment(
    x = 0, y = 43,
    xend = 0, yend = 40,
    lineend = "round", # See available arrow types in example above
    linejoin = "round",
    linewidth = 1, 
    arrow = arrow(length = unit(0.2, "cm")),
    colour = "red") + 
  coord_fixed(ratio=1, xlim=c(-1,45), ylim=c(-1,45), expand=F) +
  #coord_fixed(ratio=1, xlim=c(0,max(my_df$M)), ylim=c(0,max(my_df$F)), expand=F) +
  theme_classic()+
  theme(legend.position="none") 


# add border to significant tiles
Significant <- my_df[(my_df$Signif == "True")&(my_df$Markers > 0),]

dat <- merge(my_df, Significant)

trop_GhE_RADsex <- trop_GhE_RADsex + geom_tile(data=dat, aes(M, F), fill="transparent", colour="red", linewidth=0.25)

trop_GhE_RADsex


# trop_GhW_combined ----
setwd("/Users/Shared/Previously Relocated Items/Security/projects/2022_GBS_lotsof_Xennies/2023_RADsex/trop/GhW")

my_df <- read.table("GhW_dist.tsv", header = T)
my_df$groups <- cut(my_df$Markers,               # Add group column
                    breaks = c(-1,0, 1, 5, 10, 20, 35, 100, 500,10000000000),
                    labels = c("0","1", "1-4", "5-9", "10-20", "21-35","36-100","101-500",">500"))

#colors <- colorRampPalette(c("white","blue", "yellow", "red"))(8)
#colors <-paletteer_d("grDevices::blues9")
#colors <- colorRampPalette(c("white","snow2","blue"))(8)
#colors <-paletteer_d("RColorBrewer::Greys")
#colors <-c("#FFFFFFFF", "#F0F0F0FF", "#D9D9D9FF", "#BDBDBDFF", "#969696FF", "#737373FF", "#525252FF", "#252525FF", "#000000FF","black")
colors <-c("#FFFFFFFF","#F4FAFCFF", "#DEE8EBFF", "#C8D6DBFF", "#B2C4CAFF", "#9DB2B9FF", "#87A0A8FF", "#718E98FF", "#5B7C87FF", "#456A76FF")

max_n <- min(max(my_df$M),max(my_df$F))+0.5;max_n

# plot
trop_GhW_RADsex <- ggplot(my_df, aes(M, F, fill= groups)) + 
  geom_tile()+
  scale_fill_manual(breaks = levels(my_df$groups),
                    values = colors) +
  annotate("text", x=15, y=35, label=as.character(expression(italic("X. tropicalis"))), size=6, parse = T) +
  geom_segment(
    x = 38, y = 3,
    xend = 35, yend = 3,
    lineend = "round", # See available arrow types in example above
    linejoin = "round",
    linewidth = 1, 
    arrow = arrow(length = unit(0.2, "cm")),
    colour = "red") + 
  geom_segment(
    x = 3, y = 29,
    xend = 3, yend = 26,
    lineend = "round", # See available arrow types in example above
    linejoin = "round",
    linewidth = 1, 
    arrow = arrow(length = unit(0.2, "cm")),
    colour = "red") +
  coord_fixed(ratio=1, xlim=c(-1,45), ylim=c(-1,45), expand=F) +
  #coord_fixed(ratio=1, xlim=c(0,max(my_df$M)), ylim=c(0,max(my_df$F)), expand=F) +
  theme_classic()+
  theme(legend.position="none")+ 
  theme(axis.title.y=element_blank(),axis.text.y = element_blank())


# add border to significant tiles
Significant <- my_df[(my_df$Signif == "True")&(my_df$Markers > 0),]

dat <- merge(my_df, Significant)

trop_GhW_RADsex <- trop_GhW_RADsex + geom_tile(data=dat, aes(M, F), fill="transparent", colour="red", linewidth=0.25)

trop_GhW_RADsex


# trop_C659 ----
setwd("/Users/Shared/Previously Relocated Items/Security/projects/2022_GBS_lotsof_Xennies/2023_RADsex/trop/C659")

my_df <- read.table("C659_dist.tsv", header = T)
my_df$groups <- cut(my_df$Markers,               # Add group column
                    breaks = c(-1,0, 1, 5, 10, 20, 35, 100, 500,10000000000),
                    labels = c("0","1", "1-4", "5-9", "10-20", "21-35","36-100","101-500",">500"))

#colors <- colorRampPalette(c("white","blue", "yellow", "red"))(8)
#colors <-paletteer_d("grDevices::blues9")
#colors <- colorRampPalette(c("white","snow2","blue"))(8)
#colors <-paletteer_d("RColorBrewer::Greys")
#colors <-c("#FFFFFFFF", "#F0F0F0FF", "#D9D9D9FF", "#BDBDBDFF", "#969696FF", "#737373FF", "#525252FF", "#252525FF", "#000000FF","black")
colors <-c("#FFFFFFFF","#F4FAFCFF", "#DEE8EBFF", "#C8D6DBFF", "#B2C4CAFF", "#9DB2B9FF", "#87A0A8FF", "#718E98FF", "#5B7C87FF", "#456A76FF")

max_n <- min(max(my_df$M),max(my_df$F))+0.5;max_n

# plot
trop_C659_RADsex <- ggplot(my_df, aes(M, F, fill= groups)) + 
  geom_tile()+
  scale_fill_manual(breaks = levels(my_df$groups),
                    values = colors) +
  #geom_segment(aes(x = -0.5, y = -0.5, xend = max_n, yend = max_n)) +
  theme_classic()

# add border to significant tiles
Significant <- my_df[(my_df$Signif == "True")&(my_df$Markers > 0),]

dat <- merge(my_df, Significant)

trop_C659_RADsex <- trop_C659_RADsex + geom_tile(data=dat, aes(M, F), fill="transparent", colour="red", linewidth=0.25)

trop_C659_RADsex


# trop_C660 ----
setwd("/Users/Shared/Previously Relocated Items/Security/projects/2022_GBS_lotsof_Xennies/2023_RADsex/trop/C660")

my_df <- read.table("C660_dist.tsv", header = T)
my_df$groups <- cut(my_df$Markers,               # Add group column
                    breaks = c(-1,0, 1, 5, 10, 20, 35, 100, 500,10000000000),
                    labels = c("0","1", "1-4", "5-9", "10-20", "21-35","36-100","101-500",">500"))

#colors <- colorRampPalette(c("white","blue", "yellow", "red"))(8)
#colors <-paletteer_d("grDevices::blues9")
#colors <- colorRampPalette(c("white","snow2","blue"))(8)
#colors <-paletteer_d("RColorBrewer::Greys")
#colors <-c("#FFFFFFFF", "#F0F0F0FF", "#D9D9D9FF", "#BDBDBDFF", "#969696FF", "#737373FF", "#525252FF", "#252525FF", "#000000FF","black")
colors <-c("#FFFFFFFF","#F4FAFCFF", "#DEE8EBFF", "#C8D6DBFF", "#B2C4CAFF", "#9DB2B9FF", "#87A0A8FF", "#718E98FF", "#5B7C87FF", "#456A76FF")

max_n <- min(max(my_df$M),max(my_df$F))+0.5;max_n

# plot
trop_C660_RADsex <- ggplot(my_df, aes(M, F, fill= groups)) + 
  geom_tile()+
  scale_fill_manual(breaks = levels(my_df$groups),
                    values = colors) +
  #geom_segment(aes(x = -0.5, y = -0.5, xend = max_n, yend = max_n)) +
  theme_classic()

# add border to significant tiles
Significant <- my_df[(my_df$Signif == "True")&(my_df$Markers > 0),]

dat <- merge(my_df, Significant)

trop_C660_RADsex <- trop_C660_RADsex + geom_tile(data=dat, aes(M, F), fill="transparent", colour="red", linewidth=0.25)

trop_C660_RADsex

# SI Fig 2 ----
setwd("/Users/Shared/Previously Relocated Items/Security/projects/2022_GBS_lotsof_Xennies/2023_RADsex/")
margin = theme(plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm"))
plotz <- list(XLfam_RADsex,pygm_RADsex,
              allo_RADsex,lend_RADsex,
              fisc_RADsex,muel_RADsex,
              bore_RADsex, trop_GhW_RADsex)


jpeg("./SI_Fig2_RADsex_plotz_8_species.jpg",w=6, h=11.0, units ="in", bg="transparent", res = 500)
  grid.arrange(grobs = lapply(plotz, "+", margin), ncol = 2, nrow=4)
  #grid.arrange(XLfam_RADsex,pygm_RADsex,
  #           allo_RADsex,lend_RADsex,
  #           fisc_RADsex,muel_RADsex,
  #           bore_RADsex, trop_GhW_RADsex, ncol = 2, nrow=4)
dev.off()


```
