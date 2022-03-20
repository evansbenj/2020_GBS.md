# Final XL file

In this directory:
```
/home/ben/projects/rrg-ben/ben/2020_GBS_muel_fish_allo_cliv_laev/raw_data/cutaddapted_by_species_across_three_plates/laevis/vcfs_after_filtering_and_removal_with_XG
```

# For PCA, get rid of non-XL samples and divide file up by subgenome:
```
bcftools concat DB_newnew_chr1L_out.vcf_filtered.vcf.gz_filtered_removed.vcf DB_newnew_chr1S_out.vcf_filtered.vcf.gz_filtered_removed.vcf DB_newnew_chr2L_out.vcf_filtered.vcf.gz_filtered_removed.vcf DB_newnew_chr2S_out.vcf_filtered.vcf.gz_filtered_removed.vcf DB_newnew_chr3L_out.vcf_filtered.vcf.gz_filtered_removed.vcf DB_newnew_chr3S_out.vcf_filtered.vcf.gz_filtered_removed.vcf DB_newnew_chr4L_out.vcf_filtered.vcf.gz_filtered_removed.vcf DB_newnew_chr4S_out.vcf_filtered.vcf.gz_filtered_removed.vcf DB_newnew_chr5L_out.vcf_filtered.vcf.gz_filtered_removed.vcf DB_newnew_chr5S_out.vcf_filtered.vcf.gz_filtered_removed.vcf DB_newnew_chr6L_out.vcf_filtered.vcf.gz_filtered_removed.vcf DB_newnew_chr6S_out.vcf_filtered.vcf.gz_filtered_removed.vcf DB_newnew_chr7L_out.vcf_filtered.vcf.gz_filtered_removed.vcf DB_newnew_chr7S_out.vcf_filtered.vcf.gz_filtered_removed.vcf DB_newnew_chr8L_out.vcf_filtered.vcf.gz_filtered_removed.vcf DB_newnew_chr8S_out.vcf_filtered.vcf.gz_filtered_removed.vcf DB_newnew_chr9_10L_out.vcf_filtered.vcf.gz_filtered_removed.vcf DB_newnew_chr9_10S_out.vcf_filtered.vcf.gz_filtered_removed.vcf -o all_XL_withXG.vcf.gz -O z   
```
```
bcftools view -Ou --samples "^RT5_Botsw_GGATTGGT_cuttrim_sorted.bam,XGL713_123.fq.gz_sorted.bam,amnh17260_Nigeria_GTGAGGGT_cuttrim_sorted.bam,BJE2897_Lendu_TTCCTGGA_cuttrim_sorted.bam,BJE3252_Cameroon_TATCGGGA_cuttrim_sorted.bam" all_XL_withXG.vcf.gz -o all_XL_onlyXL.vcf.gz -O z
```
```
bcftools view all_XL_onlyXL.vcf.gz --regions chr1L,chr2L,chr3L,chr4L,chr5L,chr6L,chr7L,chr8L,chr9_10L -o all_XL_onlyXL_only_Lsubgenome.vcf.gz -O z
```
```
bcftools view all_XL_onlyXL.vcf.gz --regions chr1S,chr2S,chr3S,chr4S,chr5S,chr6S,chr7S,chr8S,chr9_10S -o all_XL_onlyXL_only_Ssubgenome.vcf.gz -O z
```
# PCA and plotting
```R
library("devtools")
library(gdsfmt)
library(SNPRelate)
library("ggpubr")

setwd("/Users/Shared/Previously Relocated Items/Security/projects/2020_XL_GBS/SNPRelate")

vcf.fn <- "all_XL_onlyXL_only_Lsubgenome.vcf.gz"
#vcf.fn <- "DB_new_chr9_10L_out.vcf_filtered.vcf.gz_filtered_removed.vcf.gz"
snpgdsVCF2GDS(vcf.fn, "test.gds", method="biallelic.only")

genofile = openfn.gds("test.gds", readonly=FALSE)

samp.annot<-data.frame(pop.group = c("2014_Inhaca_10","2014_Inhaca_150","2014_Inhaca_152",
                                     "2014_Inhaca_24","2014_Inhaca_38","2014_Inhaca_52",
                                     "2014_Inhaca_65","946_Draken","993_Draken","BJE3508_DeDorn","BJE3509_DeDorn",
                                     "BJE3510_DeDorn","BJE3511_DeDorn","BJE3512_DeDorn","BJE3513_DeDorn",
                                     "BJE3514_DeDorn","BJE3515_DeDorn","BJE3525_Laigns","BJE3526_Laigns",
                                     "BJE3527_Laigns","BJE3528_Laigns","BJE3529_Laigns","BJE3530_Laigns",
                                     "BJE3531_Laigns","BJE3532_Laigns","BJE3533_Laigns","BJE3534_BW",
                                     "BJE3535_BW","BJE3536_BW","BJE3537_BW","BJE3538_BW","BJE3539_BW",
                                     "BJE3541_BW","BJE3542_BW","BJE3543_BW","BJE3544_BW","BJE3545_BW",
                                     "BJE3546_BW","BJE3547_GRNP","BJE3548_GRNP","BJE3549_GRNP","BJE3550_GRNP",
                                     "BJE3551_GRNP","BJE3552_GRNP","BJE3553_GRNP","BJE3554_GRNP",
                                     "BJE3573_VicW","BJE3574_VicW","BJE3575_Kimber","BJE3576_Kimber",
                                     "BJE3578_Kimber","BJE3579_Kimber","BJE3580_Kimber",
                                     "BJE3581_Kimber","BJE3582_Kimber","BJE3633_Niewou",
                                     "BJE3640_Niewou","BJE3641_Niewou","BJE3642_Niewou","BJE3644_Niewou",
                                     "BJE3645_Niewou","BJE3647_Niewou","BJE3654_ThreeSis","BJE3655_ThreeSis",
                                     "BJE3656_ThreeSis","BJE3657_ThreeSis","BJE3658_ThreeSis",
                                     "BJE3659_ThreeSis","BJE3660_ThreeSis","BJE3661_ThreeSis",
                                     "BJE3662_ThreeSis","BJE3663_ThreeSis","BJE3664_ThreeSis",
                                     "BJE3665_ThreeSis","BJE3666_ThreeSis","BJE3667_Citrus",
                                     "BJE3668_Citrus","BJE3669_Citrus","BJE3670_Citrus","BJE3671_Citrus",
                                     "BJE3672_Citrus","BJE3673_Citrus","BJE3674_Citrus","BJE3675_Citrus",
                                     "BJE3676_Citrus","BJE3677_Citrus","BJE3678_Citrus","JM_no_label1_Draken",
                                     "JM_no_label2_Draken","Vred_8_Vred"))

add.gdsn(genofile, "sample.annot", samp.annot)

snpgdsSummary("test.gds")

pca <- snpgdsPCA(genofile, num.thread=2)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))


tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)

library(ggplot2)
library(ggrepel)

#pdf("PCA_plot_XL_Lsubgenome.pdf",w=8, h=8, version="1.4", bg="transparent")
tab$Species <- c("laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis")
tab$samp.color <- c("yellow","yellow","yellow","yellow","yellow","yellow","yellow","yellow","yellow","blue","blue","blue","blue","blue","blue","blue","blue","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","blue","blue","blue","blue","blue","blue","blue","blue","yellow","yellow","yellow","yellow","yellow","yellow","yellow","yellow","yellow","purple","purple","purple","purple","purple","purple","purple","yellow","yellow","yellow","yellow","yellow","yellow","yellow","yellow","yellow","yellow","yellow","yellow","yellow","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","yellow","yellow","blue")
tab$samp.fieldid <- c("2014_Inhaca_10","2014_Inhaca_150","2014_Inhaca_152",
                      "2014_Inhaca_24","2014_Inhaca_38","2014_Inhaca_52",
                      "2014_Inhaca_65","946_Draken","993_Draken","BJE3508_DeDorn","BJE3509_DeDorn",
                      "BJE3510_DeDorn","BJE3511_DeDorn","BJE3512_DeDorn","BJE3513_DeDorn",
                      "BJE3514_DeDorn","BJE3515_DeDorn","BJE3525_Laigns","BJE3526_Laigns",
                      "BJE3527_Laigns","BJE3528_Laigns","BJE3529_Laigns","BJE3530_Laigns",
                      "BJE3531_Laigns","BJE3532_Laigns","BJE3533_Laigns","BJE3534_BW",
                      "BJE3535_BW","BJE3536_BW","BJE3537_BW","BJE3538_BW","BJE3539_BW",
                      "BJE3541_BW","BJE3542_BW","BJE3543_BW","BJE3544_BW","BJE3545_BW",
                      "BJE3546_BW","BJE3547_GRNP","BJE3548_GRNP","BJE3549_GRNP","BJE3550_GRNP",
                      "BJE3551_GRNP","BJE3552_GRNP","BJE3553_GRNP","BJE3554_GRNP",
                      "BJE3573_VicW","BJE3574_VicW","BJE3575_Kimber","BJE3576_Kimber",
                      "BJE3578_Kimber","BJE3579_Kimber","BJE3580_Kimber",
                      "BJE3581_Kimber","BJE3582_Kimber","BJE3633_Niewou",
                      "BJE3640_Niewou","BJE3641_Niewou","BJE3642_Niewou","BJE3644_Niewou",
                      "BJE3645_Niewou","BJE3647_Niewou","BJE3654_ThreeSis","BJE3655_ThreeSis",
                      "BJE3656_ThreeSis","BJE3657_ThreeSis","BJE3658_ThreeSis",
                      "BJE3659_ThreeSis","BJE3660_ThreeSis","BJE3661_ThreeSis",
                      "BJE3662_ThreeSis","BJE3663_ThreeSis","BJE3664_ThreeSis",
                      "BJE3665_ThreeSis","BJE3666_ThreeSis","BJE3667_Citrus",
                      "BJE3668_Citrus","BJE3669_Citrus","BJE3670_Citrus","BJE3671_Citrus",
                      "BJE3672_Citrus","BJE3673_Citrus","BJE3674_Citrus","BJE3675_Citrus",
                      "BJE3676_Citrus","BJE3677_Citrus","BJE3678_Citrus","JM_no_label1_Draken",
                      "JM_no_label2_Draken","Vred_8_Vred")

d<-ggplot(data=tab, aes(x=EV2,y=EV1, label = samp.fieldid, color = samp.color)) +
  # label axis 
  labs(x=expression("Eigenvector 2"), y=expression("Eigenvector 1")) +
  # legend details
  scale_colour_manual(name="Population", values = c("yellow"="yellow","blue"="blue","purple"="purple","green"="green"),
                      breaks=c("yellow","blue","purple","green"),
                      labels=c("NW","SE","Admix-NE","Admix-SW"))+
  # add points and fieldID labels
  # geom_text_repel(aes(EV2,EV1, label=(samp.fieldid)), max.overlaps = 1000) + 
  geom_point(size=4) + 
  # change to cleaner theme
  theme_classic(base_size = 16) +
  # make it clean
  theme_bw()+ theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
  # italicize species names
  #theme(legend.text = element_text(face="italic"))+ 
  # move the legend
  theme(legend.position = c(.2, .45)) +
  #theme(legend.position = "none") +
  # add a title
  ggtitle("L subgenome")  
  # remove boxes around legend symbols
  #theme(legend.key = element_blank())
d
#dev.off()


vcf.fn <- "all_XL_onlyXL_only_Ssubgenome.vcf.gz"
#vcf.fn <- "DB_new_chr9_10L_out.vcf_filtered.vcf.gz_filtered_removed.vcf.gz"
snpgdsVCF2GDS(vcf.fn, "test1.gds", method="biallelic.only")

genofile1 = openfn.gds("test1.gds", readonly=FALSE)

samp.annot<-data.frame(pop.group = c("2014_Inhaca_10","2014_Inhaca_150","2014_Inhaca_152",
                                     "2014_Inhaca_24","2014_Inhaca_38","2014_Inhaca_52",
                                     "2014_Inhaca_65","946_Draken","993_Draken","BJE3508_DeDorn","BJE3509_DeDorn",
                                     "BJE3510_DeDorn","BJE3511_DeDorn","BJE3512_DeDorn","BJE3513_DeDorn",
                                     "BJE3514_DeDorn","BJE3515_DeDorn","BJE3525_Laigns","BJE3526_Laigns",
                                     "BJE3527_Laigns","BJE3528_Laigns","BJE3529_Laigns","BJE3530_Laigns",
                                     "BJE3531_Laigns","BJE3532_Laigns","BJE3533_Laigns","BJE3534_BW",
                                     "BJE3535_BW","BJE3536_BW","BJE3537_BW","BJE3538_BW","BJE3539_BW",
                                     "BJE3541_BW","BJE3542_BW","BJE3543_BW","BJE3544_BW","BJE3545_BW",
                                     "BJE3546_BW","BJE3547_GRNP","BJE3548_GRNP","BJE3549_GRNP","BJE3550_GRNP",
                                     "BJE3551_GRNP","BJE3552_GRNP","BJE3553_GRNP","BJE3554_GRNP",
                                     "BJE3573_VicW","BJE3574_VicW","BJE3575_Kimber","BJE3576_Kimber",
                                     "BJE3578_Kimber","BJE3579_Kimber","BJE3580_Kimber",
                                     "BJE3581_Kimber","BJE3582_Kimber","BJE3633_Niewou",
                                     "BJE3640_Niewou","BJE3641_Niewou","BJE3642_Niewou","BJE3644_Niewou",
                                     "BJE3645_Niewou","BJE3647_Niewou","BJE3654_ThreeSis","BJE3655_ThreeSis",
                                     "BJE3656_ThreeSis","BJE3657_ThreeSis","BJE3658_ThreeSis",
                                     "BJE3659_ThreeSis","BJE3660_ThreeSis","BJE3661_ThreeSis",
                                     "BJE3662_ThreeSis","BJE3663_ThreeSis","BJE3664_ThreeSis",
                                     "BJE3665_ThreeSis","BJE3666_ThreeSis","BJE3667_Citrus",
                                     "BJE3668_Citrus","BJE3669_Citrus","BJE3670_Citrus","BJE3671_Citrus",
                                     "BJE3672_Citrus","BJE3673_Citrus","BJE3674_Citrus","BJE3675_Citrus",
                                     "BJE3676_Citrus","BJE3677_Citrus","BJE3678_Citrus","JM_no_label1_Draken",
                                     "JM_no_label2_Draken","Vred_8_Vred"))

add.gdsn(genofile1, "sample.annot", samp.annot)

snpgdsSummary("test1.gds")

pca <- snpgdsPCA(genofile1, num.thread=2)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))


tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)

library(ggplot2)
library(ggrepel)

#pdf("PCA_plot_XL_Lsubgenome.pdf",w=8, h=8, version="1.4", bg="transparent")
tab$Species <- c("laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis")
tab$samp.color <- c("yellow","yellow","yellow","yellow","yellow","yellow","yellow","yellow","yellow","blue","blue","blue","blue","blue","blue","blue","blue","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","blue","blue","blue","blue","blue","blue","blue","blue","yellow","yellow","yellow","yellow","yellow","yellow","yellow","yellow","yellow","purple","purple","purple","purple","purple","purple","purple","yellow","yellow","yellow","yellow","yellow","yellow","yellow","yellow","yellow","yellow","yellow","yellow","yellow","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","yellow","yellow","blue")
tab$samp.fieldid <- c("2014_Inhaca_10","2014_Inhaca_150","2014_Inhaca_152",
                      "2014_Inhaca_24","2014_Inhaca_38","2014_Inhaca_52",
                      "2014_Inhaca_65","946_Draken","993_Draken","BJE3508_DeDorn","BJE3509_DeDorn",
                      "BJE3510_DeDorn","BJE3511_DeDorn","BJE3512_DeDorn","BJE3513_DeDorn",
                      "BJE3514_DeDorn","BJE3515_DeDorn","BJE3525_Laigns","BJE3526_Laigns",
                      "BJE3527_Laigns","BJE3528_Laigns","BJE3529_Laigns","BJE3530_Laigns",
                      "BJE3531_Laigns","BJE3532_Laigns","BJE3533_Laigns","BJE3534_BW",
                      "BJE3535_BW","BJE3536_BW","BJE3537_BW","BJE3538_BW","BJE3539_BW",
                      "BJE3541_BW","BJE3542_BW","BJE3543_BW","BJE3544_BW","BJE3545_BW",
                      "BJE3546_BW","BJE3547_GRNP","BJE3548_GRNP","BJE3549_GRNP","BJE3550_GRNP",
                      "BJE3551_GRNP","BJE3552_GRNP","BJE3553_GRNP","BJE3554_GRNP",
                      "BJE3573_VicW","BJE3574_VicW","BJE3575_Kimber","BJE3576_Kimber",
                      "BJE3578_Kimber","BJE3579_Kimber","BJE3580_Kimber",
                      "BJE3581_Kimber","BJE3582_Kimber","BJE3633_Niewou",
                      "BJE3640_Niewou","BJE3641_Niewou","BJE3642_Niewou","BJE3644_Niewou",
                      "BJE3645_Niewou","BJE3647_Niewou","BJE3654_ThreeSis","BJE3655_ThreeSis",
                      "BJE3656_ThreeSis","BJE3657_ThreeSis","BJE3658_ThreeSis",
                      "BJE3659_ThreeSis","BJE3660_ThreeSis","BJE3661_ThreeSis",
                      "BJE3662_ThreeSis","BJE3663_ThreeSis","BJE3664_ThreeSis",
                      "BJE3665_ThreeSis","BJE3666_ThreeSis","BJE3667_Citrus",
                      "BJE3668_Citrus","BJE3669_Citrus","BJE3670_Citrus","BJE3671_Citrus",
                      "BJE3672_Citrus","BJE3673_Citrus","BJE3674_Citrus","BJE3675_Citrus",
                      "BJE3676_Citrus","BJE3677_Citrus","BJE3678_Citrus","JM_no_label1_Draken",
                      "JM_no_label2_Draken","Vred_8_Vred")

d1<-ggplot(data=tab, aes(x=EV2,y=EV1, label = samp.fieldid, color = samp.color)) +
  # label axis 
  labs(x=expression("Eigenvector 2"), y=expression("Eigenvector 1")) +
  # legend details
  scale_colour_manual(name="Population", values = c("yellow"="yellow","blue"="blue","purple"="purple","green"="green"),
                      breaks=c("yellow","blue","purple","green"),
                      labels=c("NW","SE","Admix-NE","Admix-SW"))+
  # add points and fieldID labels
  #geom_text_repel(aes(EV2,EV1, label=(samp.fieldid)), max.overlaps = 1000) + 
  geom_point(size=4) + 
  # change to cleaner theme
  theme_classic(base_size = 16) +
  # make it clean
  theme_bw()+ theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
  # italicize species names
  #theme(legend.text = element_text(face="italic"))+ 
  # move the legend
  #theme(legend.position = c(.84, .88)) +
  theme(legend.position = "none") +
  # add a title
  ggtitle("S subgenome") + 
  # remove boxes around legend symbols
  theme(legend.key = element_blank())+
  annotate("text", x = 0.02, y = -0.01, label = "Niewoudtville") + #, angle = 45)
  annotate("text", x = -0.08, y = 0.17, label = "Kimberly") +
  annotate("text", x = -0.1, y = 0.11, label = "Inhaca") +
  annotate("text", x = -0.01, y = 0.15, label = "Victoria West") +
  annotate("text", x = 0.07, y = 0.11, label = "Three Sisters") +
  annotate("text", x = 0.15, y = 0, label = "Beaufort West") +
  annotate("text", x = 0.15, y = -0.05, label = "Laignsburg") +
  annotate("text", x = 0, y = -0.1, label = "Garden Route") +
  annotate("text", x = -0.02, y = -0.13, label = "DeDoorns") +
  annotate("text", x = -0.09, y = -0.10, label = "Citrusdale", angle = 45) +
  theme(legend.key = element_blank())
d1
#dev.off()



pdf("PCA_plot_XL_LandS_subgenome.pdf",w=8, h=4, version="1.4", bg="transparent")
ggarrange(d, d1, heights = c(2, 0.7),
          ncol = 2, nrow = 1)
dev.off()











# Dendrogram based on dissimilarity matrix
set.seed(100)
ibs.hc <- snpgdsHCluster(snpgdsIBS(genofile , sample.id=NULL, autosome.only=TRUE, 
                                   remove.monosnp=TRUE,  maf=NaN, missing.rate=NaN, 
                                   num.thread=2, verbose=TRUE))
rv <- snpgdsCutTree(ibs.hc)
plot(rv$dendrogram, leaflab="none")

# https://benbowlab.github.io/Benin-NGS/
#This step performs a clustering analysis similar to above but 
# with a different equation.The next line creates a tree file based on 
# dissimilarity rather than relatedness.
dissMatrix  =  snpgdsIBS(genofile , sample.id=NULL, autosome.only=TRUE,
                         remove.monosnp=TRUE,  maf=NaN, missing.rate=NaN, 
                         num.thread=2, verbose=TRUE)
snpHCluster =  snpgdsHCluster(dissMatrix, sample.id=NULL, need.mat=TRUE, hang=0.01)
cutTree = snpgdsCutTree(snpHCluster, z.threshold=15, outlier.n=5, 
                        n.perm = 5000, samp.group=NULL, 
                        col.outlier="red", col.list=NULL, pch.outlier=4, 
                        pch.list=NULL,
                        label.H=FALSE, label.Z=FALSE, verbose=TRUE)
cutTree

snpgdsDrawTree(cutTree, main = "Phylogenetic Tree",edgePar=list(col=rgb(0.5,0.5,0.5,0.75),t.col="black"),
               y.label.kinship=T,leaflab="perpendicular")

# distance matrix - use IBS
dissMatrix  =  snpgdsIBS(genofile , sample.id=NULL, snp.id=snp.id, autosome.only=TRUE, 
                         remove.monosnp=TRUE,  maf=NaN, missing.rate=NaN, num.thread=2, verbose=TRUE)
snpgdsClose(genofile)

snpHCluster =  snpgdsHCluster(dissMatrix, sample.id=NULL, need.mat=TRUE, hang=0.01)

cutTree = snpgdsCutTree(snpHCluster, z.threshold=15, outlier.n=2, n.perm = 5000, samp.group=NULL, 
                        col.outlier="red", col.list=NULL, pch.outlier=2, pch.list=NULL,label.H=FALSE, label.Z=FALSE, 
                        verbose=TRUE)


library(dendextend)
library(dplyr)

par(mar = c(2,2,2,10))
cutTree %>% 
  set("labels_cex",1.2) %>% 
  plot(main="d1 (good margins)", horiz = TRUE)


pdf("Tree_plot_XL.pdf",w=8, h=8, version="1.4", bg="transparent")
    par(mar = c(15, 1, 1, 1),
        cex.axis = 1, cex.lab = 1) # Set the margin on all sides to 5
    par(oma = c(10, 1, 2, 3) + 0.1)
    #par(omi = c(5, 1, 0, 0)) # Outer margins
    snpgdsDrawTree(cutTree, main = "Dataset 1",
                   edgePar=list(col=rgb(0.5,0.5,0.5,0.75),t.col="black"),
               y.label.kinship=T,leaflab="perpendicular",
               mar = c(15, 5, 5, 5),
               cex.axis = 0.1, cex.lab = 0.1, pch=1)
dev.off()
```
