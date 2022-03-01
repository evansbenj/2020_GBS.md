# PCA

```R
library("devtools")
library(gdsfmt)
library(SNPRelate)

setwd("/Users/Shared/Previously Relocated Items/Security/projects/2020_XL_GBS/SNPRelate")

vcf.fn <- "DB_new_chr9_10L_out.vcf_filtered.vcf.gz_filtered_removed.vcf.gz"
snpgdsVCF2GDS(vcf.fn, "test.gds", method="biallelic.only")

genofile = openfn.gds("test.gds", readonly=FALSE)

samp.annot<-data.frame(pop.group = c("2014_Inhaca_10","2014_Inhaca_150","2014_Inhaca_152",
                                     "2014_Inhaca_24","2014_Inhaca_38","2014_Inhaca_52",
                                     "2014_Inhaca_65","946_Draken","993_Draken","BJE2897_Lendu_vict",
                                     "BJE3252_Cameroon_pow","BJE3508_DeDorn","BJE3509_DeDorn",
                                     "BJE3510_DeDorn","BJE3511_DeDorn","BJE3512_DeDorn","BJE3513_DeDorn",
                                     "BJE3514_DeDorn","BJE3515_DeDorn","BJE3525_Laigns","BJE3526_Laigns",
                                     "BJE3527_Laigns","BJE3528_Laigns","BJE3529_Laigns","BJE3530_Laigns",
                                     "BJE3531_Laigns","BJE3532_Laigns","BJE3533_Laigns","BJE3534_BW",
                                     "BJE3535_BW","BJE3536_BW","BJE3537_BW","BJE3538_BW","BJE3539_BW",
                                     "BJE3541_BW","BJE3542_BW","BJE3543_BW","BJE3544_BW","BJE3545_BW",
                                     "BJE3546_BW","BJE3547_GRNP","BJE3548_GRNP","BJE3549_GRNP","BJE3550_GRNP",
                                     "BJE3551_GRNP","BJE3552_GRNP","BJE3553_GRNP","BJE3554_GRNP",
                                     "BJE3573_VicW","BJE3574_VicW","BJE3575_Kimber","BJE3576_Kimber",
                                     "BJE3577_Kimber","BJE3578_Kimber","BJE3579_Kimber","BJE3580_Kimber",
                                     "BJE3581_Kimber","BJE3582_Kimber","BJE3632_Niewou","BJE3633_Niewou",
                                     "BJE3640_Niewou","BJE3641_Niewou","BJE3642_Niewou","BJE3644_Niewou",
                                     "BJE3645_Niewou","BJE3647_Niewou","BJE3654_ThreeSis","BJE3655_ThreeSis",
                                     "BJE3656_ThreeSis","BJE3657_ThreeSis","BJE3658_ThreeSis",
                                     "BJE3659_ThreeSis","BJE3660_ThreeSis","BJE3661_ThreeSis",
                                     "BJE3662_ThreeSis","BJE3663_ThreeSis","BJE3664_ThreeSis",
                                     "BJE3665_ThreeSis","BJE3666_ThreeSis","BJE3667_Citrus",
                                     "BJE3668_Citrus","BJE3669_Citrus","BJE3670_Citrus","BJE3671_Citrus",
                                     "BJE3672_Citrus","BJE3673_Citrus","BJE3674_Citrus","BJE3675_Citrus",
                                     "BJE3676_Citrus","BJE3677_Citrus","BJE3678_Citrus","JM_no_label1_Draken",
                                     "JM_no_label2_Draken","RT5_Botsw","Vred_8_Vred","amnh17260_Nigeria_pow"))

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

pdf("PCA_plot_XL.pdf",w=8, h=8, version="1.4", bg="transparent")
tab$Species <- c("laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","vict","pow","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","laevis","pow","laevis","pow")
tab$samp.color <- c("yellow","yellow","yellow","yellow","yellow","yellow","yellow","yellow","yellow","pink","black","blue","blue","blue","blue","blue","blue","blue","blue","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","green","blue","blue","blue","blue","blue","blue","blue","blue","green","green","yellow","yellow","yellow","yellow","yellow","yellow","yellow","yellow","purple","purple","purple","purple","purple","purple","purple","purple","green","green","green","green","green","green","green","green","green","green","green","green","green","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","blue","yellow","yellow","black","blue","black")
tab$samp.fieldid <- c("2014_Inhaca_10","2014_Inhaca_150","2014_Inhaca_152",
                      "2014_Inhaca_24","2014_Inhaca_38","2014_Inhaca_52",
                      "2014_Inhaca_65","946_Draken","993_Draken","BJE2897_Lendu_vict",
                      "BJE3252_Cameroon_pow","BJE3508_DeDorn","BJE3509_DeDorn",
                      "BJE3510_DeDorn","BJE3511_DeDorn","BJE3512_DeDorn","BJE3513_DeDorn",
                      "BJE3514_DeDorn","BJE3515_DeDorn","BJE3525_Laigns","BJE3526_Laigns",
                      "BJE3527_Laigns","BJE3528_Laigns","BJE3529_Laigns","BJE3530_Laigns",
                      "BJE3531_Laigns","BJE3532_Laigns","BJE3533_Laigns","BJE3534_BW",
                      "BJE3535_BW","BJE3536_BW","BJE3537_BW","BJE3538_BW","BJE3539_BW",
                      "BJE3541_BW","BJE3542_BW","BJE3543_BW","BJE3544_BW","BJE3545_BW",
                      "BJE3546_BW","BJE3547_GRNP","BJE3548_GRNP","BJE3549_GRNP","BJE3550_GRNP",
                      "BJE3551_GRNP","BJE3552_GRNP","BJE3553_GRNP","BJE3554_GRNP",
                      "BJE3573_VicW","BJE3574_VicW","BJE3575_Kimber","BJE3576_Kimber",
                      "BJE3577_Kimber","BJE3578_Kimber","BJE3579_Kimber","BJE3580_Kimber",
                      "BJE3581_Kimber","BJE3582_Kimber","BJE3632_Niewou","BJE3633_Niewou",
                      "BJE3640_Niewou","BJE3641_Niewou","BJE3642_Niewou","BJE3644_Niewou",
                      "BJE3645_Niewou","BJE3647_Niewou","BJE3654_ThreeSis","BJE3655_ThreeSis",
                      "BJE3656_ThreeSis","BJE3657_ThreeSis","BJE3658_ThreeSis",
                      "BJE3659_ThreeSis","BJE3660_ThreeSis","BJE3661_ThreeSis",
                      "BJE3662_ThreeSis","BJE3663_ThreeSis","BJE3664_ThreeSis",
                      "BJE3665_ThreeSis","BJE3666_ThreeSis","BJE3667_Citrus",
                      "BJE3668_Citrus","BJE3669_Citrus","BJE3670_Citrus","BJE3671_Citrus",
                      "BJE3672_Citrus","BJE3673_Citrus","BJE3674_Citrus","BJE3675_Citrus",
                      "BJE3676_Citrus","BJE3677_Citrus","BJE3678_Citrus","JM_no_label1_Draken",
                      "JM_no_label2_Draken","RT5_Botsw","Vred_8_Vred","amnh17260_Nigeria_pow")

d<-ggplot(data=tab, aes(x=EV1,y=-EV2, label = samp.fieldid, color = samp.color)) +
  # label axis 
  labs(x=expression("Eigenvector 1"), y=expression("-Eigenvector 2")) +
  # legend details
  scale_colour_manual(name="Species", values = c("green"="green","purple"="purple","pink"="pink", "blue"="blue","black" = "black"),breaks=c("green","purple","pink","blue","black"),labels=c("Admixted", "Niewoudtville", "X. vict","SW Cape","X. poweri"))+
  # add points and fieldID labels
  geom_text_repel(aes(EV1,-EV2, label=(samp.fieldid)), max.overlaps = 1000) + geom_point(size=4) + 
  # change to cleaner theme
  theme_classic(base_size = 16) +
  # make it clean
  theme_bw()+ theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
  # italicize species names
  theme(legend.text = element_text(face="italic"))+ 
  # move the legend
  theme(legend.position = c(.14, .28)) +
  # add a title
  ggtitle("Principal Components Analsis") + 
  # remove boxes around legend symbols
  theme(legend.key = element_blank())
d
dev.off()
```
