# JoinMap

I found the results from OneMap unsatisfying for several reasons. The main one is that I was unable to replicate findings from Bredreson etal. using the same data (they used JoinMap). As well the map distances I calculated seemed way off (too big).

So I am trying again with JoinMap. I'm using code here to convert vcf to loc format:
```
https://github.com/tomkurowski/vcf2loc
```

# Filtering
The vcf files have been filtered using these criteria:
```
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 20.0" --filter-name "QUAL20" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
```

Using vcftools I am going to filter a bit more to (1) ensure that >80% of the samples have data and (2) ensure that there is at least 8X coverage per site

```
module load StdEnv/2020 vcftools/0.1.16
vcftools --vcf XX_Chr1_removed.vcf --max-missing 0.8 --min-meanDP 8 --recode --out XX_Chr1_removed_JoinMap
```

I also need to change the length of the sample names so that they are less than 20 characters each (the '-e' is needed to make this work on OSX):
```
sed -i -e 's/__sorted.bam//g' Mitros_C659_Chr10_removed_JoinMap.recode.vcf
```

Save only GT field
```
module load bcftools/1.11   StdEnv/2020 intel/2020.1.217
bcftools annotate -x INFO,^FORMAT/GT Mitros_C659_Chr10_removed_JoinMap.recode.vcf -Ov -o Mitros_C659_Chr10_removed_JoinMap_onlyGT.vcf
```
get rid of bars
```
sed -i 's/|/\//g' Mitros_C659_Chr10_removed_JoinMap_onlyGT.vcf
```
now the python script (usually) works well:
```
/home/ben/projects/rrg-ben/ben/2022_GBS_lotsofxennies/individual_gvcfs_by_species/2023_Mitros_trop/JoinMap_C659/vcf2loc/vcf2loc.py -t CP -a SRR8704355 -b SRR8704354 -o Mitros_C659_Chr10_removed_JoinMap.recode.loc Mitros_C659_Chr10_removed_JoinMap.recode_noGQ.vcf
/home/ben/projects/rrg-ben/ben/2022_GBS_lotsofxennies/individual_gvcfs_by_species/2023_Mitros_trop/JoinMap_C659/vcf2loc/vcf2loc.py -t CP -a SRR8704346_and_7 -b SRR8704348_and_9 -o temp_new.vcf.loc temp_new_new.vcf
```
(this stalls when there are data only in the parents but not in any offspring, but this should not occur after vcftools filtering
(it also stalls when one parent is 0/1 and the other is 1/1). These sites need to be manually removed


# Potential problems (deprecated; not used)
Sometimes there are still sites in the vcf file with this field `./.:0,0:.:.`
You can see them like this
```
bcftools query -f '[\t%DP]\n' Mitros_C660_removed_Chr7_min80DP9_noGQ.vcf | grep '\.'
```
and you can locate them like this:
```
bcftools query -f '[\t%PL]\n' Mitros_C660_removed_Chr7_min80DP9_noGQ.vcf | grep '\.'
grep '0,12,141' Mitros_C660_removed_Chr7_min80DP6_noGQ.vcf | grep '0,21,230' | grep '0,36,412' | grep '0,33,405'
```
Easiest thing to do is to use sed to replace this
```
sed -i -e 's/\.:0\,0:\./0:0\,0:0/g' XL_Chr3L_JoinMap_maxmiss80DP8_noGQ.vcf
```
and maybe also:
```
sed -i -e 's/\.\/0:0\,0:0:\./0\/0:0\,0:0:0/g' XL_Chr3L_JoinMap_maxmiss80DP8_noGQ.vcf
```
and maybe this also:
```
sed -i -e 's/\.:0\,0\,0:\.:/0:0\,0\,0:0:/g' XL_Chr9_10S_JoinMap_maxmiss80DP8_noGQ.vcf
```
and maybe this also:
```
sed -i -e 's/\.:0\,0\,0\,0:\.:/0:0\,0\,0\,0:0:/g' XL_Chr9_10S_JoinMap_maxmiss80DP8_noGQ.vcf
```
and maybe this also:
```
sed -i -e 's/\.:0\,0\,0\,0\,0:\.:/0:0\,0\,0\,0\,0:0:/g' XL_Chr9_10S_JoinMap_maxmiss80DP8_noGQ.vcf
```
and this one too:
```
sed -i -e 's/\.:0\,0\,0\,0\,0\,0:\.:/0:0\,0\,0\,0\,0\,0:0:/g' XT_GW_Chr8_maxmiss80DP8_noGQ.vcf
```

and even after all this for one file (XL_Chr9_10S) I had to analyse parts of the file until I found two offending lines, which I removed, and then it worked for the whole file. This is because these two sites both had a unphased homoz alt genotype call:
```
GT:AD:DP:PGT:PID:PL:PS   1/1:0,54:54:.:.:2430,163,0:.
```

When I changed the '1\1' to '1|1' it worked.

I tried this and it fixed the problem for one input file:
```
sed -i -e 's/1\/1\:0\,/1\|1\:0\,/g'  Mitros_C659_Chr8_removed_JoinMap_min80DP8_noGQ.vcf
```
But in another there was still a problematic entry:
'1/1:1,14:15:395,15,0' which I changed to '1|1:1,14:15:395,15,0' and then it worked

# Fixed sites
Use cut to get a list of the fixed sites:
```
cut -d ' ' -f1 Mitros_C659_Chr10_removed_JoinMap.recode.loc > C659_Chr10_fixed.txt
```
This file then needs to be edited by removing some text in the beginning and end and replacing hard returns with spaces.


# JoinMap

The loc file can now be opened with JoinMap. 

ML function
* Go to the file menu and save a now project
* Load the loc file from the File menu "Load Data"
* After right clicking the yellow square node in the left most pane, in the "Dataset" Menu you need to select "Create New Dataset from Data Tabsheet"
* Now the data are loaded and you can "Check for Coding Errors" in the Dataset menu.
* Now right click on the yellow square icon and then select the "Locus Genot. Freq." tab and then click on the calculator icon in the toolbar below the part of the menu with words
* This will test for segregation distortion and the results (X2 values) are in the "Locus Genot Freq" pane in a column called "X2"; sort this by p value by clicking the "Signif" tab twice
* Highlight the significant ones by right clicking. I've been excluding the ones with 5 or more asterisks. Go to the "Population" menu and select "Exclude Marked Items". This will delete the ones with segregation distortion.
* You can confirm that they are excluded by clicking on the "Loci" tab and checking some of them - they should have the "exclude" checkbox selected
* Now click on the yellow population node again and select the "Groupings (tree)" pane in the right pane
* Click on the calculate icon in the icon toolbar; this will generate a tree with grouped markers with different stringencies (based on LOD scores but this is adjustable)
* Right click on a node in the Grouping tree (s) that you want to focus on. This should be a node with lots of markers in it and with a LOD score of at least 3.
* In the Population Menu, select "Create Groups Using the Grouping Tree". This creates another node in the leftmost pane called "Grouping 1" which has one or more groups within it
* Within a "Grouping" there are one or more "Groups" (hopefully only one). If you click on a Group and select the "Loci" tab you can exclude identical individuals. Start by excluding the second identical one and keeping the first.
* Now click on the "Fixed Orders" tab. Here you can load the orders beginning with an @ sign
* Then if you click on a group you can go to the "Group" menu option and click "Calculate Map"
* This generates a "Mapping" icon within the "Group" that has squiggly yellow pattern. Within this there is a map icon (purple) that contains the joint map plus one for each parent (1, 1_P1, 1_P2, green icons)
* Once this first round of mapping is done, you can click on the purple map icon that is the parent of the join and parental maps. Check the stress for each locus and exclude ones that are >100 by clicking on the yellow Group icon and excluding markers in the loci tab. Repeat the map with the reduced number of markers by going to the "Group" menu option and clicking "Calculate Map". After what may be several iterations of this, also check that the order of the markers in the map is ascending; exclude any that are out of order and redo the map until a final map is produced with sequential markers and low stress for all markers.
* After right clicking and selecting all of the rows, you can export this using the "Edit" menu and select "Export to File"

Regression and Kosambi function
* Go to the file menu and save a new project
* Load the loc file from the File menu "Load Data"
* After right clicking the yellow square node in the left most pane, in the "Dataset" Menu you need to select "Create New Dataset from Data Tabsheet". You probably need to try to right click it several times until it turns fuscha.
* Now the data are loaded and you can "Check for Coding Errors" in the Dataset menu.
* Click on the yellow icon with the squiggly stuff; in the "Group" tab click on "Regression Mapping"; in the "Regression Mapping" tab click on "Kosambi's" Mapping function. Then cluck "Save to Project"
* Click on the Dataset1 node or both and then in the "Dataset" menu select "Create Maternal and Paternal Population Nodes". This created two "root" GBS nodes (GBSP1, GBSP2). Each one can be processed as previosly
* Now right click on each yellow square icon (GBSP1, GBSP2) and select the "Locus Genot. Freq." tab and then click on the calculator icon in the toolbar below the part of the menu with words
* This will test for segregation distortion and the results (X2 values) are in the "Locus Genot Freq" pane in a column called "X2"; sort this by p value by clicking the "Signif" tab twice
* Highlight the significant ones by right clicking. I've been excluding the ones with 4 or more asterisks. Go to the "Population" menu and select "Exclude Marked Items". This will delete the ones with segregation distortion.
* You can confirm that they are excluded by clicking on the "Loci" tab and checking some of them - they should have the "exclude" checkbox selected
* Now click on the yellow population node again and select the "Groupings (tree)" pane in the right pane
* Click on the calculate icon in the icon toolbar; this will generate a tree with grouped markers with different stringencies (based on LOD scores but this is adjustable)
* Right click on a node in the Grouping tree (s) that you want to focus on. This should be a node with lots of markers in it and with a LOD score of at least 3.
* In the Population Menu, select "Create Groups Using the Grouping Tree". This creates another node in the leftmost pane called "Grouping 1" which has one or more groups within it
* Within a "Grouping" there are one or more "Groups" (hopefully only one). If you click on a Group and select the "Loci" tab you can exclude identical individuals, but this is not necessary for the regression mapping algorithm (it is for the ML algorithm)
* Now click on the "Fixed Orders" tab. Here you can load the orders beginning with an @ sign
* Then if you click on a group you can go to the "Group" menu option and click "Calculate Map"
* This generates a "Mapping" icon within the "Group" that has squiggly yellow pattern. Within this there is a map icon (purple) that contains three yellow maps. These are usually the same (I think) and I am using the 3rd one
* After right clicking and selecting all of the rows, you can export this using the "Edit" menu and select "Export to File"


# Plotting matpat maps
* XT
```R
library (ggplot2)
library(tidyverse)
library(stringr) # needed to split a column
library(reshape2) # this facilitates the overlay plot
library(ggformula) # for spline derivative
library(splines) # for splines
library(ggrcs) # for spline
library(rms) # needed for spline
library(scam) # needed for spline
library(mgcv)
####################
# Modified derivative.scam function
####################
derivative.scam.miso <- function (object, smooth.number = 1, deriv = 1) 
{
  sn <- smooth.number
  if (length(object$smooth[[sn]]$term) != 1) 
    stop("derivative.smooth() currently handles only 1D smooths")
  if (deriv != 1 & deriv != 2) 
    stop(paste("deriv can be either 1 or 2"))
  n <- length(object$y);n
  q <- object$smooth[[sn]]$bs.dim;q
  first <- object$smooth[[sn]]$first.para
  last <- object$smooth[[sn]]$last.para # mpi has 10; miso has only 8
  beta.t <- object$coefficients.t[first:last] # thus beta.t for mpi is length 9, beta.t for for miso is length 7
  Vp <- object$Vp.t[first:last, first:last] # dim(object$Vp.t) for mpi is 9x9; dim(object$Vp.t) for miso is 7x7
  if (inherits(object$smooth[[sn]], c("mpi.smooth", "mpd.smooth", # true for mpi
                                      "cv.smooth", "cx.smooth", "mdcv.smooth", "mdcx.smooth", 
                                      "micv.smooth", "micx.smooth"))) {
    Sig <- object$smooth[[sn]]$Sigma # this is a 9x9 matrix with "1" in the bottom diagonal for mpi
    if (deriv == 1) {
      P <- diff(diag(q - 1), difference = 1) # this is a 9x9 matrix with 1s in the diagonal
      Xd <- object$smooth[[sn]]$Xdf1 %*% P %*% Sig # this is a 3 way matrix multiplication of a vector of length 160 x two 9x9 matrixes
    }
    else {
      P <- diff(diag(q - 1), difference = 2)
      Xd <- object$smooth[[sn]]$Xdf2 %*% P %*% Sig
    }
  }
  if (inherits(object$smooth[[sn]], c("miso.smooth")))
  { # this would not be used by mpi but would be used by miso
    xx <- object$model[object$smooth[[sn]]$term]
    if (object$smooth[[sn]]$by != "NA") {
      by <- rep(1, n)
      newd <- data.frame(x = xx, by = by)
      names(newd) <- c(object$smooth[[sn]]$term, object$smooth[[sn]]$by)
    }
    else { # this is for miso
      newd <- data.frame(x = xx) # newd is 200x1
      names(newd) <- object$smooth[[sn]]$term # colname "x1"
    }
    X0 <- PredictMat(object$smooth[[sn]], newd) # X0 is 200x10 but first three columns are zeros
    eps <- 0.0000001
    xx <- xx + eps
    if (object$smooth[[sn]]$by != "NA") {
      newd <- data.frame(x = xx, by = by)
      names(newd) <- c(object$smooth[[sn]]$term, object$smooth[[sn]]$by)
    }
    else {
      newd <- data.frame(x = xx)
      names(newd) <- object$smooth[[sn]]$term
    }
    X1 <- PredictMat(object$smooth[[sn]], newd) # X1 is also 200x10 but first three columns are zeros
    Xd <- (X1 - X0)/eps # Xd is 200x10
    if (deriv == 2) {
      xx <- xx + eps
      if (object$smooth[[sn]]$by != "NA") {
        newd <- data.frame(x = xx, by = by)
        names(newd) <- c(object$smooth[[sn]]$term, object$smooth[[sn]]$by)
      }
      else {
        newd <- data.frame(x = xx)
        names(newd) <- object$smooth[[sn]]$term
      }
      X2 <- PredictMat(object$smooth[[sn]], newd) # X2 is 200x10
      Xd <- (X2 - 2 * X1 + X0)/eps^2 # Xd is 200x10
    }
  }
  df <- Xd[,-c(1:3)] %*% beta.t # here I have trimmed off the first three columns which are zeros
  df.sd <- rowSums(Xd[,-c(1:3)] %*% Vp * Xd[,-c(1:3)])^0.5 # here as well
  list(d = df, se.d = df.sd)
}

options(scipen=999)
# open data ----
setwd("/Users/Shared/Previously Relocated Items/Security/projects/2022_GBS_lotsof_Xennies/JoinMap/XT_GW/XT_GW_DP8_results")
# head -1 Chr1_jointmap.txt > all_joinmap.txt; awk 'FNR>1{print}' *jointmap*txt >> all_joinmap.txt
my_df_mat <- read.table("all_joinmap_regkos_mat.txt", header = T, sep = "\t")
my_df_pat <- read.table("all_joinmap_regkos_pat.txt", header = T, sep = "\t")
colnames(my_df_mat) <- c("SN","Nr","Locus","SegPat","Identical","Group","cM")
colnames(my_df_pat) <- c("SN","Nr","Locus","SegPat","Identical","Group","cM")
my_df_mat$matpat <- "mat"
my_df_pat$matpat <- "pat"
my_df <- rbind(my_df_mat,my_df_pat)

my_df[c('Chr', 'Coordinate')] <- str_split_fixed(my_df$Locus, '_', 2)

# add to df so that the entire Chr is plotted for each chr
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr1",217471165))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr2",181034960))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr3",153873356))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr4",153961318))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr5",164033574))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr6",154486311))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr7",133565929))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr8",147241509))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr9",91218943))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr10",52432565))

my_df$Chr <- factor(my_df$Chr,
                    levels = c("Chr1", "Chr2","Chr3",
                               "Chr4","Chr5","Chr6",
                               "Chr7","Chr8","Chr9",
                               "Chr10"), ordered = T)



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

# this uses the p/q arm ratios of Bredeson et al. BioRxiv
# Locations of trop centromeres
# Chr   Ave
# Chr1	89237096.5
# Chr2	67510626.5
# Chr3	16750087
# Chr4	46596006.5
# Chr5	62015161.5
# Chr6	73091979
# Chr7	60385499
# Chr8	21445719.5
# Chr9	42124650
# Chr10	21225599.5

# convert to numeric
my_df[, c(7,10)] <- sapply(my_df[, c(7,10)], as.numeric)

# make a Mb coordinate column
my_df$Coordinate_Mb <- my_df$Coordinate/1000000
# Calculate the derivatives for plotting first
# save the results to a big df called "master_derivative_df"

#create data frame with 0 rows and 3 columns
master_derivative_df <- data.frame(matrix(ncol = 5, nrow = 0))

#provide column names
colnames(master_derivative_df) <- c('Coordinate_Mb', 'predictions', 'df', 'matpat', 'Chr')

# Cycle through each Chr and each mat pat
for(i in unique(my_df$Chr)){
  for(j in c("mat","pat")){
    print(i)
    print(j)
    #i <- "Chr1"
    #j <- "mat"
    # select only the Chr1
    my_df_matpat_Chr_only <- my_df[(my_df$Chr==i)&(my_df$matpat==j),]
    # plot(my_df_matpat_Chr_only$Coordinate_Mb,my_df_matpat_Chr_only$cM)
    if((dim(my_df_matpat_Chr_only)[1] > 10)&
       (max(my_df_matpat_Chr_only$Coordinate_Mb)-min(my_df_matpat_Chr_only$Coordinate_Mb) >25)){
        # Build the model
        model <- scam(cM ~ s(Coordinate_Mb, 
                            bs = "miso"), 
                            data = my_df_matpat_Chr_only) 
        # add the predicted values to the df
        my_df_matpat_Chr_only$fitted.values <-model$fitted.values
        # plot(my_df_matpat_Chr_only$Coordinate_Mb,my_df_matpat_Chr_only$fitted.values)
        # Make predictions for length of chr every 5Mb
        lengths <- as.data.frame(seq(as.integer(min(my_df_matpat_Chr_only$Coordinate_Mb)),as.integer(max(my_df_matpat_Chr_only$Coordinate_Mb)),5))
        colnames(lengths) <- "Coordinate_Mb"
        predictions <- model %>% predict(lengths) # this sometimes produces some slightly negative values
        # above is the same as this:
        #  predictions <- predict(model,lengths)
        # predictions <- predict(model,lengths,type="lpmatrix")
        # now make a new df with the 5Mb coordinate increments and also the predicted values
        new_lengths <- cbind(lengths,predictions)
        # plot(new_lengths$Coordinate_Mb,new_lengths$predictions)
        # re-estimate the model using the fixed intervals and predicted values
        # model <- scam(predictions ~ s(Coordinate_Mb, bs = "mpi", k=knotz), data = new_lengths)
        model <- scam(predictions ~ s(Coordinate_Mb, 
                      bs = "miso"), 
                      data = new_lengths) # trying without the monotonic increase
        # try derivative.scam function
        d1 <- derivative.scam.miso(model,smooth.number=1,deriv=1)
        new_lengths$derivative <- d1$d
        new_lengths$matpat <- j
        new_lengths$Chr <- i
        master_derivative_df <- rbind(master_derivative_df,new_lengths)
    }    
  }
}  

#png(filename = "temp.png",w=1500, h=200,units = "px", bg="transparent")
#  ggplot(data=new_lengths, aes(x=Coordinate_Mb, y=predictions)) + 
#  geom_smooth(data = new_lengths, method = scam, 
#              formula = y ~ s(x, k = knotz), #, bs = "mpi"), # the "bs = "mpi" section seems to mess things up
#              se = FALSE) + geom_point(size=0.5) 
#dev.off()
#png(filename = "temp2.png",w=1500, h=200,units = "px", bg="transparent")
#  ggplot(data=new_lengths, aes(x=Coordinate_Mb, y=derivative)) + 
#  geom_line() + geom_point()
#dev.off()

# XT_Recombination ----
png(filename = "XT_Recombination_plot_miso.png",w=1500, h=200,units = "px", bg="transparent")
XT_Recombination<-ggplot(my_df %>% arrange(Chr), aes(x=as.numeric(Coordinate)/1000000, y=cM, colour = matpat)) + 
  # the "group = 1" part was needed to avoid a warning for some reason
  scale_color_manual(breaks = c("mat", "pat","end"), values=c("red","blue","white"), labels=c('Maternal','Paternal',' ')) +
  geom_point(size=0.5, alpha = 0.3) +
  geom_smooth(method = scam, 
              formula = y ~ s(x, 
              #k= round((max(my_df_matpat_Chr_only$Coordinate_Mb)-min(my_df_matpat_Chr_only$Coordinate_Mb))/5), 
              bs = "miso"), 
              se = FALSE) +
  scale_y_continuous(name="Map Units\n(cM)", limits=c(-20,200), breaks=c(0,50,100,150,200)) +
  scale_x_continuous(name=element_blank(), breaks=c(0,50,100,150,200,250)) + # name="Coordinates (Mb)", 
  facet_grid(. ~ Chr, scales = "free", space='free') +
  #facet_grid(~factor(Chr),scales="free_x",space = "free_x") +
  # facet_grid(~Chr,scales="free_x",space = "free_x") +
  #geom_point(data = data.frame(Coordinate = 89237096.5, cM = 0, Chr = "Chr1"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 67510626.5, cM = 0, Chr = "Chr2"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 16750087, cM = 0, Chr = "Chr3"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 46596006.5, cM = 0, Chr = "Chr4"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 62015161.5, cM = 0, Chr = "Chr5"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 73091979, cM = 0, Chr = "Chr6"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 60385499, cM = 0, Chr = "Chr7"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 21445719.5, cM = 0, Chr = "Chr8"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 42124650, cM = 0, Chr = "Chr9"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 21225599.5, cM = 0, Chr = "Chr10"), colour="black", size=3) +
  theme_classic() + theme(legend.position="none") +
  expand_limits(x = 0) +
  theme(legend.title=element_blank()) +
  theme(text = element_text(size = 18)) 
XT_Recombination
dev.off()

# Now plot derivatives

master_derivative_df$Chr <- factor(master_derivative_df$Chr,
                                   levels = c("Chr1", "Chr2","Chr3",
                                              "Chr4","Chr5","Chr6",
                                              "Chr7","Chr8","Chr9",
                                              "Chr10"), ordered = T)

# XT_Derivative ----
png(filename = "XT_Derivative_plot_miso.png",w=1500, h=200,units = "px", bg="transparent")
XT_Derivative<-ggplot(master_derivative_df %>% arrange(Chr), aes(x=Coordinate_Mb, y=derivative, col = matpat)) + 
  scale_color_manual(breaks = c("mat", "pat","end"), values=c("red","blue","white"), labels=c('Maternal','Paternal','end')) +
  geom_point(size=0.5) +
  geom_line() +
  scale_y_continuous(name="Recombination rate\n(cM/Mb)", limits=c(-1,8), breaks=seq(0,8,2)) +
  scale_x_continuous(name=element_blank(), breaks=c(0,50,100,150,200,250)) + # name="Coordinates (Mb)", 
  facet_grid(~factor(Chr),scales="free_x",space = "free_x") +
  #geom_vline(data=filter(my_df, Chr=="Chr1"), aes(xintercept=89.2370965), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr2"), aes(xintercept=67.5106265), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr3"), aes(xintercept=16.750087), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr4"), aes(xintercept=46.5960065), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr5"), aes(xintercept=62.0151615), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr6"), aes(xintercept=73.091979), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr7"), aes(xintercept=60.385499), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr8"), aes(xintercept=21.4457195), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr9"), aes(xintercept=42.124650), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr10"), aes(xintercept=21.2255995), colour="black") + 
  theme_classic() +  theme(legend.position="none") +
  expand_limits(x = 0) +
  theme(legend.title=element_blank()) +
  theme(text = element_text(size = 18)) # +
XT_Derivative 
dev.off()

# now scale all the chromosomes and then plot the slopes for 
# pat and mat on one overlay plot
# for XB use the XB genome lengths
# Chr1L	232529968
# Chr2L	184566230
# Chr3L	145564450
# Chr4L	156120766
# Chr5L	174499025
# Chr6L	157843503
# Chr7L	136892545
# Chr8L	123836260
# Chr9_10L	135078615
# Chr1S	196169797
# Chr2S	167897112
# Chr3S	127416163
# Chr4S	131359389
# Chr5S	139053355
# Chr6S	137668414
# Chr7S	105895007
# Chr8S	105436523
# Chr9_10S	110702965

XT_chrlengths <- data.frame(
  chromosome = c("Chr1", "Chr2","Chr3",
                 "Chr4","Chr5","Chr6",
                 "Chr7","Chr8","Chr9",
                 "Chr10"),
  length = c("217471165", "181034960","153873356",
             "153961318","164033574","154486311",
             "133565929","147241509","91218943",
             "52432565")
)
XT_chrlengths$length <- as.numeric(XT_chrlengths$length)

# for XL use XL
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

# for XT use XT
# Chr1	217471165
# Chr2	181034960
# Chr3	153873356
# Chr4	153961318
# Chr5	164033574
# Chr6	154486311
# Chr7	133565929
# Chr8	147241509
# Chr9	91218943
# Chr10	52432565

# get mat and pat total length
mat_cM_length <- 0
pat_cM_length <- 0
mat_Coordinate_length <- 0
pat_Coordinate_length <- 0

for(i in 1:length(XT_chrlengths$chromosome)){
  if(length(my_df[(my_df$matpat == "mat")&(my_df$Chr == XT_chrlengths$chromosome[i]),]$Coordinate) >0){
    mat_cM_length <- mat_cM_length + max(my_df[(my_df$matpat == "mat")&(my_df$Chr == XT_chrlengths$chromosome[i]),]$cM)
    pat_cM_length <- pat_cM_length + max(my_df[(my_df$matpat == "pat")&(my_df$Chr == XT_chrlengths$chromosome[i]),]$cM)
    mat_Coordinate_length <- mat_Coordinate_length + max(my_df[(my_df$matpat == "mat")&(my_df$Chr == XT_chrlengths$chromosome[i]),]$Coordinate)-
      min(my_df[(my_df$matpat == "mat")&(my_df$Chr == XT_chrlengths$chromosome[i]),]$Coordinate)
    pat_Coordinate_length <- pat_Coordinate_length + max(my_df[(my_df$matpat == "pat")&(my_df$Chr == XT_chrlengths$chromosome[i]),]$Coordinate)-
      min(my_df[(my_df$matpat == "pat")&(my_df$Chr == XT_chrlengths$chromosome[i]),]$Coordinate)
  }
}


mat_cM_length; pat_cM_length
mat_Coordinate_length/1000000; pat_Coordinate_length/1000000
mat_cM_length/(mat_Coordinate_length/1000000)
pat_cM_length/(pat_Coordinate_length/1000000)


# standardize coordinates for all chrs
# divide all by length of Chr1L
my_df$Standardized_Coordinate <- my_df$Coordinate/XT_chrlengths$length[XT_chrlengths$chromosome == "Chr1"]
# now update the other chrs
my_df$Standardized_Coordinate[my_df$Chr == "Chr2"] <-my_df$Coordinate[my_df$Chr == "Chr2"]/XT_chrlengths$length[XT_chrlengths$chromosome == "Chr2"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr3"] <-my_df$Coordinate[my_df$Chr == "Chr3"]/XT_chrlengths$length[XT_chrlengths$chromosome == "Chr3"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr4"] <-my_df$Coordinate[my_df$Chr == "Chr4"]/XT_chrlengths$length[XT_chrlengths$chromosome == "Chr4"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr5"] <-my_df$Coordinate[my_df$Chr == "Chr5"]/XT_chrlengths$length[XT_chrlengths$chromosome == "Chr5"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr6"] <-my_df$Coordinate[my_df$Chr == "Chr6"]/XT_chrlengths$length[XT_chrlengths$chromosome == "Chr6"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr7"] <-my_df$Coordinate[my_df$Chr == "Chr7"]/XT_chrlengths$length[XT_chrlengths$chromosome == "Chr7"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr8"] <-my_df$Coordinate[my_df$Chr == "Chr8"]/XT_chrlengths$length[XT_chrlengths$chromosome == "Chr8"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr9"] <-my_df$Coordinate[my_df$Chr == "Chr9"]/XT_chrlengths$length[XT_chrlengths$chromosome == "Chr9"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr10"] <-my_df$Coordinate[my_df$Chr == "Chr10"]/XT_chrlengths$length[XT_chrlengths$chromosome == "Chr10"]
# do not standardize the cM for each chromosome because different amounts of recombination occur in each one

# Now get the derivative of the generalized additive model (gam)

#create data frame with 0 rows and 5 columns
monster_derivative_df <- data.frame(matrix(ncol = 5, nrow = 0))

#provide column names
colnames(monster_derivative_df) <- c('Coordinate_Mb', 'predictions', 'df', 'matpat', 'Chr')

# Cycle through each Chr and each mat pat
for(i in unique(my_df$Chr)){
  for(j in unique(my_df$matpat)){
    print(i)
    print(j)
    #i <- "Chr2L"
    #j <- "pat"
    # select only one chromosme at a time
    # and only the mat or pat recombination events
    my_df_matpat_Chr_only <- my_df[(my_df$Chr==i)&(my_df$matpat==j),]
    if((dim(my_df_matpat_Chr_only)[1] > 10)&
       (max(my_df_matpat_Chr_only$Coordinate_Mb)-min(my_df_matpat_Chr_only$Coordinate_Mb) >25)){
      # Build the model
      # use the same number of knots as was used for the unscaled data for each chr
      # model <- gam(cM ~ s(Standardized_Coordinate, k=knotz), data = my_df_matpat_Chr_only)
      # this is the scam model with monotonic increase
      model <- scam(cM ~ s(Standardized_Coordinate,
                    bs = "miso"), 
                    data = my_df_matpat_Chr_only) 
      # add the predicted values to the df
      my_df_matpat_Chr_only$fitted.values <-model$fitted.values
      # Make predictions for length of chr every 0.005 units
      lengths <- as.data.frame(seq(0,1,0.05))
      colnames(lengths) <- "Standardized_Coordinate"
      predictions <- model %>% predict(lengths)
      # now make a new df with the 0.005 increments and also the predicted values
      # the length is 1 so this should be 200 predictions per chromosome
      new_lengths <- cbind(lengths,predictions)
      model <- scam(predictions ~ s(Standardized_Coordinate, 
                    bs = "miso"), 
                    data = new_lengths) # trying without the monotonic increase
      # get derivative using derivative.scam function
      d1 <- derivative.scam.miso(model,smooth.number=1,deriv=1)
      new_lengths$derivative <- d1$d
      new_lengths$matpat <- j
      new_lengths$Chr <- i
      # trim predictions from new_lengths that do not have data
      new_derivative_df <- subset(new_lengths, (Standardized_Coordinate > as.numeric(min(my_df_matpat_Chr_only$Coordinate))/(XT_chrlengths$length[XT_chrlengths$chromosome == i]))&
                                    (Standardized_Coordinate < as.numeric(max(my_df_matpat_Chr_only$Coordinate))/(XT_chrlengths$length[XT_chrlengths$chromosome == i])))
      monster_derivative_df <- rbind(monster_derivative_df,new_derivative_df)
    }    
  }
}  


png(filename = "XT_Combined_Scaled_Derivative_miso.png",w=1000, h=300,units = "px", bg="transparent")
trop_scale<-ggplot(monster_derivative_df, aes(x=Standardized_Coordinate, y=derivative, col = matpat)) + 
  scale_color_manual(breaks = c("mat", "pat"), values=c("red","blue"), labels=c('Maternal','Paternal')) +
  geom_point(size=0.5, alpha = 0.5) +
  geom_smooth() +
#  scale_y_continuous(name="Recombination rate\n(cM/scaled Coordinates)", limits=c(-50,800), breaks=seq(0,800,200)) +
#  scale_x_continuous(name="Scaled Coordinates", limits=c(0,1), breaks=c(0,.5,1)) +
  scale_y_continuous(name=" ", limits=c(-50,800), breaks=seq(0,800,200)) +
  scale_x_continuous(name=" ", limits=c(0,1), breaks=c(0,.5,1)) +
  #facet_grid(~factor(matpat)) +
  #geom_vline(data=filter(my_df, Chr=="Chr1"), aes(xintercept=89.2370965), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr2"), aes(xintercept=67.5106265), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr3"), aes(xintercept=16.750087), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr4"), aes(xintercept=46.5960065), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr5"), aes(xintercept=62.0151615), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr6"), aes(xintercept=73.091979), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr7"), aes(xintercept=60.385499), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr8"), aes(xintercept=21.4457195), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr9"), aes(xintercept=42.124650), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr10"), aes(xintercept=21.2255995), colour="black") + 
  theme_classic() +
  expand_limits(x = 0) +
  theme(legend.position="none") +
  theme(legend.title=element_blank()) +
  annotate("text", x=0.28, y=600, label=as.character(expression(italic("X. tropicalis"))), size=6, parse = T, hjust = 0) +
  theme(text = element_text(size = 23)) # +
trop_scale 
dev.off()


```

* XL
```R
library (ggplot2)
library(tidyverse)
library(stringr) # needed to split a column
library(reshape2) # this facilitates the overlay plot
library(ggformula) # for spline derivative
library(splines) # for splines
library(ggrcs) # for spline
library(rms) # needed for spline
library(scam) # needed for spline
library(mgcv)
####################
# Modified derivative.scam function
####################
derivative.scam.miso <- function (object, smooth.number = 1, deriv = 1) 
{
  sn <- smooth.number
  if (length(object$smooth[[sn]]$term) != 1) 
    stop("derivative.smooth() currently handles only 1D smooths")
  if (deriv != 1 & deriv != 2) 
    stop(paste("deriv can be either 1 or 2"))
  n <- length(object$y);n
  q <- object$smooth[[sn]]$bs.dim;q
  first <- object$smooth[[sn]]$first.para
  last <- object$smooth[[sn]]$last.para # mpi has 10; miso has only 8
  beta.t <- object$coefficients.t[first:last] # thus beta.t for mpi is length 9, beta.t for for miso is length 7
  Vp <- object$Vp.t[first:last, first:last] # dim(object$Vp.t) for mpi is 9x9; dim(object$Vp.t) for miso is 7x7
  if (inherits(object$smooth[[sn]], c("mpi.smooth", "mpd.smooth", # true for mpi
                                      "cv.smooth", "cx.smooth", "mdcv.smooth", "mdcx.smooth", 
                                      "micv.smooth", "micx.smooth"))) {
    Sig <- object$smooth[[sn]]$Sigma # this is a 9x9 matrix with "1" in the bottom diagonal for mpi
    if (deriv == 1) {
      P <- diff(diag(q - 1), difference = 1) # this is a 9x9 matrix with 1s in the diagonal
      Xd <- object$smooth[[sn]]$Xdf1 %*% P %*% Sig # this is a 3 way matrix multiplication of a vector of length 160 x two 9x9 matrixes
    }
    else {
      P <- diff(diag(q - 1), difference = 2)
      Xd <- object$smooth[[sn]]$Xdf2 %*% P %*% Sig
    }
  }
  if (inherits(object$smooth[[sn]], c("miso.smooth")))
  { # this would not be used by mpi but would be used by miso
    xx <- object$model[object$smooth[[sn]]$term]
    if (object$smooth[[sn]]$by != "NA") {
      by <- rep(1, n)
      newd <- data.frame(x = xx, by = by)
      names(newd) <- c(object$smooth[[sn]]$term, object$smooth[[sn]]$by)
    }
    else { # this is for miso
      newd <- data.frame(x = xx) # newd is 200x1
      names(newd) <- object$smooth[[sn]]$term # colname "x1"
    }
    X0 <- PredictMat(object$smooth[[sn]], newd) # X0 is 200x10 but first three columns are zeros
    eps <- 0.0000001
    xx <- xx + eps
    if (object$smooth[[sn]]$by != "NA") {
      newd <- data.frame(x = xx, by = by)
      names(newd) <- c(object$smooth[[sn]]$term, object$smooth[[sn]]$by)
    }
    else {
      newd <- data.frame(x = xx)
      names(newd) <- object$smooth[[sn]]$term
    }
    X1 <- PredictMat(object$smooth[[sn]], newd) # X1 is also 200x10 but first three columns are zeros
    Xd <- (X1 - X0)/eps # Xd is 200x10
    if (deriv == 2) {
      xx <- xx + eps
      if (object$smooth[[sn]]$by != "NA") {
        newd <- data.frame(x = xx, by = by)
        names(newd) <- c(object$smooth[[sn]]$term, object$smooth[[sn]]$by)
      }
      else {
        newd <- data.frame(x = xx)
        names(newd) <- object$smooth[[sn]]$term
      }
      X2 <- PredictMat(object$smooth[[sn]], newd) # X2 is 200x10
      Xd <- (X2 - 2 * X1 + X0)/eps^2 # Xd is 200x10
    }
  }
  df <- Xd[,-c(1:3)] %*% beta.t # here I have trimmed off the first three columns which are zeros
  df.sd <- rowSums(Xd[,-c(1:3)] %*% Vp * Xd[,-c(1:3)])^0.5 # here as well
  list(d = df, se.d = df.sd)
}

options(scipen=999)
# open data ----
setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/JoinMap/XL/XL_DP8_results_new")
# head -1 Chr1_jointmap.txt > all_joinmap.txt; awk 'FNR>1{print}' *jointmap*txt >> all_joinmap.txt
my_df_mat <- read.table("allXL_joinmap_regkos_mat.txt", header = T, sep = "\t")
my_df_pat <- read.table("allXL_joinmap_regkos_pat.txt", header = T, sep = "\t")
colnames(my_df_mat) <- c("SN","Nr","Locus","SegPat","Identical","Group","cM")
colnames(my_df_pat) <- c("SN","Nr","Locus","SegPat","Identical","Group","cM")
my_df_mat$matpat <- "mat"
my_df_pat$matpat <- "pat"
my_df <- rbind(my_df_mat,my_df_pat)

my_df[c('Chr', 'Coordinate')] <- str_split_fixed(my_df$Locus, '_', 2)

# This is needed to correctly parse Chr9_10L
good <- subset(my_df, (Chr == "Chr1L")|(Chr == "Chr2L")|(Chr == "Chr3L")|(Chr == "Chr4L")|
                 (Chr == "Chr5L")|(Chr == "Chr6L")|(Chr == "Chr7L")|(Chr == "Chr8L")|
                 (Chr == "Chr1S")|(Chr == "Chr2S")|(Chr == "Chr3S")|(Chr == "Chr4S")|
                 (Chr == "Chr5S")|(Chr == "Chr6S")|(Chr == "Chr7S")|(Chr == "Chr8S"));good
bad <- subset(my_df, (Chr == "Chr9"));bad
bad[c('Chr', 'Coordinate1')] <- str_split_fixed(bad$Coordinate, '_', 2);bad
bad$Chr1 <- paste("Chr9_",bad$Chr, sep="")
bad$Chr <- bad$Chr1
bad$Coordinate <- bad$Coordinate1
bad <- bad %>% select(c(-Chr1,-Coordinate1))

my_df <- rbind(good,bad)


# add to df so that the entire Chr is plotted for each chr
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr1L",233740090))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr2L",191000146))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr3L",161426101))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr4L",155250554))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr5L",171415384))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr6L",164223595))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr7L",139837618))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr8L",135449133))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr9_10L",137811819))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr1S",202412970))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr2S",169306100))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr3S",131962816))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr4S",132731174))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr5S",143394103))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr6S",137316286))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr7S",113060389))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr8S",103977862))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr9_10S",117266291))

my_df$Chr <- factor(my_df$Chr,
                    levels = c("Chr1L", "Chr2L","Chr3L",
                               "Chr4L","Chr5L","Chr6L",
                               "Chr7L","Chr8L","Chr9_10L",
                               "Chr1S", "Chr2S","Chr3S",
                               "Chr4S","Chr5S","Chr6S",
                               "Chr7S","Chr8S","Chr9_10S"), ordered = T)

my_df$LorS <- ifelse(my_df$Chr %in% c("Chr1L", "Chr2L","Chr3L",
                                       "Chr4L","Chr5L","Chr6L",
                                       "Chr7L","Chr8L","Chr9_10L"), "L", "S")



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

# this uses the p/q arm ratios of Bredeson et al. BioRxiv
# Locations of trop centromeres
# Chr   Ave
# Chr1	89237096.5
# Chr2	67510626.5
# Chr3	16750087
# Chr4	46596006.5
# Chr5	62015161.5
# Chr6	73091979
# Chr7	60385499
# Chr8	21445719.5
# Chr9	42124650
# Chr10	21225599.5

# convert to numeric
my_df[, c(7,10)] <- sapply(my_df[, c(7,10)], as.numeric)

# make a Mb coordinate column
my_df$Coordinate_Mb <- my_df$Coordinate/1000000
# Calculate the derivatives for plotting first
# save the results to a big df called "master_derivative_df"

#create data frame with 0 rows and 3 columns
master_derivative_df <- data.frame(matrix(ncol = 5, nrow = 0))

#provide column names
colnames(master_derivative_df) <- c('Coordinate_Mb', 'predictions', 'df', 'matpat', 'Chr')

# Cycle through each Chr and each mat pat
for(i in unique(my_df$Chr)){
  for(j in c("mat","pat")){
    print(i)
    print(j)
    #i <- "Chr1L"
    #j <- "mat"
    # select only the Chr1
    my_df_matpat_Chr_only <- my_df[(my_df$Chr==i)&(my_df$matpat==j),]
    # plot(my_df_matpat_Chr_only$Coordinate_Mb,my_df_matpat_Chr_only$cM)
    if((dim(my_df_matpat_Chr_only)[1] > 10)&
       (max(my_df_matpat_Chr_only$Coordinate_Mb)-min(my_df_matpat_Chr_only$Coordinate_Mb) >25)){
        # Build the model
        #knotz <- round((max(my_df_matpat_Chr_only$Coordinate_Mb)-min(my_df_matpat_Chr_only$Coordinate_Mb))/5);knotz
        #model <- gam(cM ~ s(Coordinate_Mb, k=knotz), data = my_df_matpat_Chr_only)
        model <- scam(cM ~ s(Coordinate_Mb, 
                            # k=knotz,
                            bs = "miso"), 
                            data = my_df_matpat_Chr_only) 
        # add the predicted values to the df
        my_df_matpat_Chr_only$fitted.values <-model$fitted.values
        # plot(my_df_matpat_Chr_only$Coordinate_Mb,my_df_matpat_Chr_only$fitted.values)
        # Make predictions for length of chr every 5Mb
        lengths <- as.data.frame(seq(as.integer(min(my_df_matpat_Chr_only$Coordinate_Mb)),as.integer(max(my_df_matpat_Chr_only$Coordinate_Mb)),5))
        colnames(lengths) <- "Coordinate_Mb"
        predictions <- model %>% predict(lengths) # this sometimes produces some slightly negative values
        # above is the same as this:
        #  predictions <- predict(model,lengths)
        # predictions <- predict(model,lengths,type="lpmatrix")
        # now make a new df with the 5Mb coordinate increments and also the predicted values
        new_lengths <- cbind(lengths,predictions)
        # plot(new_lengths$Coordinate_Mb,new_lengths$predictions)
        # re-estimate the model using the fixed intervals and predicted values
        # model <- scam(predictions ~ s(Coordinate_Mb, bs = "mpi", k=knotz), data = new_lengths)
        model <- scam(predictions ~ s(Coordinate_Mb, 
                      bs = "miso"), 
                      data = new_lengths) # trying without the monotonic increase
        # try derivative.scam function
        d1 <- derivative.scam.miso(model,smooth.number=1,deriv=1)
        new_lengths$derivative <- d1$d
        #x1 <- new_lengths$Coordinate_Mb
        #f1 <- new_lengths$predictions
        
        # ggplot(data=new_lengths, aes(x=Coordinate_Mb, y=predictions)) + 
        #    geom_smooth(data = new_lengths, method = scam, 
        #              formula = y ~ s(x), #, bs = "mpi"), # the "bs = "mpi" section seems to mess things up
        #              se = FALSE) + geom_point(size=0.5) 
        # ggplot(data=new_lengths, aes(x=Coordinate_Mb, y=derivative)) + 
        #  geom_line() + geom_point()
        # plot(new_lengths$Coordinate_Mb,new_lengths$predictions)
        # plot(new_lengths$Coordinate_Mb,d1$d)
        new_lengths$matpat <- j
        new_lengths$Chr <- i
        master_derivative_df <- rbind(master_derivative_df,new_lengths)
    }    
  }
}  

#png(filename = "temp.png",w=1500, h=200,units = "px", bg="transparent")
#  ggplot(data=new_lengths, aes(x=Coordinate_Mb, y=predictions)) + 
#  geom_smooth(data = new_lengths, method = scam, 
#              formula = y ~ s(x, k = knotz), #, bs = "mpi"), # the "bs = "mpi" section seems to mess things up
#              se = FALSE) + geom_point(size=0.5) 
#dev.off()
#png(filename = "temp2.png",w=1500, h=200,units = "px", bg="transparent")
#  ggplot(data=new_lengths, aes(x=Coordinate_Mb, y=derivative)) + 
#  geom_line() + geom_point()
#dev.off()

# XL_Recombination_L ----
# subset L
my_df_Lonly <- my_df[my_df$LorS == "L",]
png(filename = "XL_Recombination_plot_miso_L.png",w=1500, h=200,units = "px", bg="transparent")
XL_Recombination_L<-ggplot(my_df_Lonly %>% arrange(Chr), aes(x=as.numeric(Coordinate)/1000000, y=cM, colour = matpat)) + 
  # the "group = 1" part was needed to avoid a warning for some reason
  scale_color_manual(breaks = c("mat", "pat","end"), values=c("red","blue","white"), labels=c('Maternal','Paternal',' ')) +
  geom_point(size=0.5, alpha = 0.3) +
  geom_smooth(method = scam, 
              formula = y ~ s(x, 
              #k= round((max(my_df_matpat_Chr_only$Coordinate_Mb)-min(my_df_matpat_Chr_only$Coordinate_Mb))/5), 
              bs = "miso"), 
              se = FALSE) +
  scale_y_continuous(name="Map Units\n(cM)", limits=c(-20,200), breaks=c(0,50,100,150,200)) +
  scale_x_continuous(name=element_blank(), breaks=c(0,50,100,150,200,250)) +
  facet_grid(. ~ Chr, scales = "free", space='free') +
  #facet_grid(~factor(Chr),scales="free_x",space = "free_x") +
  # facet_grid(~Chr,scales="free_x",space = "free_x") +
  #geom_point(data = data.frame(Coordinate = 89237096.5, cM = 0, Chr = "Chr1"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 67510626.5, cM = 0, Chr = "Chr2"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 16750087, cM = 0, Chr = "Chr3"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 46596006.5, cM = 0, Chr = "Chr4"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 62015161.5, cM = 0, Chr = "Chr5"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 73091979, cM = 0, Chr = "Chr6"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 60385499, cM = 0, Chr = "Chr7"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 21445719.5, cM = 0, Chr = "Chr8"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 42124650, cM = 0, Chr = "Chr9"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 21225599.5, cM = 0, Chr = "Chr10"), colour="black", size=3) +
  theme_classic() + theme(legend.position="none") +
  expand_limits(x = 0) +
  theme(legend.title=element_blank()) +
  theme(text = element_text(size = 18)) 
XL_Recombination_L 
dev.off()

# XL_Recombination_S ----
# subset S
my_df_Sonly <- my_df[my_df$LorS == "S",]
png(filename = "XL_Recombination_plot_miso_S.png",w=1500, h=200,units = "px", bg="transparent")
XL_Recombination_S<-ggplot(my_df_Sonly %>% arrange(Chr), aes(x=as.numeric(Coordinate)/1000000, y=cM, colour = matpat)) + 
  # the "group = 1" part was needed to avoid a warning for some reason
  scale_color_manual(breaks = c("mat", "pat","end"), values=c("red","blue","white"), labels=c('Maternal','Paternal','end')) +
  geom_point(size=0.5, alpha = 0.3) +
  geom_smooth(method = scam, 
              formula = y ~ s(x, 
             # k= round((max(my_df_matpat_Chr_only$Coordinate_Mb)-min(my_df_matpat_Chr_only$Coordinate_Mb))/5), 
              bs = "miso"), 
              se = FALSE) +
  scale_y_continuous(name="Map Units\n(cM)", limits=c(-20,200), breaks=c(0,50,100,150,200)) +
  scale_x_continuous(name=element_blank(), breaks=c(0,50,100,150,200,250)) +
  facet_grid(. ~ Chr, scales = "free", space='free') +
  #facet_grid(~factor(Chr),scales="free_x",space = "free_x") +
  # facet_grid(~Chr,scales="free_x",space = "free_x") +
  #geom_point(data = data.frame(Coordinate = 89237096.5, cM = 0, Chr = "Chr1"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 67510626.5, cM = 0, Chr = "Chr2"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 16750087, cM = 0, Chr = "Chr3"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 46596006.5, cM = 0, Chr = "Chr4"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 62015161.5, cM = 0, Chr = "Chr5"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 73091979, cM = 0, Chr = "Chr6"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 60385499, cM = 0, Chr = "Chr7"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 21445719.5, cM = 0, Chr = "Chr8"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 42124650, cM = 0, Chr = "Chr9"), colour="black", size=3) +
#geom_point(data = data.frame(Coordinate = 21225599.5, cM = 0, Chr = "Chr10"), colour="black", size=3) +
theme_classic() +
  expand_limits(x = 0) +
  theme(legend.title=element_blank()) +  theme(legend.position="none") +
  theme(text = element_text(size = 18)) # +
# guides(color = guide_legend(override.aes = list(size = 5)))
# Get rid of the legend
#theme(legend.position = "none")
XL_Recombination_S 
dev.off()


# Now plot derivatives

master_derivative_df$Chr <- factor(master_derivative_df$Chr,
                    levels = c("Chr1L", "Chr2L","Chr3L",
                               "Chr4L","Chr5L","Chr6L",
                               "Chr7L","Chr8L","Chr9_10L",
                               "Chr1S", "Chr2S","Chr3S",
                               "Chr4S","Chr5S","Chr6S",
                               "Chr7S","Chr8S","Chr9_10S"), ordered = T)

master_derivative_df$LorS <- ifelse(master_derivative_df$Chr %in% 
                                      c("Chr1L", "Chr2L","Chr3L",
                                      "Chr4L","Chr5L","Chr6L",
                                      "Chr7L","Chr8L","Chr9_10L"), "L", "S")

# XL_Derivative_L ----
# subset L
master_derivative_df_Lonly <- master_derivative_df[master_derivative_df$LorS == "L",]

png(filename = "XL_Derivative_plot_miso_Lonly.png",w=1500, h=200,units = "px", bg="transparent")
XL_Derivative_L<-ggplot(master_derivative_df_Lonly %>% arrange(Chr), aes(x=Coordinate_Mb, y=derivative, col = matpat)) + 
  scale_color_manual(breaks = c("mat", "pat","end"), values=c("red","blue","white"), labels=c('Maternal','Paternal','end')) +
  geom_point(size=0.5) +
  geom_line() +
  scale_y_continuous(name="Recombination rate\n(cM/Mb)", limits=c(-1,8), breaks=seq(0,8,2)) +
  scale_x_continuous(name=element_blank(), breaks=c(0,50,100,150,200,250)) +
  facet_grid(~factor(Chr),scales="free_x",space = "free_x") +
  #geom_vline(data=filter(my_df, Chr=="Chr1"), aes(xintercept=89.2370965), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr2"), aes(xintercept=67.5106265), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr3"), aes(xintercept=16.750087), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr4"), aes(xintercept=46.5960065), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr5"), aes(xintercept=62.0151615), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr6"), aes(xintercept=73.091979), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr7"), aes(xintercept=60.385499), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr8"), aes(xintercept=21.4457195), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr9"), aes(xintercept=42.124650), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr10"), aes(xintercept=21.2255995), colour="black") + 
  theme_classic() +  theme(legend.position="none") +
  expand_limits(x = 0) +
  theme(legend.title=element_blank()) +
  theme(text = element_text(size = 18)) # +
XL_Derivative_L
dev.off()

# XL_Derivative_S ----
# subset S
master_derivative_df_Sonly <- master_derivative_df[master_derivative_df$LorS == "S",]

png(filename = "XL_Derivative_plot_miso_Sonly.png",w=1500, h=200,units = "px", bg="transparent")
XL_Derivative_S<-ggplot(master_derivative_df_Sonly %>% arrange(Chr), aes(x=Coordinate_Mb, y=derivative, col = matpat)) + 
  scale_color_manual(breaks = c("mat", "pat","end"), values=c("red","blue","white"), labels=c('Maternal','Paternal','end')) +
  geom_point(size=0.5) +
  geom_line() +
  scale_y_continuous(name="Recombination rate\n(cM/Mb)", limits=c(-1,8), breaks=seq(0,8,2)) +
  scale_x_continuous(name=element_blank(), breaks=c(0,50,100,150,200,250)) +
  facet_grid(~factor(Chr),scales="free_x",space = "free_x") +
  #geom_vline(data=filter(my_df, Chr=="Chr1"), aes(xintercept=89.2370965), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr2"), aes(xintercept=67.5106265), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr3"), aes(xintercept=16.750087), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr4"), aes(xintercept=46.5960065), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr5"), aes(xintercept=62.0151615), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr6"), aes(xintercept=73.091979), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr7"), aes(xintercept=60.385499), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr8"), aes(xintercept=21.4457195), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr9"), aes(xintercept=42.124650), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr10"), aes(xintercept=21.2255995), colour="black") + 
  theme_classic() +  theme(legend.position="none") +
  expand_limits(x = 0) +
  theme(legend.title=element_blank()) +
  theme(text = element_text(size = 18)) # +
XL_Derivative_S
dev.off()

# now scale all the chromosomes and then plot the slopes for 
# pat and mat on one overlay plot
# for XB use the XB genome lengths
# Chr1L	232529968
# Chr2L	184566230
# Chr3L	145564450
# Chr4L	156120766
# Chr5L	174499025
# Chr6L	157843503
# Chr7L	136892545
# Chr8L	123836260
# Chr9_10L	135078615
# Chr1S	196169797
# Chr2S	167897112
# Chr3S	127416163
# Chr4S	131359389
# Chr5S	139053355
# Chr6S	137668414
# Chr7S	105895007
# Chr8S	105436523
# Chr9_10S	110702965

XL_chrlengths <- data.frame(
  chromosome = c("Chr1L", "Chr2L","Chr3L",
    "Chr4L","Chr5L","Chr6L",
    "Chr7L","Chr8L","Chr9_10L",
    "Chr1S", "Chr2S","Chr3S",
    "Chr4S","Chr5S","Chr6S",
    "Chr7S","Chr8S","Chr9_10S"),
  length = c("233740090", "191000146","161426101",
    "155250554","171415384","164223595",
    "139837618","135449133","137811819",
    "202412970", "169306100","131962816",
    "132731174","143394103","137316286",
    "113060389","103977862","117266291")
  )
XL_chrlengths$length <- as.numeric(XL_chrlengths$length)

# for XL use XL
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

# for XT use XT
# Chr1	217471165
# Chr2	181034960
# Chr3	153873356
# Chr4	153961318
# Chr5	164033574
# Chr6	154486311
# Chr7	133565929
# Chr8	147241509
# Chr9	91218943
# Chr10	52432565

# get mat and pat total length
mat_cM_length <- 0
pat_cM_length <- 0
mat_Coordinate_length <- 0
pat_Coordinate_length <- 0

mat_cM_length_L <- 0
pat_cM_length_L <- 0
mat_Coordinate_length_L <- 0
pat_Coordinate_length_L <- 0

mat_cM_length_S <- 0
pat_cM_length_S <- 0
mat_Coordinate_length_S <- 0
pat_Coordinate_length_S <- 0

for(i in 1:length(XL_chrlengths$chromosome)){
  if(length(my_df[(my_df$matpat == "mat")&(my_df$Chr == XL_chrlengths$chromosome[i]),]$Coordinate) >0){
    mat_cM_length <- mat_cM_length + max(my_df[(my_df$matpat == "mat")&(my_df$Chr == XL_chrlengths$chromosome[i]),]$cM)
    pat_cM_length <- pat_cM_length + max(my_df[(my_df$matpat == "pat")&(my_df$Chr == XL_chrlengths$chromosome[i]),]$cM)
    mat_Coordinate_length <- mat_Coordinate_length + max(my_df[(my_df$matpat == "mat")&(my_df$Chr == XL_chrlengths$chromosome[i]),]$Coordinate)-
      min(my_df[(my_df$matpat == "mat")&(my_df$Chr == XL_chrlengths$chromosome[i]),]$Coordinate)
    pat_Coordinate_length <- pat_Coordinate_length + max(my_df[(my_df$matpat == "pat")&(my_df$Chr == XL_chrlengths$chromosome[i]),]$Coordinate)-
      min(my_df[(my_df$matpat == "pat")&(my_df$Chr == XL_chrlengths$chromosome[i]),]$Coordinate)
    if(i %in% 1:9){
        mat_cM_length_L <- mat_cM_length_L + max(my_df[(my_df$matpat == "mat")&(my_df$Chr == XL_chrlengths$chromosome[i]),]$cM)
        pat_cM_length_L <- pat_cM_length_L + max(my_df[(my_df$matpat == "pat")&(my_df$Chr == XL_chrlengths$chromosome[i]),]$cM)
        mat_Coordinate_length_L <- mat_Coordinate_length_L + max(my_df[(my_df$matpat == "mat")&(my_df$Chr == XL_chrlengths$chromosome[i]),]$Coordinate)-
          min(my_df[(my_df$matpat == "mat")&(my_df$Chr == XL_chrlengths$chromosome[i]),]$Coordinate)
        pat_Coordinate_length_L <- pat_Coordinate_length_L + max(my_df[(my_df$matpat == "pat")&(my_df$Chr == XL_chrlengths$chromosome[i]),]$Coordinate)-
          min(my_df[(my_df$matpat == "pat")&(my_df$Chr == XL_chrlengths$chromosome[i]),]$Coordinate)
    }
    else if(i %in% 10:18){
      mat_cM_length_S <- mat_cM_length_S + max(my_df[(my_df$matpat == "mat")&(my_df$Chr == XL_chrlengths$chromosome[i]),]$cM)
      pat_cM_length_S <- pat_cM_length_S + max(my_df[(my_df$matpat == "pat")&(my_df$Chr == XL_chrlengths$chromosome[i]),]$cM)
      mat_Coordinate_length_S <- mat_Coordinate_length_S + max(my_df[(my_df$matpat == "mat")&(my_df$Chr == XL_chrlengths$chromosome[i]),]$Coordinate)-
        min(my_df[(my_df$matpat == "mat")&(my_df$Chr == XL_chrlengths$chromosome[i]),]$Coordinate)
      pat_Coordinate_length_S <- pat_Coordinate_length_S + max(my_df[(my_df$matpat == "pat")&(my_df$Chr == XL_chrlengths$chromosome[i]),]$Coordinate)-
        min(my_df[(my_df$matpat == "pat")&(my_df$Chr == XL_chrlengths$chromosome[i]),]$Coordinate)
    }
  }
}


mat_cM_length; pat_cM_length
mat_Coordinate_length/1000000; pat_Coordinate_length/1000000
mat_cM_length/(mat_Coordinate_length/1000000)
pat_cM_length/(pat_Coordinate_length/1000000)

mat_cM_length_L; pat_cM_length_L
mat_Coordinate_length_L/1000000; pat_Coordinate_length_L/1000000
mat_cM_length_L/(mat_Coordinate_length_L/1000000)
pat_cM_length_L/(pat_Coordinate_length_L/1000000)

mat_cM_length_S; pat_cM_length_S
mat_Coordinate_length_S/1000000; pat_Coordinate_length_S/1000000
mat_cM_length_S/(mat_Coordinate_length_S/1000000)
pat_cM_length_S/(pat_Coordinate_length_S/1000000)


# standardize coordinates for all chrs
# divide all by length of Chr1L
my_df$Standardized_Coordinate <- my_df$Coordinate/XL_chrlengths$length[XL_chrlengths$chromosome == "Chr1L"]
# now update the other chrs
my_df$Standardized_Coordinate[my_df$Chr == "Chr2L"] <-my_df$Coordinate[my_df$Chr == "Chr2L"]/XL_chrlengths$length[XL_chrlengths$chromosome == "Chr2L"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr3L"] <-my_df$Coordinate[my_df$Chr == "Chr3L"]/XL_chrlengths$length[XL_chrlengths$chromosome == "Chr3L"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr4L"] <-my_df$Coordinate[my_df$Chr == "Chr4L"]/XL_chrlengths$length[XL_chrlengths$chromosome == "Chr4L"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr5L"] <-my_df$Coordinate[my_df$Chr == "Chr5L"]/XL_chrlengths$length[XL_chrlengths$chromosome == "Chr5L"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr6L"] <-my_df$Coordinate[my_df$Chr == "Chr6L"]/XL_chrlengths$length[XL_chrlengths$chromosome == "Chr6L"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr7L"] <-my_df$Coordinate[my_df$Chr == "Chr7L"]/XL_chrlengths$length[XL_chrlengths$chromosome == "Chr7L"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr8L"] <-my_df$Coordinate[my_df$Chr == "Chr8L"]/XL_chrlengths$length[XL_chrlengths$chromosome == "Chr8L"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr9_10L"] <-my_df$Coordinate[my_df$Chr == "Chr9_10L"]/XL_chrlengths$length[XL_chrlengths$chromosome == "Chr9_10L"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr1S"] <-my_df$Coordinate[my_df$Chr == "Chr1S"]/XL_chrlengths$length[XL_chrlengths$chromosome == "Chr1S"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr2S"] <-my_df$Coordinate[my_df$Chr == "Chr2S"]/XL_chrlengths$length[XL_chrlengths$chromosome == "Chr2S"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr3S"] <-my_df$Coordinate[my_df$Chr == "Chr3S"]/XL_chrlengths$length[XL_chrlengths$chromosome == "Chr3S"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr4S"] <-my_df$Coordinate[my_df$Chr == "Chr4S"]/XL_chrlengths$length[XL_chrlengths$chromosome == "Chr4S"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr5S"] <-my_df$Coordinate[my_df$Chr == "Chr5S"]/XL_chrlengths$length[XL_chrlengths$chromosome == "Chr5S"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr6S"] <-my_df$Coordinate[my_df$Chr == "Chr6S"]/XL_chrlengths$length[XL_chrlengths$chromosome == "Chr6S"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr7S"] <-my_df$Coordinate[my_df$Chr == "Chr7S"]/XL_chrlengths$length[XL_chrlengths$chromosome == "Chr7S"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr8S"] <-my_df$Coordinate[my_df$Chr == "Chr8S"]/XL_chrlengths$length[XL_chrlengths$chromosome == "Chr8S"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr9_10S"] <-my_df$Coordinate[my_df$Chr == "Chr9_10S"]/XL_chrlengths$length[XL_chrlengths$chromosome == "Chr9_10S"]

# do not standardize the cM for each chromosome because different amounts of recombination occur in each one

# Now get the derivative of the generalized additive model (gam)

#create data frame with 0 rows and 5 columns
monster_derivative_df <- data.frame(matrix(ncol = 5, nrow = 0))

#provide column names
colnames(monster_derivative_df) <- c('Coordinate_Mb', 'predictions', 'df', 'matpat', 'Chr')

# Cycle through each Chr and each mat pat
for(i in unique(my_df$Chr)){
  for(j in unique(my_df$matpat)){
    print(i)
    print(j)
    #i <- "Chr2L"
    #j <- "pat"
    # select only one chromosme at a time
    # and only the mat or pat recombination events
    my_df_matpat_Chr_only <- my_df[(my_df$Chr==i)&(my_df$matpat==j),]
    if((dim(my_df_matpat_Chr_only)[1] > 10)&
       (max(my_df_matpat_Chr_only$Coordinate_Mb)-min(my_df_matpat_Chr_only$Coordinate_Mb) >25)){
      # Build the model
      # use the same number of knots as was used for the unscaled data for each chr
      # model <- gam(cM ~ s(Standardized_Coordinate, k=knotz), data = my_df_matpat_Chr_only)
      # this is the scam model with monotonic increase
      model <- scam(cM ~ s(Standardized_Coordinate,
                    bs = "miso"), 
                    data = my_df_matpat_Chr_only) 
      # add the predicted values to the df
      my_df_matpat_Chr_only$fitted.values <-model$fitted.values
      # Make predictions for length of chr every 0.005 units
      lengths <- as.data.frame(seq(0,1,0.05))
      colnames(lengths) <- "Standardized_Coordinate"
      predictions <- model %>% predict(lengths)
      # now make a new df with the 0.005 increments and also the predicted values
      # the length is 1 so this should be 200 predictions per chromosome
      new_lengths <- cbind(lengths,predictions)
      model <- scam(predictions ~ s(Standardized_Coordinate, 
                    bs = "miso"), 
                    data = new_lengths) # trying without the monotonic increase
      # get derivative using derivative.scam function
      d1 <- derivative.scam.miso(model,smooth.number=1,deriv=1)
      new_lengths$derivative <- d1$d
      new_lengths$matpat <- j
      new_lengths$Chr <- i
      # trim predictions from new_lengths that do not have data
      new_derivative_df <- subset(new_lengths, (Standardized_Coordinate > as.numeric(min(my_df_matpat_Chr_only$Coordinate))/(XL_chrlengths$length[XL_chrlengths$chromosome == i]))&
                                    (Standardized_Coordinate < as.numeric(max(my_df_matpat_Chr_only$Coordinate))/(XL_chrlengths$length[XL_chrlengths$chromosome == i])))
      monster_derivative_df <- rbind(monster_derivative_df,new_derivative_df)
    }    
  }
}  


png(filename = "XL_Combined_Scaled_Derivative_miso.png",w=300, h=300,units = "px", bg="transparent")
XL_scale<-ggplot(monster_derivative_df, aes(x=Standardized_Coordinate, y=derivative, col = matpat)) + 
  scale_color_manual(breaks = c("mat", "pat"), values=c("red","blue"), labels=c('Maternal','Paternal')) +
  geom_point(size=0.5, alpha = 0.5) +
  geom_smooth() +
#  scale_y_continuous(name="Recombination rate\n(cM/scaled Coordinates)", limits=c(-50,800), breaks=seq(0,800,200)) +
#  scale_x_continuous(name="Scaled Coordinates", limits=c(0,1), breaks=c(0,.5,1)) +
  scale_y_continuous(name=" ", limits=c(-50,800), breaks=seq(0,800,200)) +
  scale_x_continuous(name=" ", limits=c(0,1), breaks=c(0,.5,1)) +
  #facet_grid(~factor(matpat)) +
  #geom_vline(data=filter(my_df, Chr=="Chr1"), aes(xintercept=89.2370965), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr2"), aes(xintercept=67.5106265), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr3"), aes(xintercept=16.750087), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr4"), aes(xintercept=46.5960065), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr5"), aes(xintercept=62.0151615), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr6"), aes(xintercept=73.091979), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr7"), aes(xintercept=60.385499), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr8"), aes(xintercept=21.4457195), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr9"), aes(xintercept=42.124650), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr10"), aes(xintercept=21.2255995), colour="black") + 
  theme_classic() +
  expand_limits(x = 0) +
  theme(legend.position="none") +
  theme(legend.title=element_blank()) +
  theme(axis.text.y = element_blank()) +
  annotate("text", x=0.32, y=600, label=as.character(expression(italic("X. laevis"))), size=6, parse = T, hjust = 0) +
  theme(text = element_text(size = 23)) # +
XL_scale
dev.off()


# Fig 4 combining plots ----
library(gridExtra)
title1=text_grob("Recombination rate\n(cM/scaled Coordinates)", size = 24, rot = 90)  
title2=text_grob("Scaled Coordinates", size =24)  
#title3=text_grob(expression(italic("X. tropicalis")), size =20)  
#title4=text_grob(expression(italic("X. laevis")), size =20)  
#title5=text_grob(expression(italic("X. borealis")), size =20)  
jpeg("./Fig4_Scaled_recombination.jpg",w=12, h=4, units ="in", bg="transparent", res = 200)
  grid.arrange(arrangeGrob(trop_scale),#, top=title3),
             arrangeGrob(XL_scale),#, top=title4), 
             arrangeGrob(XB_scale),#, top=title5),
             ncol = 3, left=title1, bottom=title2, widths=c(1.15,1,1))
dev.off()

# SI Fig 8 combining plots ----
library(gridExtra)
title1=text_grob(expression(italic("X. tropicalis")), size =20, x = 0.09, hjust = 0)  
title2=text_grob(expression(paste(italic("X. laevis")," L")), size =20, x = 0.09, hjust = 0)  
title3=text_grob(expression(paste(italic("X. laevis")," S")), size =20, x = 0.09, hjust = 0)  
title4=text_grob(expression(paste(italic("X. borealis")," L")), size =20, x = 0.09, hjust = 0)  
title5=text_grob(expression(paste(italic("X. borealis")," S")), size =20, x = 0.09, hjust = 0)  
jpeg("./SI_Fig8_Scaled_recombination.jpg",w=14, h=12, units ="in", bg="transparent", res = 200)
  grid.arrange(arrangeGrob(XT_Recombination, top=title1),
             arrangeGrob(XL_Recombination_L, top=title2), 
             arrangeGrob(XL_Recombination_S, top=title3),
             arrangeGrob(XB_Recombination_L, top=title4), 
             arrangeGrob(XB_Recombination_S, top=title5),
             ncol = 1, nrow = 5)#, left=title1, bottom=title2)
dev.off()

# SI Fig 9 combining plots ----
library(gridExtra)
title1=text_grob(expression(italic("X. tropicalis")), size =20, x = 0.09, hjust = 0)  
title2=text_grob(expression(paste(italic("X. laevis")," L")), size =20, x = 0.09, hjust = 0)  
title3=text_grob(expression(paste(italic("X. laevis")," S")), size =20, x = 0.09, hjust = 0)  
title4=text_grob(expression(paste(italic("X. borealis")," L")), size =20, x = 0.09, hjust = 0)  
title5=text_grob(expression(paste(italic("X. borealis")," S")), size =20, x = 0.09, hjust = 0)  
jpeg("./SI_Fig9_Derivative.jpg",w=14, h=12, units ="in", bg="transparent", res = 200)
grid.arrange(arrangeGrob(XT_Derivative, top=title1),
             arrangeGrob(XL_Derivative_L, top=title2), 
             arrangeGrob(XL_Derivative_S, top=title3),
             arrangeGrob(XB_Derivative_L, top=title4), 
             arrangeGrob(XB_Derivative_S, top=title5),
             ncol = 1, nrow = 5)#, left=title1, bottom=title2)
dev.off()

```
* XB
```R
library (ggplot2)
library(tidyverse)
library(stringr) # needed to split a column
library(reshape2) # this facilitates the overlay plot
library(ggformula) # for spline derivative
library(splines) # for splines
library(ggrcs) # for spline
library(rms) # needed for spline
library(scam) # needed for spline
library(mgcv)
####################
# Modified derivative.scam function
####################
derivative.scam.miso <- function (object, smooth.number = 1, deriv = 1) 
{
  sn <- smooth.number
  if (length(object$smooth[[sn]]$term) != 1) 
    stop("derivative.smooth() currently handles only 1D smooths")
  if (deriv != 1 & deriv != 2) 
    stop(paste("deriv can be either 1 or 2"))
  n <- length(object$y);n
  q <- object$smooth[[sn]]$bs.dim;q
  first <- object$smooth[[sn]]$first.para
  last <- object$smooth[[sn]]$last.para # mpi has 10; miso has only 8
  beta.t <- object$coefficients.t[first:last] # thus beta.t for mpi is length 9, beta.t for for miso is length 7
  Vp <- object$Vp.t[first:last, first:last] # dim(object$Vp.t) for mpi is 9x9; dim(object$Vp.t) for miso is 7x7
  if (inherits(object$smooth[[sn]], c("mpi.smooth", "mpd.smooth", # true for mpi
                                      "cv.smooth", "cx.smooth", "mdcv.smooth", "mdcx.smooth", 
                                      "micv.smooth", "micx.smooth"))) {
    Sig <- object$smooth[[sn]]$Sigma # this is a 9x9 matrix with "1" in the bottom diagonal for mpi
    if (deriv == 1) {
      P <- diff(diag(q - 1), difference = 1) # this is a 9x9 matrix with 1s in the diagonal
      Xd <- object$smooth[[sn]]$Xdf1 %*% P %*% Sig # this is a 3 way matrix multiplication of a vector of length 160 x two 9x9 matrixes
    }
    else {
      P <- diff(diag(q - 1), difference = 2)
      Xd <- object$smooth[[sn]]$Xdf2 %*% P %*% Sig
    }
  }
  if (inherits(object$smooth[[sn]], c("miso.smooth")))
  { # this would not be used by mpi but would be used by miso
    xx <- object$model[object$smooth[[sn]]$term]
    if (object$smooth[[sn]]$by != "NA") {
      by <- rep(1, n)
      newd <- data.frame(x = xx, by = by)
      names(newd) <- c(object$smooth[[sn]]$term, object$smooth[[sn]]$by)
    }
    else { # this is for miso
      newd <- data.frame(x = xx) # newd is 200x1
      names(newd) <- object$smooth[[sn]]$term # colname "x1"
    }
    X0 <- PredictMat(object$smooth[[sn]], newd) # X0 is 200x10 but first three columns are zeros
    eps <- 0.0000001
    xx <- xx + eps
    if (object$smooth[[sn]]$by != "NA") {
      newd <- data.frame(x = xx, by = by)
      names(newd) <- c(object$smooth[[sn]]$term, object$smooth[[sn]]$by)
    }
    else {
      newd <- data.frame(x = xx)
      names(newd) <- object$smooth[[sn]]$term
    }
    X1 <- PredictMat(object$smooth[[sn]], newd) # X1 is also 200x10 but first three columns are zeros
    Xd <- (X1 - X0)/eps # Xd is 200x10
    if (deriv == 2) {
      xx <- xx + eps
      if (object$smooth[[sn]]$by != "NA") {
        newd <- data.frame(x = xx, by = by)
        names(newd) <- c(object$smooth[[sn]]$term, object$smooth[[sn]]$by)
      }
      else {
        newd <- data.frame(x = xx)
        names(newd) <- object$smooth[[sn]]$term
      }
      X2 <- PredictMat(object$smooth[[sn]], newd) # X2 is 200x10
      Xd <- (X2 - 2 * X1 + X0)/eps^2 # Xd is 200x10
    }
  }
  df <- Xd[,-c(1:3)] %*% beta.t # here I have trimmed off the first three columns which are zeros
  df.sd <- rowSums(Xd[,-c(1:3)] %*% Vp * Xd[,-c(1:3)])^0.5 # here as well
  list(d = df, se.d = df.sd)
}
# open data ----
options(scipen=999)
setwd("/Users/Shared/Previously Relocated Items/Security/projects/2022_GBS_lotsof_Xennies/JoinMap/XB/XB_DP8_results_new")
# head -1 Chr1_jointmap.txt > all_joinmap.txt; awk 'FNR>1{print}' *jointmap*txt >> all_joinmap.txt
my_df_mat <- read.table("XB_all_joinmap_regkosDP8_mat.txt", header = T, sep = "\t")
my_df_pat <- read.table("XB_all_joinmap_regkosDP8_pat.txt", header = T, sep = "\t")
colnames(my_df_mat) <- c("SN","Nr","Locus","SegPat","Identical","Group","cM")
colnames(my_df_pat) <- c("SN","Nr","Locus","SegPat","Identical","Group","cM")
my_df_mat$matpat <- "mat"
my_df_pat$matpat <- "pat"
my_df <- rbind(my_df_mat,my_df_pat)

my_df[c('Chr', 'Coordinate')] <- str_split_fixed(my_df$Locus, '_', 2)

# This is needed to correctly parse Chr9_10L
good <- subset(my_df, (Chr == "Chr1L")|(Chr == "Chr2L")|(Chr == "Chr3L")|(Chr == "Chr4L")|
                 (Chr == "Chr5L")|(Chr == "Chr6L")|(Chr == "Chr7L")|(Chr == "Chr8L")|
                 (Chr == "Chr1S")|(Chr == "Chr2S")|(Chr == "Chr3S")|(Chr == "Chr4S")|
                 (Chr == "Chr5S")|(Chr == "Chr6S")|(Chr == "Chr7S")|(Chr == "Chr8S"));good
bad <- subset(my_df, (Chr == "Chr9"));bad
bad[c('Chr', 'Coordinate1')] <- str_split_fixed(bad$Coordinate, '_', 2);bad
bad$Chr1 <- paste("Chr9_",bad$Chr, sep="")
bad$Chr <- bad$Chr1
bad$Coordinate <- bad$Coordinate1
bad <- bad %>% select(c(-Chr1,-Coordinate1))

my_df <- rbind(good,bad)

# add to df so that the entire Chr is plotted for each chr
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr1L",232529968))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr2L",184566230))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr3L",145564450))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr4L",156120766))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr5L",174499025))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr6L",157843503))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr7L",136892545))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr8L",123836260))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr9_10L",135078615))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr1S",196169797))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr2S",167897112))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr3S",127416163))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr4S",131359389))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr5S",139053355))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr6S",137668414))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr7S",105895007))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr8S",105436523))
my_df <- rbind(my_df,c("NA","NA","NA","NA","NA","NA","NA","end","Chr9_10S",110702965))



my_df$Chr <- factor(my_df$Chr,
                    levels = c("Chr1L", "Chr2L","Chr3L",
                               "Chr4L","Chr5L","Chr6L",
                               "Chr7L","Chr8L","Chr9_10L",
                               "Chr1S", "Chr2S","Chr3S",
                               "Chr4S","Chr5S","Chr6S",
                               "Chr7S","Chr8S","Chr9_10S"), ordered = T)

my_df$LorS <- ifelse(my_df$Chr %in% c("Chr1L", "Chr2L","Chr3L",
                                       "Chr4L","Chr5L","Chr6L",
                                       "Chr7L","Chr8L","Chr9_10L"), "L", "S")



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

# this uses the p/q arm ratios of Bredeson et al. BioRxiv
# Locations of trop centromeres
# Chr   Ave
# Chr1	89237096.5
# Chr2	67510626.5
# Chr3	16750087
# Chr4	46596006.5
# Chr5	62015161.5
# Chr6	73091979
# Chr7	60385499
# Chr8	21445719.5
# Chr9	42124650
# Chr10	21225599.5

# convert to numeric
my_df[, c(7,10)] <- sapply(my_df[, c(7,10)], as.numeric)

# make a Mb coordinate column
my_df$Coordinate_Mb <- my_df$Coordinate/1000000
# Calculate the derivatives for plotting first
# save the results to a big df called "master_derivative_df"


#create data frame with 0 rows and 3 columns
master_derivative_df <- data.frame(matrix(ncol = 5, nrow = 0))

#provide column names
colnames(master_derivative_df) <- c('Coordinate_Mb', 'predictions', 'df', 'matpat', 'Chr')


# Cycle through each Chr and each mat pat
for(i in unique(my_df$Chr)){
  for(j in c("mat","pat")){
    print(i)
    print(j)
    #i <- "Chr1L"
    #j <- "mat"
    # select only the Chr1
    my_df_matpat_Chr_only <- my_df[(my_df$Chr==i)&(my_df$matpat==j),]
    # plot(my_df_matpat_Chr_only$Coordinate_Mb,my_df_matpat_Chr_only$cM)
    if((dim(my_df_matpat_Chr_only)[1] > 10)&
       (max(my_df_matpat_Chr_only$Coordinate_Mb)-min(my_df_matpat_Chr_only$Coordinate_Mb) >30)){
        # Build the model
        knotz <- round((max(my_df_matpat_Chr_only$Coordinate_Mb)-min(my_df_matpat_Chr_only$Coordinate_Mb))/6);knotz
        #model <- gam(cM ~ s(Coordinate_Mb, k=knotz), data = my_df_matpat_Chr_only)
        model <- scam(cM ~ s(Coordinate_Mb, bs = "miso"), 
                     # k = knotz,
                      data = my_df_matpat_Chr_only) 
        # add the predicted values to the df
        my_df_matpat_Chr_only$fitted.values <-model$fitted.values
        # plot(my_df_matpat_Chr_only$Coordinate_Mb,my_df_matpat_Chr_only$fitted.values)
        # Make predictions for length of chr every 5Mb
        lengths <- as.data.frame(seq(as.integer(min(my_df_matpat_Chr_only$Coordinate_Mb)),as.integer(max(my_df_matpat_Chr_only$Coordinate_Mb)),5))
        colnames(lengths) <- "Coordinate_Mb"
        predictions <- model %>% predict(lengths) # this sometimes produces some slightly negative values
        # now make a new df with the 5Mb coordinate increments and also the predicted values
        new_lengths <- cbind(lengths,predictions)
        # plot(new_lengths$Coordinate_Mb,new_lengths$predictions)
        # re-estimate the model using the fixed intervals and predicted values
        # model <- scam(predictions ~ s(Coordinate_Mb, bs = "mpi", k=knotz), data = new_lengths)
        model <- scam(predictions ~ s(Coordinate_Mb, 
                     # k=round((max(my_df_matpat_Chr_only$Coordinate_Mb)-min(my_df_matpat_Chr_only$Coordinate_Mb))/6), 
                      bs = "miso"), data = new_lengths) # trying without the monotonic increase
        # try derivative.scam function
        d1 <- derivative.scam.miso(model,smooth.number=1,deriv=1)
        new_lengths$derivative <- d1$d
        #x1 <- new_lengths$Coordinate_Mb
        #f1 <- new_lengths$predictions
        
        # ggplot(data=new_lengths, aes(x=Coordinate_Mb, y=predictions)) + 
        #    geom_smooth(data = new_lengths, method = scam, 
        #              formula = y ~ s(x), #, bs = "mpi"), # the "bs = "mpi" section seems to mess things up
        #              se = FALSE) + geom_point(size=0.5) 
        # ggplot(data=new_lengths, aes(x=Coordinate_Mb, y=derivative)) + 
        #  geom_line() + geom_point()
        # plot(new_lengths$Coordinate_Mb,new_lengths$predictions)
        # plot(new_lengths$Coordinate_Mb,d1$d)
        new_lengths$matpat <- j
        new_lengths$Chr <- i
        master_derivative_df <- rbind(master_derivative_df,new_lengths)
    }    
  }
}  

#png(filename = "temp.png",w=1500, h=200,units = "px", bg="transparent")
#  ggplot(data=new_lengths, aes(x=Coordinate_Mb, y=predictions)) + 
#  geom_smooth(data = new_lengths, method = scam, 
#              formula = y ~ s(x, k = knotz), #, bs = "mpi"), # the "bs = "mpi" section seems to mess things up
#              se = FALSE) + geom_point(size=0.5) 
#dev.off()
#png(filename = "temp2.png",w=1500, h=200,units = "px", bg="transparent")
#  ggplot(data=new_lengths, aes(x=Coordinate_Mb, y=derivative)) + 
#  geom_line() + geom_point()
#dev.off()

# subset L
my_df_Lonly <- my_df[my_df$LorS == "L",]
png(filename = "XB_Recombination_plot_miso_L.png",w=1500, h=200,units = "px", bg="transparent")
XB_Recombination_L<-ggplot(my_df_Lonly %>% arrange(Chr), aes(x=as.numeric(Coordinate)/1000000, y=cM, colour = matpat)) + 
  # the "group = 1" part was needed to avoid a warning for some reason
  scale_color_manual(breaks = c("mat", "pat","end"), values=c("red","blue","white"), labels=c('Maternal','Paternal',' ')) +
  geom_point(size=0.5, alpha = 0.3) +
  geom_smooth(method = scam, 
              formula = y ~ s(x, 
             #  k= round((max(my_df_matpat_Chr_only$Coordinate_Mb)-min(my_df_matpat_Chr_only$Coordinate_Mb))/6), 
              bs = "miso"), 
              se = FALSE,
              data = subset(my_df_Lonly, matpat == "pat")) +
  geom_smooth(method = scam, 
              formula = y ~ s(x, 
             # k= round((max(my_df_matpat_Chr_only$Coordinate_Mb)-min(my_df_matpat_Chr_only$Coordinate_Mb))/6), 
              bs = "miso"), 
              se = FALSE,
              data = subset(my_df_Lonly, ((matpat == "mat")&(Chr != "Chr8L")))) +
  scale_y_continuous(name="Map Units\n(cM)", limits=c(-20,200), breaks=c(0,50,100,150,200)) +
  scale_x_continuous(name=element_blank(), breaks=c(0,50,100,150,200,250)) +
  #geom_hline(yintercept=0) +
  #geom_hline(yintercept=c(-0.5,0.5), linetype='dashed', color=c('black', 'black'))+
  # get rid of gray background
  #facet_wrap(~Chr, nrow=2, scales = "free_x")+
  facet_grid(. ~ Chr, scales = "free", space='free') +
  #facet_grid(~factor(Chr),scales="free_x",space = "free_x") +
  # facet_grid(~Chr,scales="free_x",space = "free_x") +
  #geom_point(data = data.frame(Coordinate = 89237096.5, cM = 0, Chr = "Chr1"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 67510626.5, cM = 0, Chr = "Chr2"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 16750087, cM = 0, Chr = "Chr3"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 46596006.5, cM = 0, Chr = "Chr4"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 62015161.5, cM = 0, Chr = "Chr5"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 73091979, cM = 0, Chr = "Chr6"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 60385499, cM = 0, Chr = "Chr7"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 21445719.5, cM = 0, Chr = "Chr8"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 42124650, cM = 0, Chr = "Chr9"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 21225599.5, cM = 0, Chr = "Chr10"), colour="black", size=3) +
  theme_classic() + theme(legend.position="none") +
  expand_limits(x = 0) +
  theme(legend.title=element_blank()) +
  theme(text = element_text(size = 18)) # +
  # guides(color = guide_legend(override.aes = list(size = 5)))
# Get rid of the legend
#theme(legend.position = "none")
XB_Recombination_L 
dev.off()

# subset S
my_df_Sonly <- my_df[my_df$LorS == "S",]
png(filename = "XB_Recombination_plot_miso_S.png",w=1500, h=200,units = "px", bg="transparent")
XB_Recombination_S<-ggplot(my_df_Sonly %>% arrange(Chr), aes(x=as.numeric(Coordinate)/1000000, y=cM, colour = matpat)) + 
  # the "group = 1" part was needed to avoid a warning for some reason
  scale_color_manual(breaks = c("mat", "pat","end"), values=c("red","blue","white"), labels=c('Maternal','Paternal','end')) +
  geom_point(size=0.5, alpha = 0.3) +
  geom_smooth(method = scam, 
              formula = y ~ s(x, 
              #k= round((max(my_df_matpat_Chr_only$Coordinate_Mb)-min(my_df_matpat_Chr_only$Coordinate_Mb))/5), 
              bs = "miso"), 
              se = FALSE,
              data = subset(my_df_Sonly, matpat == "pat")) +
  geom_smooth(method = scam, 
              formula = y ~ s(x, 
              #k= round((max(my_df_matpat_Chr_only$Coordinate_Mb)-min(my_df_matpat_Chr_only$Coordinate_Mb))/5), 
              bs = "miso"), 
              se = FALSE,
              data = subset(my_df_Sonly, ((matpat == "mat")&(Chr != "Chr5S")))) +
  scale_y_continuous(name="Map Units\n(cM)", limits=c(-20,200), breaks=c(0,50,100,150,200)) +
  scale_x_continuous(name="Coordinates (Mb)", breaks=c(0,50,100,150,200,250)) +
  #geom_hline(yintercept=0) +
  #geom_hline(yintercept=c(-0.5,0.5), linetype='dashed', color=c('black', 'black'))+
  # get rid of gray background
  #facet_wrap(~Chr, nrow=2, scales = "free_x")+
  facet_grid(. ~ Chr, scales = "free", space='free') +
  #facet_grid(~factor(Chr),scales="free_x",space = "free_x") +
  # facet_grid(~Chr,scales="free_x",space = "free_x") +
  #geom_point(data = data.frame(Coordinate = 89237096.5, cM = 0, Chr = "Chr1"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 67510626.5, cM = 0, Chr = "Chr2"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 16750087, cM = 0, Chr = "Chr3"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 46596006.5, cM = 0, Chr = "Chr4"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 62015161.5, cM = 0, Chr = "Chr5"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 73091979, cM = 0, Chr = "Chr6"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 60385499, cM = 0, Chr = "Chr7"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 21445719.5, cM = 0, Chr = "Chr8"), colour="black", size=3) +
  #geom_point(data = data.frame(Coordinate = 42124650, cM = 0, Chr = "Chr9"), colour="black", size=3) +
#geom_point(data = data.frame(Coordinate = 21225599.5, cM = 0, Chr = "Chr10"), colour="black", size=3) +
theme_classic() +
  expand_limits(x = 0) +
  theme(legend.title=element_blank()) +  theme(legend.position="none") +
  theme(text = element_text(size = 18)) # +
# guides(color = guide_legend(override.aes = list(size = 5)))
# Get rid of the legend
#theme(legend.position = "none")
XB_Recombination_S 
dev.off()

# try a subset with no facet
# i <- "Chr1L"
# my_df_matpat_Chr_only <- my_df[(my_df$Chr==i),]
# p<-ggplot(my_df_matpat_Chr_only, aes(x=as.numeric(Coordinate)/1000000, y=cM, colour = matpat)) + 
#   geom_point(size=0.5, alpha = 0.2) +
  #geom_smooth(method = gam, 
  #            formula = y ~ x, 
  #            se = FALSE) + 
#   geom_smooth(method = scam, 
#               data = my_df_matpat_Chr_only,
#               formula = y ~ s(x), 
#               se = FALSE) +
  #scale_y_continuous(name="Map Units (cM)", limits=c(0,200), breaks=c(0,50,100,150,200)) +
  #scale_x_continuous(name="Coordinates (Mb)", breaks=c(0,50,100,150,200,250)) +
  #facet_grid(. ~ Chr, scales = "free", space='free') +
#   theme_classic() 
  #expand_limits(x = 0) +
  #theme(legend.title=element_blank()) +
  #theme(text = element_text(size = 12)) # +
# p 

# Now plot derivatives

master_derivative_df$Chr <- factor(master_derivative_df$Chr,
                    levels = c("Chr1L", "Chr2L","Chr3L",
                               "Chr4L","Chr5L","Chr6L",
                               "Chr7L","Chr8L","Chr9_10L",
                               "Chr1S", "Chr2S","Chr3S",
                               "Chr4S","Chr5S","Chr6S",
                               "Chr7S","Chr8S","Chr9_10S"), ordered = T)

master_derivative_df$LorS <- ifelse(master_derivative_df$Chr %in% 
                                      c("Chr1L", "Chr2L","Chr3L",
                                      "Chr4L","Chr5L","Chr6L",
                                      "Chr7L","Chr8L","Chr9_10L"), "L", "S")

# subset L
master_derivative_df_Lonly <- master_derivative_df[master_derivative_df$LorS == "L",]

png(filename = "XB_Derivative_plot_miso_Lonly.png",w=1500, h=200,units = "px", bg="transparent")
XB_Derivative_L<-ggplot(master_derivative_df_Lonly %>% arrange(Chr), aes(x=Coordinate_Mb, y=derivative, col = matpat)) + 
  scale_color_manual(breaks = c("mat", "pat","end"), values=c("red","blue","white"), labels=c('Maternal','Paternal','end')) +
  geom_point(size=0.5) +
  geom_line() +
  scale_y_continuous(name="Recombination rate\n(cM/Mb)", limits=c(-1,8), breaks=seq(0,8,2)) +
  scale_x_continuous(name=element_blank(), breaks=c(0,50,100,150,200,250)) +
  facet_grid(~factor(Chr),scales="free_x",space = "free_x") +
  #geom_vline(data=filter(my_df, Chr=="Chr1"), aes(xintercept=89.2370965), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr2"), aes(xintercept=67.5106265), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr3"), aes(xintercept=16.750087), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr4"), aes(xintercept=46.5960065), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr5"), aes(xintercept=62.0151615), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr6"), aes(xintercept=73.091979), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr7"), aes(xintercept=60.385499), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr8"), aes(xintercept=21.4457195), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr9"), aes(xintercept=42.124650), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr10"), aes(xintercept=21.2255995), colour="black") + 
  theme_classic() +  theme(legend.position="none") +
  expand_limits(x = 0) +
  theme(legend.title=element_blank()) +
  theme(text = element_text(size = 18)) # +
XB_Derivative_L
dev.off()

# subset S
master_derivative_df_Sonly <- master_derivative_df[master_derivative_df$LorS == "S",]

png(filename = "XB_Derivative_plot_miso_Sonly.png",w=1500, h=200,units = "px", bg="transparent")
XB_Derivative_S<-ggplot(master_derivative_df_Sonly %>% arrange(Chr), aes(x=Coordinate_Mb, y=derivative, col = matpat)) + 
  scale_color_manual(breaks = c("mat", "pat","end"), values=c("red","blue","white"), labels=c('Maternal','Paternal','end')) +
  geom_point(size=0.5) +
  geom_line() +
  scale_y_continuous(name="Recombination rate\n(cM/Mb)", limits=c(-1,8), breaks=seq(0,8,2)) +
  scale_x_continuous(name="Coordinates (Mb)", breaks=c(0,50,100,150,200,250)) +
  facet_grid(~factor(Chr),scales="free_x",space = "free_x") +
  #geom_vline(data=filter(my_df, Chr=="Chr1"), aes(xintercept=89.2370965), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr2"), aes(xintercept=67.5106265), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr3"), aes(xintercept=16.750087), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr4"), aes(xintercept=46.5960065), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr5"), aes(xintercept=62.0151615), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr6"), aes(xintercept=73.091979), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr7"), aes(xintercept=60.385499), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr8"), aes(xintercept=21.4457195), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr9"), aes(xintercept=42.124650), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr10"), aes(xintercept=21.2255995), colour="black") + 
  theme_classic() +  theme(legend.position="none") +
  expand_limits(x = 0) +
  theme(legend.title=element_blank()) +
  theme(text = element_text(size = 18)) # +
XB_Derivative_S
dev.off()

# now scale all the chromosomes and then plot the slopes for 
# pat and mat on one overlay plot
# for XB use the XB genome lengths
# Chr1L	232529968
# Chr2L	184566230
# Chr3L	145564450
# Chr4L	156120766
# Chr5L	174499025
# Chr6L	157843503
# Chr7L	136892545
# Chr8L	123836260
# Chr9_10L	135078615
# Chr1S	196169797
# Chr2S	167897112
# Chr3S	127416163
# Chr4S	131359389
# Chr5S	139053355
# Chr6S	137668414
# Chr7S	105895007
# Chr8S	105436523
# Chr9_10S	110702965

XB_chrlengths <- data.frame(
  chromosome = c("Chr1L", "Chr2L","Chr3L",
                 "Chr4L","Chr5L","Chr6L",
                 "Chr7L","Chr8L","Chr9_10L",
                 "Chr1S", "Chr2S","Chr3S",
                 "Chr4S","Chr5S","Chr6S",
                 "Chr7S","Chr8S","Chr9_10S"),
  length = c("232529968", "184566230","145564450",
             "156120766","174499025","157843503",
             "136892545","123836260","135078615",
             "196169797", "167897112","127416163",
             "131359389","139053355","137668414",
             "105895007","105436523","110702965")
)
XB_chrlengths$length <- as.numeric(XB_chrlengths$length)

# for XL use XL
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

# for XT use XT
# Chr1	217471165
# Chr2	181034960
# Chr3	153873356
# Chr4	153961318
# Chr5	164033574
# Chr6	154486311
# Chr7	133565929
# Chr8	147241509
# Chr9	91218943
# Chr10	52432565

# get mat and pat total length
mat_cM_length <- 0
pat_cM_length <- 0
mat_Coordinate_length <- 0
pat_Coordinate_length <- 0

mat_cM_length_L <- 0
pat_cM_length_L <- 0
mat_Coordinate_length_L <- 0
pat_Coordinate_length_L <- 0

mat_cM_length_S <- 0
pat_cM_length_S <- 0
mat_Coordinate_length_S <- 0
pat_Coordinate_length_S <- 0

for(i in 1:length(XB_chrlengths$chromosome)){
  if(length(my_df[(my_df$matpat == "mat")&(my_df$Chr == XB_chrlengths$chromosome[i]),]$Coordinate) >0){
    mat_cM_length <- mat_cM_length + max(my_df[(my_df$matpat == "mat")&(my_df$Chr == XB_chrlengths$chromosome[i]),]$cM)
    pat_cM_length <- pat_cM_length + max(my_df[(my_df$matpat == "pat")&(my_df$Chr == XB_chrlengths$chromosome[i]),]$cM)
    mat_Coordinate_length <- mat_Coordinate_length + max(my_df[(my_df$matpat == "mat")&(my_df$Chr == XB_chrlengths$chromosome[i]),]$Coordinate)-
      min(my_df[(my_df$matpat == "mat")&(my_df$Chr == XB_chrlengths$chromosome[i]),]$Coordinate)
    pat_Coordinate_length <- pat_Coordinate_length + max(my_df[(my_df$matpat == "pat")&(my_df$Chr == XB_chrlengths$chromosome[i]),]$Coordinate)-
      min(my_df[(my_df$matpat == "pat")&(my_df$Chr == XB_chrlengths$chromosome[i]),]$Coordinate)
    if(i %in% 1:9){
        mat_cM_length_L <- mat_cM_length_L + max(my_df[(my_df$matpat == "mat")&(my_df$Chr == XB_chrlengths$chromosome[i]),]$cM)
        pat_cM_length_L <- pat_cM_length_L + max(my_df[(my_df$matpat == "pat")&(my_df$Chr == XB_chrlengths$chromosome[i]),]$cM)
        mat_Coordinate_length_L <- mat_Coordinate_length_L + max(my_df[(my_df$matpat == "mat")&(my_df$Chr == XB_chrlengths$chromosome[i]),]$Coordinate)-
          min(my_df[(my_df$matpat == "mat")&(my_df$Chr == XB_chrlengths$chromosome[i]),]$Coordinate)
        pat_Coordinate_length_L <- pat_Coordinate_length_L + max(my_df[(my_df$matpat == "pat")&(my_df$Chr == XB_chrlengths$chromosome[i]),]$Coordinate)-
          min(my_df[(my_df$matpat == "pat")&(my_df$Chr == XB_chrlengths$chromosome[i]),]$Coordinate)
    }
    else if(i %in% 10:18){
      mat_cM_length_S <- mat_cM_length_S + max(my_df[(my_df$matpat == "mat")&(my_df$Chr == XB_chrlengths$chromosome[i]),]$cM)
      pat_cM_length_S <- pat_cM_length_S + max(my_df[(my_df$matpat == "pat")&(my_df$Chr == XB_chrlengths$chromosome[i]),]$cM)
      mat_Coordinate_length_S <- mat_Coordinate_length_S + max(my_df[(my_df$matpat == "mat")&(my_df$Chr == XB_chrlengths$chromosome[i]),]$Coordinate)-
        min(my_df[(my_df$matpat == "mat")&(my_df$Chr == XB_chrlengths$chromosome[i]),]$Coordinate)
      pat_Coordinate_length_S <- pat_Coordinate_length_S + max(my_df[(my_df$matpat == "pat")&(my_df$Chr == XB_chrlengths$chromosome[i]),]$Coordinate)-
        min(my_df[(my_df$matpat == "pat")&(my_df$Chr == XB_chrlengths$chromosome[i]),]$Coordinate)
    }
  }
}


mat_cM_length; pat_cM_length
mat_Coordinate_length/1000000; pat_Coordinate_length/1000000
mat_cM_length/(mat_Coordinate_length/1000000)
pat_cM_length/(pat_Coordinate_length/1000000)

mat_cM_length_L; pat_cM_length_L
mat_Coordinate_length_L/1000000; pat_Coordinate_length_L/1000000
mat_cM_length_L/(mat_Coordinate_length_L/1000000)
pat_cM_length_L/(pat_Coordinate_length_L/1000000)

mat_cM_length_S; pat_cM_length_S
mat_Coordinate_length_S/1000000; pat_Coordinate_length_S/1000000
mat_cM_length_S/(mat_Coordinate_length_S/1000000)
pat_cM_length_S/(pat_Coordinate_length_S/1000000)


# standardize coordinates for all chrs
# divide all by length of Chr1L
my_df$Standardized_Coordinate <- my_df$Coordinate/XB_chrlengths$length[XB_chrlengths$chromosome == "Chr1L"]
# now update the other chrs
my_df$Standardized_Coordinate[my_df$Chr == "Chr2L"] <-my_df$Coordinate[my_df$Chr == "Chr2L"]/XB_chrlengths$length[XB_chrlengths$chromosome == "Chr2L"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr3L"] <-my_df$Coordinate[my_df$Chr == "Chr3L"]/XB_chrlengths$length[XB_chrlengths$chromosome == "Chr3L"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr4L"] <-my_df$Coordinate[my_df$Chr == "Chr4L"]/XB_chrlengths$length[XB_chrlengths$chromosome == "Chr4L"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr5L"] <-my_df$Coordinate[my_df$Chr == "Chr5L"]/XB_chrlengths$length[XB_chrlengths$chromosome == "Chr5L"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr6L"] <-my_df$Coordinate[my_df$Chr == "Chr6L"]/XB_chrlengths$length[XB_chrlengths$chromosome == "Chr6L"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr7L"] <-my_df$Coordinate[my_df$Chr == "Chr7L"]/XB_chrlengths$length[XB_chrlengths$chromosome == "Chr7L"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr8L"] <-my_df$Coordinate[my_df$Chr == "Chr8L"]/XB_chrlengths$length[XB_chrlengths$chromosome == "Chr8L"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr9_10L"] <-my_df$Coordinate[my_df$Chr == "Chr9_10L"]/XB_chrlengths$length[XB_chrlengths$chromosome == "Chr9_10L"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr1S"] <-my_df$Coordinate[my_df$Chr == "Chr1S"]/XB_chrlengths$length[XB_chrlengths$chromosome == "Chr1S"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr2S"] <-my_df$Coordinate[my_df$Chr == "Chr2S"]/XB_chrlengths$length[XB_chrlengths$chromosome == "Chr2S"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr3S"] <-my_df$Coordinate[my_df$Chr == "Chr3S"]/XB_chrlengths$length[XB_chrlengths$chromosome == "Chr3S"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr4S"] <-my_df$Coordinate[my_df$Chr == "Chr4S"]/XB_chrlengths$length[XB_chrlengths$chromosome == "Chr4S"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr5S"] <-my_df$Coordinate[my_df$Chr == "Chr5S"]/XB_chrlengths$length[XB_chrlengths$chromosome == "Chr5S"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr6S"] <-my_df$Coordinate[my_df$Chr == "Chr6S"]/XB_chrlengths$length[XB_chrlengths$chromosome == "Chr6S"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr7S"] <-my_df$Coordinate[my_df$Chr == "Chr7S"]/XB_chrlengths$length[XB_chrlengths$chromosome == "Chr7S"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr8S"] <-my_df$Coordinate[my_df$Chr == "Chr8S"]/XB_chrlengths$length[XB_chrlengths$chromosome == "Chr8S"]
my_df$Standardized_Coordinate[my_df$Chr == "Chr9_10S"] <-my_df$Coordinate[my_df$Chr == "Chr9_10S"]/XB_chrlengths$length[XB_chrlengths$chromosome == "Chr9_10S"]

# do not standardize the cM for each chromosome because different amounts of recombination occur in each one

# Now get the derivative of the generalized additive model (gam)

#create data frame with 0 rows and 5 columns
monster_derivative_df <- data.frame(matrix(ncol = 5, nrow = 0))

#provide column names
colnames(monster_derivative_df) <- c('Coordinate_Mb', 'predictions', 'df', 'matpat', 'Chr')

# Cycle through each Chr and each mat pat
for(i in unique(my_df$Chr)){
  for(j in unique(my_df$matpat)){
    print(i)
    print(j)
    #i <- "Chr2L"
    #j <- "pat"
    # select only one chromosme at a time
    # and only the mat or pat recombination events
    my_df_matpat_Chr_only <- my_df[(my_df$Chr==i)&(my_df$matpat==j),]
    if((dim(my_df_matpat_Chr_only)[1] > 10)&
       (max(my_df_matpat_Chr_only$Coordinate_Mb)-min(my_df_matpat_Chr_only$Coordinate_Mb) >25)){
      # Build the model
      # use the same number of knots as was used for the unscaled data for each chr
      # this is the scam model with monotonic increase
      model <- scam(cM ~ s(Standardized_Coordinate,
                          # k= round((max(my_df_matpat_Chr_only$Coordinate_Mb)-min(my_df_matpat_Chr_only$Coordinate_Mb))/6),
                           bs = "miso"), 
                          data = my_df_matpat_Chr_only) 
      # add the predicted values to the df
      my_df_matpat_Chr_only$fitted.values <-model$fitted.values
      # Make predictions for length of chr every 0.005 units
      lengths <- as.data.frame(seq(0,1,0.05))
      colnames(lengths) <- "Standardized_Coordinate"
      predictions <- model %>% predict(lengths)
      # now make a new df with the 0.005 increments and also the predicted values
      # the length is 1 so this should be 200 predictions per chromosome
      new_lengths <- cbind(lengths,predictions)
      model <- scam(predictions ~ s(Standardized_Coordinate, 
                   # k= round((max(my_df_matpat_Chr_only$Coordinate_Mb)-min(my_df_matpat_Chr_only$Coordinate_Mb))/6),
                    bs = "miso"), 
                    data = new_lengths) # trying without the monotonic increase
      # get derivative using derivative.scam function
      d1 <- derivative.scam.miso(model,smooth.number=1,deriv=1)
      new_lengths$derivative <- d1$d
      new_lengths$matpat <- j
      new_lengths$Chr <- i
      # trim predictions from new_lengths that do not have data
      new_derivative_df <- subset(new_lengths, (Standardized_Coordinate > as.numeric(min(my_df_matpat_Chr_only$Coordinate))/(XB_chrlengths$length[XB_chrlengths$chromosome == i]))&
                                    (Standardized_Coordinate < as.numeric(max(my_df_matpat_Chr_only$Coordinate))/(XB_chrlengths$length[XB_chrlengths$chromosome == i])))
      monster_derivative_df <- rbind(monster_derivative_df,new_derivative_df)
    }    
  }
}  


png(filename = "XB_Combined_Scaled_Derivative_miso.png",w=300, h=300,units = "px", bg="transparent")
XB_scale<-ggplot(monster_derivative_df, aes(x=Standardized_Coordinate, y=derivative, col = matpat)) + 
  scale_color_manual(breaks = c("mat", "pat"), values=c("red","blue"), labels=c('Maternal','Paternal')) +
  geom_point(size=0.5, alpha = 0.5) +
  geom_smooth() +
#  scale_y_continuous(name="Recombination rate\n(cM/scaled Coordinates)", limits=c(-50,800), breaks=seq(0,800,200)) +
#  scale_x_continuous(name="Scaled Coordinates", limits=c(0,1), breaks=c(0,.5,1)) +
  scale_y_continuous(name=" ", limits=c(-50,800), breaks=seq(0,800,200)) +
  scale_x_continuous(name=" ", limits=c(0,1), breaks=c(0,.5,1)) +
  #facet_grid(~factor(matpat)) +
  #geom_vline(data=filter(my_df, Chr=="Chr1"), aes(xintercept=89.2370965), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr2"), aes(xintercept=67.5106265), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr3"), aes(xintercept=16.750087), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr4"), aes(xintercept=46.5960065), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr5"), aes(xintercept=62.0151615), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr6"), aes(xintercept=73.091979), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr7"), aes(xintercept=60.385499), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr8"), aes(xintercept=21.4457195), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr9"), aes(xintercept=42.124650), colour="black") + 
  #geom_vline(data=filter(my_df, Chr=="Chr10"), aes(xintercept=21.2255995), colour="black") + 
  theme_classic() +
  expand_limits(x = 0) +
  theme(legend.position="none") +
  theme(legend.title=element_blank()) +
  theme(axis.text.y = element_blank()) +
  annotate("text", x=0.32, y=600, label=as.character(expression(italic("X. borealis"))), size=6, parse = T, hjust = 0) +
  theme(text = element_text(size = 23)) # +
XB_scale 
dev.off()



```
