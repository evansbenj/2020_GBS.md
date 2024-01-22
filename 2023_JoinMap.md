# JoinMap

I found the results from OneMap unsatisfying for several reasons. The main one is that I was unable to replicate findings from Bredreson etal. using the same data (they used JoinMap). As well the map distances I calculated seemed way off (too big).

So I am trying again with JoinMap.

# Filtering
The vcf files have been filtered using these criteria:
```
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 20.0" --filter-name "QUAL20" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
```

Using vcftools I am going to filter a bit more to (1) ensure that >80% of the samples have data and (2) ensure that there is at least 10X coverage per site

```
module load StdEnv/2020 vcftools/0.1.16
vcftools --vcf XX_Chr1_removed.vcf --max-missing 0.8 --min-meanDP 8 --recode --out XX_Chr1_removed_JoinMap
```

I also need to change the length of the sample names so that they are less than 20 characters each (the '-e' is needed to make this work on OSX):
```
sed -i -e 's/__sorted.bam//g' Mitros_C659_Chr10_removed_JoinMap.recode.vcf
```

# Remove GQ field
Now I want to use a python script (https://github.com/tomkurowski/vcf2loc) to convert the vcf file to a JoinMap input file (loc). But there was an error because the GQ field has non-numeric values ('.') sometimes. I'm using bcftools to remove this field in hopes this will fix the problem:
```
module load bcftools/1.11   StdEnv/2020 intel/2020.1.217
bcftools annotate -x FORMAT/GQ Mitros_C659_Chr10_removed_JoinMap.recode.vcf -Ov -o Mitros_C659_Chr10_removed_JoinMap.recode_noGQ.vcf
```

now the python script works well:
```
/home/ben/projects/rrg-ben/ben/2022_GBS_lotsofxennies/individual_gvcfs_by_species/2023_Mitros_trop/JoinMap_C659/vcf2loc/vcf2loc.py -t CP -a SRR8704355 -b SRR8704354 -o Mitros_C659_Chr10_removed_JoinMap.recode.loc Mitros_C659_Chr10_removed_JoinMap.recode_noGQ.vcf
```

# Potential problems
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
options(scipen=999)
setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2022_GBS_lotsof_Xennies/JoinMap/C660/regkos_results_allin")
setwd("/Users/Shared/Previously Relocated Items/Security/projects/2022_GBS_lotsof_Xennies/JoinMap/C660/DP8_results_allin")
setwd("/Users/Shared/Previously Relocated Items/Security/projects/2022_GBS_lotsof_Xennies/JoinMap/XT_GW/results_maxmiss80DP8_segdist_removed")
# concatenate joinmap files from each chr like this:
# head -1 Chr1_jointmap.txt > all_joinmap.txt; awk 'FNR>1{print}' *jointmap*txt >> all_joinmap.txt
# head -1 C660_Chr1_regkos_matmap.txt > all_joinmap_regkos_mat.txt; awk 'FNR>1{print}' *regkos_matmap*txt >> all_joinmap_regkos_mat.txt
# head -1 C660_Chr1_regkos_patmap.txt > all_joinmap_regkos_pat.txt; awk 'FNR>1{print}' *regkos_patmap*txt >> all_joinmap_regkos_pat.txt
my_df_mat <- read.table("all_GW_joinmap_regkos_mat.txt", header = T, sep = "\t")
my_df_pat <- read.table("all_GW_joinmap_regkos_pat.txt", header = T, sep = "\t")
colnames(my_df_mat) <- c("SN","Nr","Locus","SegPat","Identical","Group","cM")
colnames(my_df_pat) <- c("SN","Nr","Locus","SegPat","Identical","Group","cM")
my_df_mat$matpat <- "mat"
my_df_pat$matpat <- "pat"
my_df <- rbind(my_df_mat,my_df_pat)
my_df[c('Chr', 'Coordinate')] <- str_split_fixed(my_df$Locus, '_', 2)

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
my_df$Coordinate_Mb <- my_df$Coordinate/1000000
# Calculate the derivatives for plotting first
# save the results to a big df called "master_derivative_df"

#create data frame with 0 rows and 3 columns
master_derivative_df <- data.frame(matrix(ncol = 5, nrow = 0))

#provide column names
colnames(master_derivative_df) <- c('Coordinate_Mb', 'predictions', 'df', 'matpat', 'Chr')

library(mgcv)

# Cycle through each Chr and each mat pat
for(i in unique(my_df$Chr)){
  for(j in unique(my_df$matpat)){
    print(i)
    print(j)
    # select only the Chr1
    my_df_matpat_Chr_only <- my_df[(my_df$Chr==i)&(my_df$matpat==j),]
    # Build the model
    model <- gam(cM ~ s(Coordinate_Mb), data = my_df_matpat_Chr_only)
    # add the predicted values to the df
    my_df_matpat_Chr_only <-cbind(my_df_matpat_Chr_only,model$fitted.values)
    colnames(my_df_matpat_Chr_only)[12] <- "fitted_values"
    # Make predictions for length of chr every 5Mb
    lengths <- as.data.frame(seq(as.integer(min(my_df_matpat_Chr_only$Coordinate_Mb)),as.integer(max(my_df_matpat_Chr_only$Coordinate_Mb)),5))
    colnames(lengths) <- "Coordinate_Mb"
    predictions <- model %>% predict(lengths)
    # now make a new df with the 5Mb coordinate increments and also the predicted values
    new_lengths <- cbind(lengths,predictions)
    
    ## now evaluate derivatives of smooths with associated standard 
    ## errors, by finite differencing...
    X0 <- predict(model,as.data.frame(new_lengths),type="lpmatrix") 
    
    eps <- 1e-7 ## finite difference interval
    x.mesh <- new_lengths$Coordinate_Mb + eps ## shift the evaluation mesh
    newd <- data.frame(Coordinate_Mb = x.mesh)
    X1 <- predict(model,newd,type="lpmatrix")
    
    Xp <- (X1-X0)/eps ## maps coefficients to (fd approx.) derivatives
    #colnames(Xp)      ## can check which cols relate to which smooth
    
    #par(mfrow=c(2,2))
    ## calculate derivatives and corresponding CIs
    Xi <- Xp*0 
    Xi[,1:9+1] <- Xp[,1:9+1] ## Xi%*%coef(b) = smooth deriv i
    # df has the derivative for each value in "x.mesh"
    df <- Xi%*%coef(model)              ## ith smooth derivative 
    df.sd <- rowSums(Xi%*%model$Vp*Xi)^.5 ## cheap diag(Xi%*%b$Vp%*%t(Xi))^.5
    #plot(x.mesh,df,type="l",ylim=range(c(df+2*df.sd,df-2*df.sd)))
    #lines(x.mesh,df+2*df.sd,lty=2);lines(x.mesh,df-2*df.sd,lty=2)
    derivative_df <-cbind(new_lengths,df)
    derivative_df$matpat <- j
    derivative_df$Chr <- i
    master_derivative_df <- rbind(master_derivative_df,derivative_df)
  }
}  



png(filename = "Recombination_plot.png",w=1000, h=200,units = "px", bg="transparent")
p<-ggplot(my_df %>% arrange(Chr), aes(x=Coordinate/1000000, y=cM, col = matpat)) + 
  scale_color_manual(breaks = c("mat", "pat"), values=c("red","blue"), labels=c('Maternal','Paternal')) +
  geom_point(size=0.5) +
  geom_smooth(method = "scam", 
              # b-spline monotonic increase
              # see ?shape.constrained.smooth.terms
              formula = y ~ s(x, k = 10, bs = "mpi"), 
              se = FALSE) +
  #geom_line() +
  #smooth.spline(x=as.numeric(Coordinate)/1000000) +
  #stat_smooth(method = smP, se= FALSE, colour='red', formula = y~x) +
  #stat_spline() +
  #geom_spline(aes(x =Coordinate/1000000, y = cM, colour = matpat), nknots = 5) +
  #stat_smooth(method = lm, formula = y ~ splines::bs(x, df = 2)) +
  # stat_smooth(method = "lm", formula = y~ns(x,knots = seq(0,250000000,500000)),  lty = 1, col = "red", se=FALSE) +
  #geom_line(aes(colour = "red"), linetype = 1) +
  #geom_line(aes(my_df$Coordinate,my_df$cM), data=data.frame(spline(my_df$Coordinate,my_df$cM, n=100))) +
  scale_y_continuous(name="Map Units (cM)", limits=c(0,260), breaks=c(0,50,100,150)) +
  scale_x_continuous(name="Coordinates (Mb)", breaks=c(0,50,100,150,200,250)) +
  #geom_hline(yintercept=0) +
  #geom_hline(yintercept=c(-0.5,0.5), linetype='dashed', color=c('black', 'black'))+
  # get rid of gray background
  facet_grid(~factor(Chr),scales="free_x",space = "free_x") +
  geom_vline(data=filter(my_df, Chr=="Chr1"), aes(xintercept=89.2370965), colour="black") + 
  geom_vline(data=filter(my_df, Chr=="Chr2"), aes(xintercept=67.5106265), colour="black") + 
  geom_vline(data=filter(my_df, Chr=="Chr3"), aes(xintercept=16.750087), colour="black") + 
  geom_vline(data=filter(my_df, Chr=="Chr4"), aes(xintercept=46.5960065), colour="black") + 
  geom_vline(data=filter(my_df, Chr=="Chr5"), aes(xintercept=62.0151615), colour="black") + 
  geom_vline(data=filter(my_df, Chr=="Chr6"), aes(xintercept=73.091979), colour="black") + 
  geom_vline(data=filter(my_df, Chr=="Chr7"), aes(xintercept=60.385499), colour="black") + 
  geom_vline(data=filter(my_df, Chr=="Chr8"), aes(xintercept=21.4457195), colour="black") + 
  geom_vline(data=filter(my_df, Chr=="Chr9"), aes(xintercept=42.124650), colour="black") + 
  geom_vline(data=filter(my_df, Chr=="Chr10"), aes(xintercept=21.2255995), colour="black") + 
  theme_classic() +
  expand_limits(x = 0) +
  theme(legend.title=element_blank()) +
  theme(text = element_text(size = 12)) # +
  # guides(color = guide_legend(override.aes = list(size = 5)))
# Get rid of the legend
#theme(legend.position = "none")
p 
dev.off()


# Now plot derivatives
master_derivative_df$Chr <- factor(master_derivative_df$Chr,
                    levels = c("Chr1", "Chr2","Chr3",
                               "Chr4","Chr5","Chr6",
                               "Chr7","Chr8","Chr9",
                               "Chr10"), ordered = T)

png(filename = "Derivative_plot.png",w=1000, h=200,units = "px", bg="transparent")
p<-ggplot(master_derivative_df %>% arrange(Chr), aes(x=Coordinate_Mb, y=df, col = matpat)) + 
  scale_color_manual(breaks = c("mat", "pat"), values=c("red","blue"), labels=c('Maternal','Paternal')) +
  #geom_point(size=0.5) +
  geom_line() +
  scale_y_continuous(name="Recombination rate (cM/Mb)", limits=c(-1,10), breaks=seq(0,10,2)) +
  scale_x_continuous(name="Coordinates (Mb)", breaks=c(0,50,100,150,200,250)) +
  facet_grid(~factor(Chr),scales="free_x",space = "free_x") +
  geom_vline(data=filter(my_df, Chr=="Chr1"), aes(xintercept=89.2370965), colour="black") + 
  geom_vline(data=filter(my_df, Chr=="Chr2"), aes(xintercept=67.5106265), colour="black") + 
  geom_vline(data=filter(my_df, Chr=="Chr3"), aes(xintercept=16.750087), colour="black") + 
  geom_vline(data=filter(my_df, Chr=="Chr4"), aes(xintercept=46.5960065), colour="black") + 
  geom_vline(data=filter(my_df, Chr=="Chr5"), aes(xintercept=62.0151615), colour="black") + 
  geom_vline(data=filter(my_df, Chr=="Chr6"), aes(xintercept=73.091979), colour="black") + 
  geom_vline(data=filter(my_df, Chr=="Chr7"), aes(xintercept=60.385499), colour="black") + 
  geom_vline(data=filter(my_df, Chr=="Chr8"), aes(xintercept=21.4457195), colour="black") + 
  geom_vline(data=filter(my_df, Chr=="Chr9"), aes(xintercept=42.124650), colour="black") + 
  geom_vline(data=filter(my_df, Chr=="Chr10"), aes(xintercept=21.2255995), colour="black") + 
  theme_classic() +
  expand_limits(x = 0) +
  theme(legend.title=element_blank()) +
  theme(text = element_text(size = 12)) # +
p 
dev.off()

# get mat and pat total length

mat_cM_length <- 0
pat_cM_length <- 0
mat_Coordinate_length <- 0
pat_Coordinate_length <- 0

chrs <- c("Chr1", "Chr2","Chr3",
          "Chr4","Chr5","Chr6",
          "Chr7","Chr8","Chr9",
          "Chr10")
for(i in 1:length(chrs)){
  mat_cM_length <- mat_cM_length + max(my_df[(my_df$matpat == "mat")&(my_df$Chr == chrs[i]),]$cM)
  pat_cM_length <- pat_cM_length + max(my_df[(my_df$matpat == "pat")&(my_df$Chr == chrs[i]),]$cM)
  mat_Coordinate_length <- mat_Coordinate_length + max(my_df[(my_df$matpat == "mat")&(my_df$Chr == chrs[i]),]$Coordinate)-
                                                  min(my_df[(my_df$matpat == "mat")&(my_df$Chr == chrs[i]),]$Coordinate)
  pat_Coordinate_length <- pat_Coordinate_length + max(my_df[(my_df$matpat == "pat")&(my_df$Chr == chrs[i]),]$Coordinate)-
                                                  min(my_df[(my_df$matpat == "pat")&(my_df$Chr == chrs[i]),]$Coordinate)
}

mat_cM_length; pat_cM_length
mat_Coordinate_length; pat_Coordinate_length

mat_cM_length/mat_Coordinate_length*1000000
pat_cM_length/pat_Coordinate_length*1000000
# get the sum of 

# https://stackoverflow.com/questions/6356665/how-do-i-plot-the-first-derivative-of-the-smoothing-function
require(splines) #thx @Chase for the notice
model <- lm(as.numeric(Position)~as.numeric(Pos),data=my_df)

dY <- diff(as.numeric(Pos))/diff(as.numeric(Position))  # the derivative of your function
dX <- rowMeans(embed(X$t,2)) # centers the X values for plotting
plot(dX,dY,type="l",main="Derivative") #check

# https://stats.stackexchange.com/questions/147733/equation-of-a-fitted-smooth-spline-and-its-analytical-derivative
# subset Chr1 only
Chr1_only <- my_df[my_df$Chr == "Chr1",];Chr1_only
spline <- interpSpline(as.numeric(Chr1_only$Pos),as.numeric(Chr1_only$Position))
plot(spline)
points(x,y)
```

