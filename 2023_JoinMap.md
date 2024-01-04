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
vcftools --vcf XX_Chr1_removed.vcf --max-missing 0.8 --min-meanDP 10 --recode --out XX_Chr1_removed_JoinMap
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
./vcf2loc.py -t CP -a SRR8704355 -b SRR8704354 -o Mitros_C659_Chr10_removed_JoinMap.recode.loc Mitros_C659_Chr10_removed_JoinMap.recode_noGQ.vcf
```
# Fixed sites
Use cut to get a list of the fixed sites:
```
cut -d ' ' -f1 Mitros_C659_Chr10_removed_JoinMap.recode.loc > C659_Chr10_fixed.txt
```
This file then needs to be edited by removing some text in the beginning and end and replacing hard returns with spaces.

# JoinMap

The loc file can now be opened with JoinMap. 

* Go to the file menu and save a now project
* Load the loc file from the File menu "Load Data"
* After right clicking the yellow square node in the left most pane, in the "Dataset" Menu you need to select "Create New Dataset from Data Tabsheet"
* Now the data are loaded and you can "Check for Coding Errors" in the Dataset menu.
* Now right click on the yellow square icon and then select the "Locus Genot. Freq." tab and then click on the calculator icon in the toolbar below the part of the menu with words
* This will test for segregation distortion and the results (X2 values) are in the "Locus Genot Freq" pane in a column called "X2"; sort this by p value by clicking the "Signif" tab twice
* Highlight the significant ones by right clicking. I've been excluding the ones with 4 or more asterisks. Go to the "Population" menu and select "Exclude Marked Items". This will delete the ones with segregation distortion.
* You can confirm that they are excluded by clicking on the "Loci" tab and checking some of them - they should have the "exclude" checkbox selected
* Now click on the yellow population node again and select the "Groupings (tree)" pane in the right pane
* Click on the calculate icon in the icon toolbar; this will generate a tree with grouped markers with different stringencies (based on LOD scores but this is adjustable)
* Right click on a node in the Grouping tree (s) that you want to focus on. This should be a node with lots of markers in it and with a LOD score of at least 4.
* In the Population Menu, select "Create Groups Using the Grouping Tree". This creates another node in the leftmost pane called "Grouping 1" which has one or more groups within it
* Within a "Grouping" there are one or more "Groups" (hopefully only one). If you click on a Group and select the "Loci" tab you can exclude identical individuals. Start by excluding the second identical one and keeping the first.
* Now click on the "Fixed Orders" tab. Here you can load the orders beginning with an @ sign
* Then if you click on a group you can go to the "Group" menu option and click "Calculate Map"
* This generates a "Mapping" icon within the "Group" that has squiggly yellow pattern. Within this there is a map icon (purple) that contains the joint map plus one for each parent (1, 1_P1, 1_P2, green icons)
* Once this first round of mapping is done, you can click on the purple map icon that is the parent of the join and parental maps. Check the stress for each locus and exclude ones that are >100 by clicking on the yellow Group icon and excluding markers in the loci tab. Repeat the map with the reduced number of markers by going to the "Group" menu option and clicking "Calculate Map". After what may be several iterations of this, also check that the order of the markers in the map is ascending; exclude any that are out of order and redo the map until a final map is produced with sequential markers and low stress for all markers.
* After right clicking and selecting all of the rows, you can export this using the "Edit" menu and select "Export to File"

  

