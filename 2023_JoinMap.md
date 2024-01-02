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
cut -d ' ' -f1 Mitros_C659_Chr10_removed_JoinMap.recode.loc
```

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
* Click on the Yellow "Calculation Options" icon and select "Show weak linkages with a rec. freq. larger than:" and then click "Save to Project"
* Now click on the yellow population node again and select the "Grounpings (tree)" pane in the right pane
* Click on the calculate icon in the icon toolbar; this will generate a tree with grouped markers with different stringencies (based on LOD scores but this is adjustable)
* Right click on a node in the Grouping tree (s) that you want to focus on. This should be a node with lots of markers in it and a reasonably high LOD score. Use some criteria such as LOD to select notes
* In the Population Menu, select "Create Groups Using the Grouping Tree". This creates another node in the leftmost pane called "Grouping 1" which has one or more groups within it
* Within a "Grouping" there are one or more "Groups". You can right click on one of them and a tab appears on the right that has "Fixed Orders". Here you can load the orders beginning with an @ sign
* You probably will have to delete markers that are identical or redundant to another marker; hopefully there is a way to automate this?
* The if you click on a group you can go to the "Group" menu option and click "Calculate Map"
* This generates a "Mapping" icon within the "Group" that has squiggly yellow pattern. Within this there is a map icon (purple) that contains the joint map plus one for each parent (P1, P2, green icons)
* If you click on one of the green icons, you can then select the "Map" tab and see the locus names and the position in cM
* After right clicking and selecting all of the rows, you can export this using the "Edit" menu and select "Export to File"

* Below is in progress
* Once there are maps for multiple groups (e.g. the mat and pat sites) you can right click on each group (a square yellow icon with 9 dots in it) and then go to the "Join" menue and select "Combine Groups for Map Integration". I think this makes a joint matpat map.

  

