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
sed -i -e 's/__sorted.bam//g' Mitros_C659_Chr10_removed_JoinMap.recode_noGQ.vcf
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

# JoinMap

The loc file can now be opened with JoinMap. After selecting the node in the left most pane, in the "Dataset" Menu you need to select "Creat New Dataset from Dada Tabsheet"

* Now the data are loaded and you can "Check for Coding Errors" in the Dataset menu.
* Now click on the "Create Population Node" option in the Dataset Menu
* With the newly created node selected, click on the calculator icon in the toolbar below the part of the menu with words
* This will test for segregation distortion and the results (X2 values) are in the "Locus Genot Freq" pane in a column called "X2"
* *** I need to figure out how to delete the loci with segregation distortion
* Now click on the population node again and select the "Grounpings (tree)" pane in the right pane
* Click on the calculate icon in the icon toolbar; this will generate a tree with grouped markers with different stringencies (based on LOD scores but this is adjustable)
* Right click on a node in the tree that you want to focus on. This should be a node with lots of markers in it and a reasonably high LOD score. Also right click on other nodes with the same level of stringency
* In the Population Menu, select "Create Groups Using the Grouping Tree". This creates another node in the leftmost pane called "Grouping 1" which has one or more groups within it
* The if you click on a group you can go to the "Group" menu option and click "Calculate Map"

  

