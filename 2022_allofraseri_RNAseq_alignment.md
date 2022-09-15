# Alignment of RNAseq data to XL genome using STAR

I'm interested in designing primers for X. allofraseri (and other species) and want to take advantage of the RNAseq data we generated for X. allofraseri many years ago.

I'm using STAR, which is a splice aware aligner to align it to XL v 9.2, which is the same version that I aligned the allofraseri GBS data to.

First I need to format the genome for STAR. Scripts are here on graham:

```
/home/ben/projects/rrg-ben/ben/2021_Austin_XB_genome/ben_scripts
```
THis is the script: 2022_STAR_index_genome.sh
```
#!/bin/sh
#SBATCH --job-name=STAR_index
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=256gb
#SBATCH --output=STAR_index.%J.out
#SBATCH --error=STAR_index.%J.err
#SBATCH --account=def-ben

module load StdEnv/2020 star/2.7.9a

STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir /home/ben/projects/rrg-ben/ben/2020_XL_v9.2_refgenome/ \
--genomeFastaFiles /home/ben/projects/rrg-ben/ben/2020_XL_v9.2_refgenome/XENLA_9.2_genome.fa \
--sjdbGTFfile /home/ben/projects/rrg-ben/ben/2020_XL_v9.2_refgenome/XLv9.2_xenbase_annotations.gff \
--sjdbOverhang 99 \
--limitGenomeGenerateRAM=124544990592
```

# map the reads using STAR, which is splice aware (2022_STAR_map_reads.sh)
```
#!/bin/sh
#SBATCH --job-name=STAR_map
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=1
#SBATCH --time=4:00:00
#SBATCH --mem=64gb
#SBATCH --output=STAR_map.%J.out
#SBATCH --error=STAR_map.%J.err
#SBATCH --account=def-ben

module load StdEnv/2020 star/2.7.9a

STAR --genomeDir /home/ben/projects/rrg-ben/ben/2020_XL_v9.2_refgenome/ \
--runThreadN 6 \
--readFilesIn ${1} ${2} \
--outFileNamePrefix ${3} \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--readFilesCommand zcat
```

# Need to fix the cigars in the bam file using GATK (2021_SplitNCigarReads.sh):
```
#!/bin/sh
#SBATCH --job-name=SplitNCigarReads
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=96:00:00
#SBATCH --mem=32gb
#SBATCH --output=SplitNCigarReads.%J.out
#SBATCH --error=SplitNCigarReads.%J.err
#SBATCH --account=def-ben

# sbatch 2021_SplitNCigarReads.sh ref bamfile

module load nixpkgs/16.09 gatk/4.1.0.0

commandline="gatk --java-options -Xmx24G SplitNCigarReads -R ${1} -I ${2} -O ${2}_SplitNCigarReads.bam"
${commandline} 
#${commandline}
```

# Calling genotypes with haplotypecaller

I am setting the threshold for map quality to be very low because I saw different thresholds for homozgous ref and alt genotypes in this paper, Fig. 3: https://acsess.onlinelibrary.wiley.com/doi/10.3835/plantgenome2019.01.0002
 https://doi.org/10.3835/plantgenome2019.01.0002
```
#!/bin/sh
#SBATCH --job-name=HaplotypeCaller
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=96:00:00
#SBATCH --mem=30gb
#SBATCH --output=HaplotypeCaller.%J.out
#SBATCH --error=HaplotypeCaller.%J.err
#SBATCH --account=def-ben


# This script will read in the *_sorted.bam file names in a directory, and 
# make and execute the GATK command "RealignerTargetCreator" on these files. 

# execute like this:
# sbatch 2021_HaplotypeCaller.sh ref bam chr
# sbatch 2021_HaplotypeCaller.sh /home/ben/projects/rrg-ben/ben/2020_XL_v9.2_refgenome/XENLA_9.2_genome.fa bam chr

module load nixpkgs/16.09 gatk/4.1.0.0

gatk --java-options -Xmx24G HaplotypeCaller  -I ${2} -R ${1} -L ${3} -O ${2}_${3}.g.vcf -ERC BP_RESOLUTION --output-mode EMIT_ALL
_CONFIDENT_SITES --read-filter MappingQualityReadFilter --minimum-mapping-quality 5
```

# Now genotype with GenotypeGVCFs
```
sbatch 2021_GenotypeGVCFs_onesample_onechr.sh ../../2020_XL_v9.2_refgenome/XENLA_9.2_genome.fa ../../2022_GBS_lotsofxennies/allofraseri_RNAseq/BJE3489_allofraseri_RNAseq_to_XLAligned_sorted_rg.bam_SplitNCigarReads.bam_chr7L.g.vcf chr7L
```

I set a very low threshold for calling genotypes (lower than the default because of issues discussed above:
```
#!/bin/sh
#SBATCH --job-name=GenotypeGVCFs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=96:00:00
#SBATCH --mem=32gb
#SBATCH --output=GenotypeGVCFs.%J.out
#SBATCH --error=GenotypeGVCFs.%J.err
#SBATCH --account=def-ben


# This script will read in the *.g.vcf file names in a directory, and 
# make and execute the GATK command "GenotypeGVCFs" on these files. 

# execute like this:
# sbatch 2021_GenotypeGVCFs.sh /home/ben/projects/rrg-ben/ben/2020_XL_v9.2_refgenome/XENLA_9.2_genome.fa /home/ben/projects/rrg-ben/ben
/2020_GBS_muel_fish_allo_cliv_laev/raw_data/cutaddapted_by_species_across_three_plates/clivii/ chr1L

module load nixpkgs/16.09 gatk/4.1.0.0

commandline="gatk --java-options -Xmx24G GenotypeGVCFs -R ${1} -V ${2} -L ${3} -all-sites --standard-min-confidence-threshold-for-calli
ng 5 --output ${2}_${3}_allsites_thresh_5.vcf"
${commandline}
```
