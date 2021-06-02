# 2020_GBS.md
Raw data are here:
```
/home/ben/projects/rrg-ben/ben/2020_GBS_muel_fish_allo_cliv_laev/raw_data/plate1/raw_data_not_demultiplexed
/home/ben/projects/rrg-ben/ben/2020_GBS_muel_fish_allo_cliv_laev/raw_data/plate2/raw_data_not_demultiplexed
/home/ben/projects/rrg-ben/ben/2020_GBS_muel_fish_allo_cliv_laev/raw_data/plate3/raw_data_not_demultiplexed
```

Raw data were demultiplexed using sabre:
```
#!/bin/sh
#SBATCH --job-name=sabre_plate1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=2gb
#SBATCH --output=sabre_plate1.%J.out
#SBATCH --error=sabre_plate1.%J.err
#SBATCH --account=def-ben


module load nixpkgs/16.09  intel/2016.4
module load module load sabre/1.00

sabre pe -f NS.1413.003.D701.Xenopus_Plate1_R1.fastq.gz -r NS.1413.003.D701.Xenopus_Plate1_R2.fastq.gz -b sabre_barcodes.txt -u
 no_bc_match_R1.fq -w no_bc_match_R2.fq
 ```
 
 Because the adapters are incorporated into the 3' end I used cutadapt to trim adapter seqs:
 ```
 #!/bin/sh
#SBATCH --job-name=cutadapt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=36:00:00
#SBATCH --mem=2gb
#SBATCH --output=cutadapt.%J.out
#SBATCH --error=cutadapt.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this
# sbatch ./2020_cutadapt.sh ../raw_data/plate1


module load python/3.7

# this allows for 20% mismatch of adapters, and trims from the 5' and 3'
# ends based on quality scores. the -B flag tells cutadapt to trim
# adaptors anywhere in the read.  The adapter sequences are the same for
# both directions because the end of the primer that would be
# incorporated into the seq in the 3' end of each read is identical.

#~/.local/bin/cutadapt -a "AGATCGGAAGAGC;max_error_rate=0.2" -A "AGATCGGAAGAGC;max_error_rate=0.2" -B -q 15,10 -o out.1.fastq -
p out.2.fastq reads.1.fastq reads.2.fastq


v=1
#  Always use for-loop, prefix glob, check if exists file.
for file in $1/*.fq.gz ; do         # Use ./* ... NEVER bare *
  if [ -e "$file" ] ; then   # Check whether file exists.
  	if [[ $v -eq 1 ]]
	then # if/then branch
	    ~/.local/bin/cutadapt -b "AGATCGGAAGAGC" -B "AGATCGGAAGAGC" -e 0.2 -q 15,10 -o ${file::-8}cut.R1.fastq -p ${file::-
8}cut.R2.fastq ${file::-8}R1.fq.gz ${file::-8}R2.fq.gz
		  v=0
	else # else branch
  		v=1
	fi
  fi
done 
```

And then I used trimmomatic to filter out low quality reads (but no adapter clipping):
```
#!/bin/sh
#SBATCH --job-name=trimmomatic
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=8gb
#SBATCH --output=trimmomatic.%J.out
#SBATCH --error=trimmomatic.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this
# sbatch ./2020_trimmomatic.sh ../raw_data/plate1


module load StdEnv/2020
module load trimmomatic/0.39

v=1
#  Always use for-loop, prefix glob, check if exists file.
for file in $1/*_cut.R1.fastq.gz ; do         # Use ./* ... NEVER bare *
  if [ -e "$file" ] ; then   # Check whether file exists.
  	if [[ $v -eq 1 ]]
	then # if/then branch
		  java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE ${file::-16}_cut.R1.fastq.gz ${file::-16}_cut.R2.fastq.g
z ${file::-16}_cuttrim.R1.fq.gz ${file::-16}_cuttrim.R2.fq.gz SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:3
		  v=0
	else # else branch
  		v=1
	fi
  fi
done 
```
To check for repetitive seqs that could (somehow) arise from sneaky adapters, I ran fastqc after trimmomatic:
```
#!/bin/sh
#SBATCH --job-name=fastqc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=36:00:00
#SBATCH --mem=2gb
#SBATCH --output=fastqc.%J.out
#SBATCH --error=fastqc.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this
# sbatch ./2020_fastqc.sh ../raw_data/plate1

# first make sure the files are all compressed

for file in $1/*_cut.R2.fastq ; do         # Use ./* ... NEVER bare *  
  if [ -e "$file" ] ; then   # Check whether file exists.                  
      echo $file # print the file name
      gzip $file # this will delete the original file once it is compressed
  fi
done

# now run fastqc on the compressed files
module load fastqc/0.11.9
for file in $1/*.fq.gz ; do 
  if [ -e "$file" ] ; then  # Check whether file exists.
      fastqc $file
  fi
done
```


I renamed and concatenated data from the same sample and organized them by species here:
```
/home/ben/projects/rrg-ben/ben/2020_GBS_muel_fish_allo_cliv_laev/raw_data/cutaddapted_by_species_across_three_plates

```
add readgroups
```
#!/bin/sh
#SBATCH --job-name=readgroups
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=8gb
#SBATCH --output=readgroups.%J.out
#SBATCH --error=readgroups.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this
# sbatch ./2021_picard_add_read_groups.sh /home/ben/projects/rrg-ben/ben/2020_GBS_muel_fish_allo_cliv_laev/raw_data/cutadd
apted_by_species_across_three_plates/clivii/ 

module load picard/2.23.3

for file in ${1}*_sorted.bam
do
    java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=${file} O=${file}_rg.bam RGID=4 RGLB=$(basename $file) RGP
L=ILLUMINA RGPU=$(basename $file) RGSM=$(basename $file)
done
```

Indel realignment is not longer recommended for GATK version 4.  Instead go directly to HaplotypeCaller and output g.vcfs
```
#!/bin/sh
#SBATCH --job-name=HaplotypeCaller
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --mem=32gb
#SBATCH --output=HaplotypeCaller.%J.out
#SBATCH --error=HaplotypeCaller.%J.err
#SBATCH --account=def-ben


# This script will read in the *_sorted.bam file names in a directory, and 
# make and execute the GATK command "RealignerTargetCreator" on these files. 

# execute like this:
# sbatch 2021_HaplotypeCaller.sh /home/ben/projects/rrg-ben/ben/2020_XL_v9.2_refgenome/XENLA_9.2_genome.fa /hom
e/ben/projects/rrg-ben/ben/2020_GBS_muel_fish_allo_cliv_laev/raw_data/cutaddapted_by_species_across_three_plate
s/clivii/ 

module load nixpkgs/16.09 gatk/4.1.0.0

for file in ${2}*_sorted.bam_rg.bam
do
    gatk --java-options -Xmx24G HaplotypeCaller  -I ${file} -R ${1} -O ${file}.g.vcf -ERC GVCF
done
```

