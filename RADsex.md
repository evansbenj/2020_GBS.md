# RADsex

I'm first going to try this approach on the raw reads: https://github.com/SexGenomicsToolkit/radsex

Issues to keep in mind are that the method really is aimed at single read data.  I can probably get around this either by concatenating forward and reverse reads just for this analysis or (and this is preferable) making the table for everything in the first step and then doing the second (distrib) step using different popmap files (one for the forward reads and one for the reverse).

RADsex is installed here:
```
/home/ben/projects/rrg-ben/ben/2023_RADsex/radsex/bin/
```

First step is to make a table of reads.  This uses up some memory:
```
#!/bin/sh
#SBATCH --job-name=radsex_tableofreads
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=6:00:00
#SBATCH --mem=128gb
#SBATCH --output=radsex_tableofreads.%J.out
#SBATCH --error=radsex_tableofreads.%J.err
#SBATCH --account=rrg-ben

# run by passing an argument like this
# sbatch ./2023_RADsex_tableofreads.sh ./
/home/ben/projects/rrg-ben/ben/2023_RADsex/radsex/bin/radsex process --input-dir ${1}/ --output-file ${1}/markers_table.tsv --threads 16 --min-depth 1
```
Once this is done I can summarize it using a file that specifies whether the individuals are male or female. I turned off Bonferrini correction (-C) and made the pvalue more stringent (-S 0.0001). This worked for fischbergi, which has a massive sex linked region.

```
#!/bin/sh
#SBATCH --job-name=radsex_distrib
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --mem=64gb
#SBATCH --output=radsex_distrib.%J.out
#SBATCH --error=radsex_distrib.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this
# sbatch ./2023_RADsex_distrib.sh marker_table outputfile sexfile
/home/ben/projects/rrg-ben/ben/2023_RADsex/radsex/bin/radsex distrib --markers-table ${1} --output-file ${2} --popmap ${3} --min-depth 1 --groups M,F -C -S 0.0001
```

Or output a fasta file with significant reads:

```
/home/ben/projects/rrg-ben/ben/2023_RADsex/radsex/bin/radsex signif --markers-table XXXmarkers_table.tsv --output-file XXX_significant_markers.fasta --popmap XX_sex_R1R2cat --min-depth 1 --groups M,F --output-fasta -C -S 0.0001
```

And this can be done for clivii, muelleri, fischbergi, and parafraseri. 

I can print the output like this:
```R
#https://github.com/SexGenomicsToolkit/radsex


setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/submitted/XB_sex_determining_gene/RADsex")

library(sgtr)
library(ggplot2)

distrib_plot <- radsex_distrib("clivii_short_distribution.tsv",
                               groups = c("M", "F"),
                               group_labels = c("Males", "Females"),
                               title = "Distribution of markers",
                               significance_color = "green",
                               bins = c(0,1,5,10,20,35,100,1000),
                               colors = c("white", "blue"))

ggsave("Rplot_clivii_short.pdf", width = 5, height = 4, dpi = 200)
```


It is also possible to map the results and print this by chr, including the probability of association (P), and the sex bias (positive for male biased, Y-linked and negative for female biased, W-linked:
```
#!/bin/sh
#SBATCH --job-name=radsex_map
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=36:00:00
#SBATCH --mem=32gb
#SBATCH --output=radsex_map.%J.out
#SBATCH --error=radsex_map.%J.err
#SBATCH --account=def-ben

/home/ben/projects/rrg-ben/ben/2023_RADsex/radsex/bin/radsex map --markers-file markers_table.tsv --output-file cliv_significant_markers_alignment_results.tsv --popmap 2023_cliv_sex.txt --genome-file ~/projects/rrg-ben/ben/2021_XL_v10_refgenome/XENLA_10.1_genome.fa.gz --min-quality 20 --min-frequency 0.01 --min-depth 1 --groups M,F

```

Get rid of scaffolds like this:
```
sed '/Sca/d' /home/ben/projects/rrg-ben/ben/2020_radsex/bin/2017_cliv_markers_alignment_all.tsv > /home/ben/projects/rrg-ben/ben/2020_radsex/bin/2017_cliv_markers_alignment_all_chronly.tsv
```

The "radsex_map_circos" in the sgtr package is broken, so I just printed the results by chr using ggplot:
```
setwd("/Users/Shared/Previously Relocated Items/Security/projects/2020_fisch_clivii_muel_allo_RADseq/RADsex")

library(ggplot2)

data_SL<-read.table("./2020_allo_catR1R2_alignment_results_chronly.tsv",header=T)
data_SL_sorted <- data_SL[order(data_SL["Contig"],data_SL["Position"]),]
head(data_SL_sorted)
data_SL_sorted$Position=as.numeric(data_SL_sorted$Position)

sp <- ggplot(data_SL_sorted, aes(x=Position, y=-log(P))) + 
  geom_point(size=1) + 
  facet_wrap( ~ Contig, ncol=2)
sp

sp <- ggplot(data_SL_sorted, aes(x=Position, y=Bias)) + 
  # Bias is positive if male specific
  # neggative if female specific
  geom_point(size=1) + 
  facet_wrap( ~ Contig, ncol=2)
sp
```

# Parsing marker file

```perl
#!/usr/bin/env perl
use strict;
use warnings;
use lib qw(~/perl_modules);
#use List::MoreUtils qw/ uniq /;

# on computecanada:
# module load perl/5.30.2
# module load gcc/9.3.0
# cpan
# install List::MoreUtils

#  This program reads in a tab delimited marker file generated
#  by RADsex and outputs the fastaseqs that are present in some minimum
#  number of males and females

# Example for XB_WGS
# perl Parse_marker_table.pl 2018_2022_allofraseri_catR1R2_markers_table.tsv 021110111000000000111110000111111111111111111111111111111111000000000000000000 out.fasta max_n_fems min_n_males WorY

# the 0,1, 2 indicates whether an individual is (0) male, (1) female, and or (2) skipped

# the max_n_fems is the maximum number of females with the sequence

# the min_n_males is the minimum number of males with the sequence (this is set up for allofraseri which may have a Y chr)

# to search for W specific seqs, the binary matrix can be reversed (use 0s for females and 1s for males) and then
# the max_n_fems refers to max number of males and the min_n_males refers to the min number of females

# 2024 fraseri W 011110101101111110111111001001101110
# ../../../Parse_RADsex_marker_table.pl markers_table.tsv 100001010010000001000000110110010001 minF7_maxM1 1 7 W

# 2024 longipes W 1001010010000000100001000001011110
# ../../../Parse_RADsex_marker_table.pl markers_table.tsv 1001010010000000100001000001011110 minF7_maxM0 0 7 W


my $inputfile = $ARGV[0];
my $input2 = $ARGV[1];
my $outputfile1 = $ARGV[2];
my $max_n_fems = $ARGV[3];
my $min_n_males = $ARGV[4];
my $WorY = $ARGV[5];

unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}

unless (open(OUTFILE, ">$outputfile1"))  {
	print "I can\'t write to $outputfile1\n";
	exit;
}
print "Creating output file: $outputfile1\n";


my @sexes = split("",$ARGV[1]);

my $males;
my $females;
my @temp;
my $counter;
my $y;
my $x;
my $number_of_male_individuals_genotyped=0;
my $number_of_female_individuals_genotyped=0;

for ($y = 0 ; $y <= $#sexes ; $y++ ) {
	if($sexes[$y] == 0){
		$number_of_male_individuals_genotyped +=1;
	}	
}	
for ($y = 0 ; $y <= $#sexes ; $y++ ) {
	if($sexes[$y] == 1){
		$number_of_female_individuals_genotyped +=1;
	}	
}
print "This includes ",$number_of_female_individuals_genotyped," female(s) and  ", $number_of_male_individuals_genotyped," males\n";

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split /[\t\/]/,$line;
	if(($temp[0] !~ /^#Number/)&&($temp[0] !~ /^id/)){
		if($#temp ne ($#sexes+2)){
			print "The number of individuals in the input line does not match the number of individuals genotyped ",
			$temp[0],"\t",$temp[1],"\t",$#temp," ",$#sexes+2,"\n";
		}
		# count the number of males and females that have this sequence
		$males=0;
		$females=0;
		$counter=0;
		for ($y = 2 ; $y <= $#temp; $y++ ) {
			if($temp[$y] ne 0){
				if($sexes[$counter] == 0){
					$males+=1;
				}
				elsif($sexes[$counter] == 1){
					$females+=1;
				}	
			}
			$counter+=1;
		}  # end for
		# check if the observed is above the minimum number and below the max number
		if(($males >= $min_n_males)&&($females <= $max_n_fems)&&($WorY eq "Y")){
			print ">Seq_",$temp[0],"_F_",$females,"_M_",$males,"\n",$temp[1],"\n";
			print OUTFILE ">",$temp[0],"_F_",$females,"_M_",$males,"\n",$temp[1],"\n";
		}	
		elsif(($males >= $min_n_males)&&($females <= $max_n_fems)&&($WorY eq "W")){
			print ">Seq_",$temp[0],"_M_",$females,"_F_",$males,"\n",$temp[1],"\n";
			print OUTFILE ">",$temp[0],"_M_",$females,"_F_",$males,"\n",$temp[1],"\n";
		}	
	} # end if
} # end while	
close OUTFILE;

```
# Blast
```
module load StdEnv/2020  gcc/9.3.0 blast+/2.14.0
```
```
blastn -query minF11_maxM0.fasta -db ../../../../../2021_XL_v10_refgenome/XENLA_10.1_genome.fa_blastable -outfmt 6 -out minF11_maxM0_to_XLv10
```

This output has multiple entries for repetitive seqs; need to quantify only unique ones
