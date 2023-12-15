# Permutations

This is a script to do permutations that evaluate whether the locations of sex-linked regions is atypically far from the centers of chromosomes or from telomeres.

The permutations randomly place the same-sized sex-linked regions for each species on the same chromosome. If the sex-linked region overlaps the beginning or the end of the chr, it is placed on the beginning or end of the chromosome.

The test stat is (SL_center - chr_center/chr_length)^2 summed over each species

This way each chromosome is scaled to be the same length so each species counts the same.

# Script
```perl
#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;


# This program will read in three files:
# (1) the coordinates of de novo sex associated regions in several Xenopus species
# (2) the length of XL chromosomes and locations of XL centromeres
# (3) the length of XT chromosomes and locations of XL centromeres

# Using this information, two permutation tests will be performed. 
# The first uses as a test statistic the squared distance of the center of the sex-linked
# region from the telomere
# The second uses as a a test statistic the squared distance of the center of the chromosome

# the permutations will randomly place all of the sex-linked regions throughout the genome and 
# calculate these statistics 999 times.


my $inputfile1 = $ARGV[0];
my $inputfile2 = $ARGV[1];
my $inputfile3 = $ARGV[2];

my %SL_Hash;
my @names;
my @temp;
my $counter=0;
my @species;

# first open up the SL coordinate file 
unless (open DATAINPUT, $inputfile1) {
	print "Can not find the input file.\n";
	exit;
}

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split('\t',$line);
	if($temp[0] ne 'subgenus'){ 
		$SL_Hash{$temp[1]}[0] = $temp[0]; # subgenus
		$SL_Hash{$temp[1]}[1] = $temp[2]; # Chr
		$SL_Hash{$temp[1]}[2] = $temp[3]; # heterogamy
		$SL_Hash{$temp[1]}[3] = $temp[4]; # SL start
		$SL_Hash{$temp[1]}[4] = $temp[5]; # SL end
		push(@species,$temp[1]);
	}
}	
close DATAINPUT;
$counter=0;
my %XL_Hash;

# now open up the XL information
unless (open DATAINPUT2, $inputfile2) {
	print "Can not find the input file.\n";
	exit;
}
while ( my $line = <DATAINPUT2>) {
	chomp($line);
	@temp=split('\t',$line);
	if($temp[0] ne 'XL_Chr'){ 
		$XL_Hash{$temp[0]}[0] = $temp[1]; # length
		$XL_Hash{$temp[0]}[1] = $temp[2]; # cent_start
		$XL_Hash{$temp[0]}[2] = $temp[3]; # cent_end
	}
}
close DATAINPUT2;

$counter=0;
my %XT_Hash;

# now open up the XT information
unless (open DATAINPUT3, $inputfile3) {
	print "Can not find the input file.\n";
	exit;
}
while ( my $line = <DATAINPUT3>) {
	chomp($line);
	@temp=split('\t',$line);
	if($temp[0] ne 'XL_Chr'){ 
		$XT_Hash{$temp[0]}[0] = $temp[1]; # length
		$XT_Hash{$temp[0]}[1] = $temp[2]; # cent_start
		$XT_Hash{$temp[0]}[2] = $temp[3]; # cent_end
	}
}	
close DATAINPUT3;


# ok now I have all the data loaded in three hashes
# print these numbers
#foreach my $key (keys %SL_Hash) {
#foreach my $key (sort { $a <=> $b} keys %SL_Hash) {
#	print $key,"\t",$SL_Hash{$key}[0],"\t";
#	print "\n";	
#}	

#print $SL_Hash{"X. allofraseri"}[1],"\n";

######################
# Test statistic
######################
my $standardized_squared_distance_to_center_teststat;
my $standardized_squared_distance_to_centromere_teststat;
my $SL_center;
foreach my $key (keys %SL_Hash) { 
	if($SL_Hash{$key}[0] eq "Xenopus"){
		# $SL_Hash{$key}[1] # SL chromosome
		# $SL_Hash{$key}[3] SL begin
		# $SL_Hash{$key}[4] SL end
		# $XL_Hash{$SL_Hash{$key}[1]}[0] # chr length
		# $XL_Hash{$SL_Hash{$key}[1]}[1] # cent begin
		# $XL_Hash{$SL_Hash{$key}[1]}[2] # cent end
		$SL_center = ($SL_Hash{$key}[4] + $SL_Hash{$key}[3])/2;
		print $key,"\t",$SL_Hash{$key}[1] ,"\t",$SL_center,"\n";
		# check if it fits within the chr
		if(($SL_center + 0.5*($SL_Hash{$key}[4]-$SL_Hash{$key}[3]) <= $XL_Hash{$SL_Hash{$key}[1]}[0])&&
		($SL_center - 0.5*($SL_Hash{$key}[4]-$SL_Hash{$key}[3]) >= 0))
			{
				# it fits, so calculate the test statistics
				# first calculate the (distance to the center of the chr / chr length)^2
				print $key,"\t",$SL_center,"\t",($XL_Hash{$SL_Hash{$key}[1]}[0]/2),"\t",$XL_Hash{$SL_Hash{$key}[1]}[0],"\t",(($SL_center-($XL_Hash{$SL_Hash{$key}[1]}[0]/2))/$XL_Hash{$SL_Hash{$key}[1]}[0])**2,"\n";
				$standardized_squared_distance_to_center_teststat += (($SL_center-($XL_Hash{$SL_Hash{$key}[1]}[0]/2))/$XL_Hash{$SL_Hash{$key}[1]}[0])**2;
				# now calculate the (distance to the centromere of the chr / chr length)^2
				$standardized_squared_distance_to_centromere_teststat += (($SL_center-(($XL_Hash{$SL_Hash{$key}[1]}[2]-$XL_Hash{$SL_Hash{$key}[1]}[1])/2))/$XL_Hash{$SL_Hash{$key}[1]}[0])**2;
			}
			else{ # the SL region doesn't fit
			print "Something is wrong with how the SL regions are defined\n";
		}
	}
	elsif($SL_Hash{$key}[0] eq "Silurana"){
		$SL_center = ($SL_Hash{$key}[4] + $SL_Hash{$key}[3])/2;
		print $key,"\t",$SL_Hash{$key}[1] ,"\t",$SL_center,"\n";
		# check if it fits within the chr
		if(($SL_center + 0.5*($SL_Hash{$key}[4]-$SL_Hash{$key}[3]) <= $XT_Hash{$SL_Hash{$key}[1]}[0])&&
		($SL_center - 0.5*($SL_Hash{$key}[4]-$SL_Hash{$key}[3]) >= 0))
			{
				# it fits, so calculate the test statistics
				# first calculate the (distance to the center of the chr / chr length)^2
				$standardized_squared_distance_to_center_teststat += (($SL_center-($XT_Hash{$SL_Hash{$key}[1]}[0]/2))/$XT_Hash{$SL_Hash{$key}[1]}[0])**2;
				# now calculate the (distance to the centromere of the chr / chr length)^2
				$standardized_squared_distance_to_centromere_teststat += (($SL_center-(($XT_Hash{$SL_Hash{$key}[1]}[2]-$XT_Hash{$SL_Hash{$key}[1]}[1])/2))/$XT_Hash{$SL_Hash{$key}[1]}[0])**2;
			}
			else{ # the SL region doesn't fit
			print "Something is wrong with how the SL regions are defined\n";
		}
	}
}

print "To center test stat ",$standardized_squared_distance_to_center_teststat,"\n";
print "To telomere test stat ", $standardized_squared_distance_to_centromere_teststat,"\n";


######################
# Permutations
######################
my @standardized_squared_distance_to_center;
my @standardized_squared_distance_to_centromere;
my $perms=1000;
my $y;
my $random_number;
my $counter=0;
for ($y = 0 ; $y < $perms; $y++ ) {
	# for each species pick a random number between 0 and the chr length of the chr that is SL
	foreach my $key (keys %SL_Hash) { # each key is a species
		# pick a random number between zero and the chromosome length
		# check which subgenus
		if($SL_Hash{$key}[0] eq "Xenopus"){
			# $SL_Hash{$key}[1] # SL chromosome
			# $SL_Hash{$key}[3] SL begin
			# $SL_Hash{$key}[4] SL end
			# $XL_Hash{$SL_Hash{$key}[1]}[0] # chr length
			# $XL_Hash{$SL_Hash{$key}[1]}[1] # cent begin
			# $XL_Hash{$SL_Hash{$key}[1]}[2] # cent end
			$random_number = int(rand($XL_Hash{$SL_Hash{$key}[1]}[0]));
			print $key,"\t",$SL_Hash{$key}[1] ,"\t",$random_number,"\n";
			# check if it fits within the chr
			if(($random_number + 0.5*($SL_Hash{$key}[4]-$SL_Hash{$key}[3]) <= $XL_Hash{$SL_Hash{$key}[1]}[0])&&
			($random_number - 0.5*($SL_Hash{$key}[4]-$SL_Hash{$key}[3]) >= 0))
				{
				# it fits, so calculate the test statistics
				# first calculate the (distance to the center of the chr / chr length)^2
				$standardized_squared_distance_to_center[$counter]+=(($random_number-($XL_Hash{$SL_Hash{$key}[1]}[0]/2))/$XL_Hash{$SL_Hash{$key}[1]}[0])**2;
				# now calculate the (distance to the centromere of the chr / chr length)^2
				$standardized_squared_distance_to_centromere[$counter]+=(($random_number-(($XL_Hash{$SL_Hash{$key}[1]}[2]-$XL_Hash{$SL_Hash{$key}[1]}[1])/2))/$XL_Hash{$SL_Hash{$key}[1]}[0])**2;
				
				print "hello ",$key,"\t",$random_number,"\t",($XL_Hash{$SL_Hash{$key}[1]}[0]/2),"\t",$XL_Hash{$SL_Hash{$key}[1]}[0],"\t",(($random_number-($XL_Hash{$SL_Hash{$key}[1]}[0]/2))/$XL_Hash{$SL_Hash{$key}[1]}[0])**2,"\n";
				
				}
			else{ # the SL region doesn't fit
				# put the region at the beginning or end of the chr, depending on which one it overlaps with
				if($random_number + 0.5*($SL_Hash{$key}[4]-$SL_Hash{$key}[3]) > $XL_Hash{$SL_Hash{$key}[1]}[0]){ # it goes over the end
					# make the center of the SLregion be as far as possible towards the end
					$random_number = $XL_Hash{$SL_Hash{$key}[1]}[0] - 0.5*($SL_Hash{$key}[4]-$SL_Hash{$key}[3]);
				}
				elsif($random_number - 0.5*($SL_Hash{$key}[4]-$SL_Hash{$key}[3]) < 0){ # it goes over the beginning
					# make the center of the SLregion be as far as possible towards the beginning
					$random_number = 0.5*($SL_Hash{$key}[4]-$SL_Hash{$key}[3]);
				}
				# Now calculate the test statistics
				# first calculate the (distance to the center of the chr / chr length)^2
				$standardized_squared_distance_to_center[$counter]+=(($random_number-($XL_Hash{$SL_Hash{$key}[1]}[0]/2))/$XL_Hash{$SL_Hash{$key}[1]}[0])**2;
				# now calculate the (distance to the centromere of the chr / chr length)^2
				$standardized_squared_distance_to_centromere[$counter]+=(($random_number-(($XL_Hash{$SL_Hash{$key}[1]}[2]-$XL_Hash{$SL_Hash{$key}[1]}[1])/2))/$XL_Hash{$SL_Hash{$key}[1]}[0])**2;
			}
		}
		elsif($SL_Hash{$key}[0] eq "Silurana"){
			# $SL_Hash{$key}[1] # SL chromosome
			# $SL_Hash{$key}[3] SL begin
			# $SL_Hash{$key}[4] SL end
			# $XT_Hash{$SL_Hash{$key}[1]}[0] # chr length
			# $XT_Hash{$SL_Hash{$key}[1]}[1] # cent begin
			# $XT_Hash{$SL_Hash{$key}[1]}[2] # cent end
			$random_number = int(rand($XT_Hash{$SL_Hash{$key}[1]}[0]));
			print $key,"\t",$SL_Hash{$key}[1] ,"\t",$random_number,"\n";
			# check if it fits within the chr
			if(($random_number + 0.5*($SL_Hash{$key}[4]-$SL_Hash{$key}[3]) <= $XT_Hash{$SL_Hash{$key}[1]}[0])&&
			($random_number - 0.5*($SL_Hash{$key}[4]-$SL_Hash{$key}[3]) >= 0))
				{
				# it fits, so calculate the test statistics
				# first calculate the (distance to the center of the chr / chr length)^2
				$standardized_squared_distance_to_center[$counter]+=(($random_number-($XT_Hash{$SL_Hash{$key}[1]}[0]/2))/$XT_Hash{$SL_Hash{$key}[1]}[0])**2;
				# now calculate the (distance to the centromere of the chr / chr length)^2
				$standardized_squared_distance_to_centromere[$counter]+=(($random_number-(($XT_Hash{$SL_Hash{$key}[1]}[2]-$XT_Hash{$SL_Hash{$key}[1]}[1])/2))/$XT_Hash{$SL_Hash{$key}[1]}[0])**2;
				}
			else{ # the SL region doesn't fit
				# put the region at the beginning or end of the chr, depending on which one it overlaps with
				if($random_number + 0.5*($SL_Hash{$key}[4]-$SL_Hash{$key}[3]) > $XT_Hash{$SL_Hash{$key}[1]}[0]){ # it goes over the end
					# make the center of the SLregion be as far as possible towards the end
					$random_number = $XT_Hash{$SL_Hash{$key}[1]}[0] - 0.5*($SL_Hash{$key}[4]-$SL_Hash{$key}[3]);
				}
				elsif($random_number - 0.5*($SL_Hash{$key}[4]-$SL_Hash{$key}[3]) < 0){ # it goes over the beginning
					# make the center of the SLregion be as far as possible towards the beginning
					$random_number = 0.5*($SL_Hash{$key}[4]-$SL_Hash{$key}[3]);
				}
				# Now calculate the test statistics
				# first calculate the (distance to the center of the chr / chr length)^2
				$standardized_squared_distance_to_center[$counter]+=(($random_number-($XT_Hash{$SL_Hash{$key}[1]}[0]/2))/$XT_Hash{$SL_Hash{$key}[1]}[0])**2;
				# now calculate the (distance to the centromere of the chr / chr length)^2
				$standardized_squared_distance_to_centromere[$counter]+=(($random_number-(($XT_Hash{$SL_Hash{$key}[1]}[2]-$XT_Hash{$SL_Hash{$key}[1]}[1])/2))/$XT_Hash{$SL_Hash{$key}[1]}[0])**2;
			}
		}
		else{
			print "Subgenus not properly defined \n;"
		}
	}
	$counter+=1;
}	
print "@standardized_squared_distance_to_centromere\n";


# We have a one-sided expectation that the test statistic will be larger than the random expectation established 
# by the permutations.


my @standardized_squared_distance_to_center = sort { $a <=> $b } @standardized_squared_distance_to_center;
my @standardized_squared_distance_to_centromere = sort { $a <=> $b } @standardized_squared_distance_to_centromere;
print "@standardized_squared_distance_to_center\n";
print "@standardized_squared_distance_to_centromere\n";

my $pval=0;
for ($y = 0 ; $y <= $#standardized_squared_distance_to_center; $y++ ) {
	if($standardized_squared_distance_to_center_teststat >= $standardized_squared_distance_to_center[$y]){ 
		$pval+=1;
	}
}	

print "Center distance test stat:",$standardized_squared_distance_to_center_teststat,"\n";
print "Center P value = ",1-$pval/$perms,"\n";

my $pval=0;
for ($y = 0 ; $y <= $#standardized_squared_distance_to_centromere; $y++ ) {
	if($standardized_squared_distance_to_centromere_teststat >= $standardized_squared_distance_to_centromere[$y]){ 
		$pval+=1;
	}
}	

print "Center distance test stat:",$standardized_squared_distance_to_centromere_teststat,"\n";
print "Center P value = ",1-$pval/$perms,"\n";
```
