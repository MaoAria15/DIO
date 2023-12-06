#! /	usr/bin/perl
use strict;
use warnings;
use utf8;
use Data::Dumper qw(Dumper);
use FindBin qw($Bin);
use lib ("$Bin/models");
use Table;

### Configuration
my $samplePrefix = qr/^meta/;
my $sampleSignificance = 0.98;
print "Please enter your desired minimum sample size: ";
my $minSampleSize = <> + 0;
print "Please enter the minimum fraction (0.0 -> 1.0) an OTU needs to be included: ";
my $sampleMinimumFraction = <> + 0.0;

#### Input Files
my $contigsFileName = "data/contigs_sizes.txt";
my $otuTableFileName = "data/h_qual_otu_table_normalized_FVT_cecum.txt";
### Output Files
my $newOtuTableFileName = "output/otu_result.tsv";
open(my $newOtuTableFile, '>', $newOtuTableFileName);
my $filteredOtusFileName = "output/otu_filtered.tsv";
open(my $filteredOtusFile, '>', $filteredOtusFileName);
print $filteredOtusFile "OTU_ID,Filtered because\n";

### Load OTU Table
my $otuTable = Table->new({
	'filePath' => $otuTableFileName,
	'headerRowIndex' => 0
});
### Load contigs file
my $contigTable = Table->new({
	'filePath' => $contigsFileName, 
	'headerRowIndex' => 0
});
# Filter by sample size
my @deletedIDs = $contigTable->deleteRowsWhereFieldIs('Size', "< ".$minSampleSize);
for my $id (@deletedIDs) {
	print $filteredOtusFile $id.",lower than ".$minSampleSize." samples\n";
	$otuTable->deleteID($id);
}

### Identify Samples Names
my @samples;
for my $fieldName (@{$otuTable->header}) {
	if($fieldName =~ $samplePrefix) { # Samples start with "truncated"
		push(@samples, $fieldName);
	}
}

### Calculate Sums of Samples
my %sampleSums;
for my $sample (@samples) {
	$sampleSums{$sample} = $otuTable->sumForField($sample);
}

### Identify 95% significance.
my %insignificance;
for my $sample (@samples) {
	# Sort OTU IDs by the sample size
	my @order = sort { $b->{$sample} <=> $a->{$sample} } $otuTable->rows;
	my $sum = 0.0;
	for my $otu (@order) {
		$sum += ($otu->{$sample}+0.0);

		if( (($sum / $sampleSums{$sample}) >= $sampleSignificance)) {
			if(!exists $insignificance{$otu->{$otuTable->idKey}}){
				$insignificance{$otu->{$otuTable->idKey}} = 0;
			}
			$insignificance{$otu->{$otuTable->idKey}} += 1;
		}
	}
}
for my $id (keys %insignificance) {
	if(($insignificance{$id} / scalar(@samples)) > $sampleMinimumFraction) {
		print $filteredOtusFile $id.",sample was not significant in ".$insignificance{$id}." samples (".($insignificance{$id} / scalar(@samples)).") > ".$sampleMinimumFraction."\n";		
		$otuTable->deleteID($id);
	}
}
### Write Results
print $newOtuTableFile join("\t", @{$otuTable->header})."\n";
for my $row ($otuTable->rows) {
	my $count = 0;
	for my $column (@{$otuTable->header}){
		print $newOtuTableFile $row->{$column};
		if($count < scalar(@{$otuTable->header})-1) {
			print $newOtuTableFile "\t";
		}
		$count++;
	}
	print $newOtuTableFile "\n";
}
close $newOtuTableFile;
close $filteredOtusFile;

