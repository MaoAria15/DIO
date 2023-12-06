#! /	usr/bin/perl
use strict;
use warnings;
use utf8;
use Data::Dumper qw(Dumper);
use FindBin qw($Bin);
use lib ("$Bin/models");
use Table;
use BlastTable;

#### Input Files
my $selectedListFileName = "data/prediction_plaque_10006.txt";
my $blastTableFileName = "data/phage_blastx_uniprot_plaque_10006.txt";

### Output Files
my $blastResultFileName = "output/blast_result.csv";
open(my $blastResultFile, '>', $blastResultFileName);
my $blastFilteredFileName = "output/blast_filtered.csv";
open(my $blastFilteredFile, '>', $blastFilteredFileName);
print $blastFilteredFile "Name,Reason\n";

### Load Blast Table
my $blastTable = BlastTable->new({
	'filePath' => $blastTableFileName,
	'headerRowIndex' => 3,
	'headerRegex' => qr/,/
});
### Load Selected List file
my $selectedListTable = Table->new({
	'filePath' => $selectedListFileName,
	'headerRowIndex' => 0,
	'headerRegex' => qr/,/
});
my @deletedIDs = $blastTable->deleteEverythingExceptIDs(keys %{$selectedListTable->data});
for my $id (@deletedIDs) {
	print $blastFilteredFile $id.",did not match the selected list\n";
}

### Filter by text
print "Filter by: ";
my $blastFilter = <>;
chomp $blastFilter;

@deletedIDs = $blastTable->filterByField('subject acc.ver', $blastFilter);
for my $id (@deletedIDs) {
	print $blastFilteredFile $id.",did not match the filter '".$blastFilter."'\n";
}

### Column filters
my $eValueMax = 1e-3;
my $alignmentMin = 50;
@deletedIDs = $blastTable->deleteRowsWhereFieldIs('evalue'," > ".$eValueMax);
for my $id (@deletedIDs) {
	print $blastFilteredFile $id.",evalue was > ".$eValueMax."\n";
}
@deletedIDs = $blastTable->deleteRowsWhereFieldIs('alignment length'," < ".$alignmentMin);
for my $id (@deletedIDs) {
	print $blastFilteredFile $id.",alignment length was < ".$alignmentMin."\n";
}

### Write Results
print $blastResultFile join(',', @{$blastTable->header})."\n";

for my $row ($blastTable->rows) {
	my $count = 0;
	for my $column (@{$blastTable->header}){
		print $blastResultFile $row->{$column};
		if($count < scalar(@{$blastTable->header})-1) {
			print $blastResultFile ",";
		}
		$count++;
	}
	print $blastResultFile "\n";
}
close $blastFilteredFile;
close $blastResultFile;

