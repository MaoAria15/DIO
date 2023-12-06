package BlastTable;

#!/usr/bin/perl
############################################################
use strict;
use warnings;
use Data::Dumper qw(Dumper);
use List::MoreUtils qw{ any };

1; # package return

sub new {
	my $class = shift;
	my ( $args ) = @_;
	my $filePath = $args->{'filePath'};
	my $headerRowIndex = $args->{'headerRowIndex'};
	my $headerRegex = $args->{'headerRegex'} || qr/\t/;
	my $skipRowRegex = $args->{'skipRowRegex'} || qr/"?#/;

	my $self = {};
	open(my $fh, "<", $filePath)
		or die "Failed to open file: $!\n";
	my $count = 0;
	while(<$fh>) { 
		chomp; 
		my @fields;
		if($count < $headerRowIndex) {
			$count++;
		} else {
			if(!exists $self->{'header'} && $_ =~ /"?# Fields: (.+)+$/) {
				@fields = split $headerRegex, $1;
				s{^\s+|\s+$}{}g foreach @fields; # Remove leading whitespaces
				@{$self->{'header'}} = @fields;
			} else {
				if(!($skipRowRegex && $_ =~ $skipRowRegex)){
					@fields = split /\t/, $_;
					my %element;
					for my $index (0 .. $#fields) {
						my $fieldName = @{$self->{'header'}}[$index];
						$element{$fieldName} = $fields[$index];						
					}
					if(!exists $self->{'data'}{$fields[0]}) {
						$self->{'data'}{$fields[0]} = [];
					}
					push(@{$self->{'data'}{$fields[0]}}, \%element);
				}
			}	
		}
	}
	close $fh;
	bless $self, $class;
	return $self;
}

sub header {
	my $self = shift;
	return $self->{'header'};
}

sub idKey {
	my $self = shift;
	return @{$self->{'header'}}[0];
}

sub ids {
	my $self = shift;
	return keys %{$self->data};
}

sub rows {
	my $self = shift;
	return map {@$_} values %{$self->data};
}

sub data {
	my $self = shift;
	return $self->{'data'};
}

sub filterByField {
	my $self = shift;
	my $field = shift;
	my $text = shift;
	my @deletedIDs;
	for my $id ($self->ids) {
		my @removedIndices;
		for my $index (0..scalar(@{$self->data->{$id}})-1) {
			if(index(lc $self->data->{$id}[$index]{$field}, lc $text) == -1) {
				push(@removedIndices, $index);
			}
		}
		for my $index (reverse @removedIndices) {
			splice @{$self->data->{$id}}, $index, 1;
		}
		if(scalar(@{$self->data->{$id}}) == 0) {
			push(@deletedIDs, $id);
			$self->deleteID($id);
		}
	}
	return @deletedIDs;
}

sub deleteRowsWhereFieldIs {
	my $self = shift;
	my $field = shift;
	my $expression = shift;
	my @deletedIDs;
	for my $id ($self->ids) {
		my @removedIndices;
		for my $index (0..scalar(@{$self->data->{$id}})-1) {
			if(eval($self->data->{$id}[$index]{$field}." ".$expression)) {
				push(@removedIndices, $index);
			}
		}
		for my $index (reverse @removedIndices) {
			splice @{$self->data->{$id}}, $index, 1;
		}
		if(scalar(@{$self->data->{$id}}) == 0) {
			push(@deletedIDs, $id);
			delete $self->data->{$id};
		}
	}
	return @deletedIDs;
}

sub deleteEverythingExceptIDs {
	my $self = shift;
	my @ids	= @_;
	my @deletedIDs;
	for my $id ( $self->ids ) {
		if(! any { $_ eq $id } @ids ) {
			push(@deletedIDs, $id);
			delete $self->data->{$id};
		}
	}
	return @deletedIDs;
}

sub deleteID {
	my $self = shift;
	my $id = shift;
	delete $self->data->{$id};
	return $self;
}