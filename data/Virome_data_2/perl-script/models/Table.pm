package Table;

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
			@fields = split $headerRegex, $_;
			s{^\s+|\s+$}{}g foreach @fields; # Remove leading whitespaces
			if(!exists $self->{'header'}) {		
				@{$self->{'header'}} = @fields;
			} else {								
				for my $index (0 .. $#fields) {
					my $fieldName = @{$self->{'header'}}[$index];
					$self->{'data'}{$fields[0]}{$fieldName} = $fields[$index];
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

sub data {
	my $self = shift;
	return $self->{'data'};
}

sub rows {
	my $self = shift;
	return values %{$self->data};
}

sub sumForField {
	my $self = shift;
	my $field = shift;
	my $sum = 0.0;
	for my $id ( $self->ids ) {
		$sum += $self->data->{$id}{$field};
	}
	return $sum;
}

sub filterByField {
	my $self = shift;
	my $field = shift;
	my $text = shift;
	my @deletedIDs;
	for my $id ($self->ids) {
		if(index(lc $self->data->{$id}{$field}, lc $text) == -1) {
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
		if(eval($self->data->{$id}{$field}." ".$expression)) {
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