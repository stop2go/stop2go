#! /usr/bin/perl

use strict;
use warnings;

my $prsfile = $ARGV[0];
my $del = quotemeta(".");

open (PARS, "<$prsfile");

while(<PARS>) {

	my $line = $_;
	chomp($line);
	my @line = split("\t", $line);
	my @temp = split($del, $line[8]);

	my $cod = $line[9];
	my $pos = $temp[$#temp];
	my $strand = $line[7];
	my $alt = $line[14];

	if ($strand eq "+") {

		substr($cod, $pos-1, 1) = $alt;

		if ($cod eq "TAA") {
			next;
		} elsif ($cod eq "TAG") {
			next;
		} elsif ($cod eq "TGA") {
			next;
		} else {
			print "$line\n";
		}

	} elsif ($strand eq "-") {
		
		$cod = reverse($cod);
		$cod =~ tr/ACGT/TGCA/;
		substr($cod, $pos-1, 1) = $alt;
		$cod = reverse($cod);
		$cod =~ tr/ACGT/TGCA/;

		if ($cod eq "TAA") {
			next;
		} elsif ($cod eq "TAG") {
			next;
		} elsif ($cod eq "TGA") {
			next;
		} else {
			print "$line\n";
		}
	}
}
