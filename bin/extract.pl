#! /usr/bin/perl

use strict;
use warnings;

my $annotfile = $ARGV[0];
my $aind = $ARGV[1];
my ($tag, $type);

if ($aind == 1) {

	$tag = $ARGV[2];
	$type = "ENSTR?";

} elsif ($aind == 2) {
	
	$tag = "transcript_id";	
	$type = $ARGV[2];

} elsif ($aind == 3) {

	$tag = $ARGV[2];
	$type = $ARGV[3];

} else {

	$tag = "transcript_id";
	$type = "ENSTR?";

}

my ($line, $trid);
my (@line, @data, @temp);

open(ANNOT, "<$annotfile");

while(<ANNOT>) { #Opening and searching the annotation file. We look for the same transcript id and the exon parts (not the CDS, start codon, stop codon, UTR, ...)

	$line = $_;

	if ($line =~ m/^##/) {
		next;
	}

	@line = split("\t", $line);
	@data = split(";", $line[8]);

	#print STDERR join("\n", @data)."\n";
		
	for (my $i = 0; $i<scalar(@data)-1; $i++) {

		@temp = split(" ", $data[$i]);

		#print STDERR "$temp[0]\n";
		if ($temp[0] =~ m/$tag/) {
			($trid) = ($temp[1] =~ m/"($type[0-9]*)\.[0-9]*"/);
		}
	}

	if ($trid) {
		print "$trid\n";
	}
}
