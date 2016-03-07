#! /usr/bin/perl

use strict;
use warnings;

my $resfile = $ARGV[0];
my $gwasfile = $ARGV[1];
my $gind = $ARGV[2];
my ($gwas_dist, $gwas_sidpos, $gwas_gnapos);

if ($gind == 1) {

	$gwas_dist = $ARGV[3]*1000;

} elsif ($gind == 2) {

	$gwas_sidpos = $ARGV[3]-1;

} elsif ($gind == 4) {

	$gwas_gnapos = $ARGV[3]-1;

} elsif ($gind == 3) {

	$gwas_dist = $ARGV[3]*1000;
	$gwas_sidpos = $ARGV[4]-1;

} elsif ($gind == 5) {

	$gwas_dist = $ARGV[3]*1000;
	$gwas_gnapos = $ARGV[4]-1;

} elsif ($gind == 6) {

	$gwas_sidpos = $ARGV[3]-1;
	$gwas_gnapos = $ARGV[4]-1;

} elsif ($gind == 7) {

	$gwas_dist = $ARGV[3]*1000;
	$gwas_sidpos = $ARGV[4]-1;
	$gwas_gnapos = $ARGV[5]-1;

} else {
	
	$gwas_dist = 5000;
	$gwas_sidpos = 4;
	$gwas_gnapos = 10;

}

my ($line, $chr, $pos, $dbSNP, $gwasline, $gwaschr, $gwaspos, $gwasdbSNP, $gwasname);
my (@line, @data, @gwasline, @gwasdata, @range);

open(RESDATA, "<$resfile");

while(<RESDATA>) {

	$line = $_;
	chomp($line);

	@line = split("\t", $line);

	$chr = $line[0];
	$chr = "chr$chr";
	$pos = $line[1];

	@range = ($pos-$gwas_dist,$pos+$gwas_dist);

	open(GWAS, "<$gwasfile");

	$dbSNP = "";

	while(<GWAS>) {

		$gwasline = $_;
		chomp($gwasline);

		if ($gwasline =~ m/#bin/) {
			next;
		}

		@gwasline = split("\t", $gwasline);

		$gwaschr = $gwasline[1];
		$gwaspos = $gwasline[3];
		$gwasdbSNP = $gwasline[$gwas_sidpos];
		$gwasname = $gwasline[$gwas_gnapos];

		if ( ($chr eq $gwaschr) && ($gwaspos >= $range[0]) && ($gwaspos <= $range[1]) ) {
			$dbSNP = $dbSNP."$gwasname($gwasdbSNP);";
		}
	}
	
	if ($dbSNP eq "") {
		print "$line\t.\n";
	} else {
		print "$line\t$dbSNP\n";
	}
}
