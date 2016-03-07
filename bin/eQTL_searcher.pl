#! /usr/bin/perl

use strict;
use warnings;

my $resfile = $ARGV[0];
my $qtlfile = $ARGV[1];
my $qind = $ARGV[2];
my ($qtl_sidpos,$qtl_gidpos,$qtl_qtlpos) = 0;

if ($qind == 1) {

	$qtl_sidpos = $ARGV[3]-1;

} elsif ($qind == 2) {

	$qtl_gidpos = $ARGV[3]-1;

} elsif ($qind == 4) {

	$qtl_qtlpos = $ARGV[3]-1;

} elsif ($qind == 3) {

	$qtl_sidpos = $ARGV[3]-1;
	$qtl_gidpos = $ARGV[4]-1;

} elsif ($qind == 5) {

	$qtl_sidpos = $ARGV[3]-1;
	$qtl_qtlpos = $ARGV[4]-1;

} elsif ($qind == 6) {

	$qtl_gidpos = $ARGV[3]-1;
	$qtl_qtlpos = $ARGV[4]-1;

} elsif ($qind == 7) {

	$qtl_sidpos = $ARGV[3]-1;
	$qtl_gidpos = $ARGV[4]-1;
	$qtl_qtlpos = $ARGV[5]-1;

} else {
	
	$qtl_sidpos = 0;
	$qtl_gidpos = 2;
	$qtl_qtlpos = 9;

}

my ($line, $qtlline, $res_snpid, $res_gid, $qtl_snpid, $qtl_gid, $qtl);
my (@line, @qtlline);

open(DATA, "<$resfile");

while(<DATA>) {

	$line = $_;
	chomp($line);
	
	@line = split("\t", $line);

	$res_snpid = $line[15];
	$res_gid = $line[4];
	$qtl = ".";

	open(QTL, "<$qtlfile");

	while(<QTL>) {

		$qtlline = $_;
		chomp($qtlline);

		if ($qtlline =~ m/SNP_ID/) {
			next;
		}

		@qtlline = split("\t", $qtlline);

		$qtl_snpid = $qtlline[$qtl_sidpos];
		($qtl_gid) = ($qtlline[$qtl_gidpos] =~ m/(ENSGR?\d+)\.\d+/);

		if ($qtl_snpid eq $res_snpid) {
			if ($qtl_gid eq $res_gid) {
				$qtl = $qtlline[$qtl_qtlpos];
			}
		}
	}

	close(QTL);

	print join("\t", @line)."\t$qtl\n";

}
