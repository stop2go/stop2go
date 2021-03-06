#! /usr/bin/perl
use List::Util 'max';

use strict;
use warnings;

my $parsfile = $ARGV[0];
my $anctags = "";
my @tags;

if ($#ARGV == 1) {
	$anctags = $ARGV[1];
} 

if ($anctags eq "") {
	@tags = ("CAnc", "GAnc", "OAnc");
}

my ($line, $ancal, $preline, $ref, $alt, $temp, $fafr, $fceu, $fasn, $famr, $nfafr, $nfceu, $nfasn, $nfamr, $max, $num);
my (@line, @data);

open(DATA, "<$parsfile");

while(<DATA>) {

	$line = $_;
	chomp($line);

	@line = split("\t", $line);

	$num = $#line;

	$ancal = ".";

	$preline = "";
	
	for (my $i=$num-20; $i<=$num-8; $i++) {
		$preline = $preline.$line[$i]."\t";
	}

		
	if ($anctags eq "") {
		for (my $i=0; $i<=$#tags; $i++) {
			if ($line =~ m/$tags[$i]=([ACGT]+)/) {
				$ancal = $1;
			}
		}
	} else {
		if ($line =~ m/$anctags=([ACGT]+)/) {
			$ancal = $1;
		}
	}

	if (length($ancal) != 1) {
		$ancal=".";
	}

	$ref = $line[$num-7];
	$alt = $line[$num-6];

	$fafr = $line[$num-4];
	$fceu = $line[$num-3];
	$fasn = $line[$num-2];
	$famr = $line[$num-1];

	if ($ancal ne ".") {

		if ($ref ne $ancal) {

			$temp = $ref;
			$ref = $alt;
			$alt = $temp;

			$nfafr = 1-$fafr;
			$nfceu = 1-$fceu;
			$nfasn = 1-$fasn;
			$nfamr = 1-$famr;

			$max = max( ($nfafr, $nfceu, $nfasn, $nfamr ) );

			$max = max( ($max-$nfafr, $max-$nfceu, $max-$nfasn, $max-$nfamr) );

			print "$preline$ancal\t$ref*\t$alt*\t$line[$num-5]\t$nfafr\t$nfceu\t$nfasn\t$nfamr\t$max\n";
			
		} else {
			print "$preline$ancal\t$ref\t$alt\t$line[$num-5]\t$fafr\t$fceu\t$fasn\t$famr\t$line[$num]\n";
		}

	} else {
		print "$preline$ancal\t$ref\t$alt\t$line[$num-5]\t$fafr\t$fceu\t$fasn\t$famr\t$line[$num]\n";;
	}
}
