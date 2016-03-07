#! /usr/bin/perl
use List::Util 'max';

use strict;
use warnings;

my $prsfile = $ARGV[0];
my $freqtags = "";
my %freqs;
my @tags;
my ($i, $j) = 0;

if ($#ARGV == 1) {
	$freqtags = $ARGV[1];
} 

if ($freqtags eq "") {
	@tags = ("AFR_AF", "EUR_AF", "ASN_AF", "AMR_AF");
} else {
	@tags = split(",", $freqtags);
}

for ($i=0; $i<=$#tags; $i++) {
	$freqs{$tags[$i]}=0;
}

open (PARS, "<$prsfile");

while(<PARS>) {

	my (@freqpars);
	my $temp;
	($i, $j) = 0;

	my $line = $_;

	chomp($line);

	my @line = split("\t", $line);

	my $num = $#line;

	#print "$line[0]\t$line[2]\t$line[3]\t$line[4]\t$line[7]\t$line[$num-2]\t$line[$num-1]\t$line[$num]\n";

	my @info = split("_", $line[$num]);

	#print join("\t", @info)."\n";

	my @freq = split(";", $line[7]);

	#print join("\n", @freq)."\n\n";

	for ($i=0; $i<=$#freq; $i++) {
		for($j=0; $j<=$#tags; $j++) {
			if ($freq[$i] =~ m/$tags[$j]=([01]\.[0-9][0-9])/) {
				$freqs{$tags[$j]} = $1;
			} 
		}
	}	

	my $max = my $maxval = my $diff = 0;

	foreach my $key (keys %freqs) {
		if ($max < $freqs{$key}) {
			$max = $freqs{$key};
		}
	}

	foreach my $key (keys %freqs) {
		$diff = $max-$freqs{$key};
		if ($maxval < $diff) {
			$maxval = $diff;
		}
	}

	for ($i=0; $i<=$#tags; $i++) {
		push(@freqpars, $freqs{$tags[$i]});
	}

	push(@freqpars, $maxval);

	print "$line[0]\t$line[$num-2]\t$line[$num-1]\t".join("\t", @info)."\t$line[3]\t$line[4]\t$line[2]\t".join("\t", @freqpars)."\n";
}
