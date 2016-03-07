#! /usr/bin/perl

use strict;
use warnings;

my $trfile = $ARGV[0];    #File where all the transcripts we want to look at are stored.
my $anfile = $ARGV[1];	  #Annotation file we want to use (the script it's ready to use a GENCODE-like one.
my $chrfolder = $ARGV[2]; #Folder where all the chr that you have to check are stored. Must be in .fa (dl from UCSC).
my $BEDname = $ARGV[3];   #Name for the BED file

#Declaring some variables.
my %genetypes = ("IG_C_gene" => "IGC", "IG_D_gene" => "IGD", "IG_J_gene" => "IGJ", "IG_V_gene" => "IGV", 
		 "TR_C_gene" => "TRC", "TR_J_gene" => "TRJ", "TR_V_gene" => "TRV", "TR_D_gene" => "TRD",
		 "IG_C_pseudogene" => "IGCps", "IG_J_pseudogene" => "IGJps", "IG_V_pseudogene" => "IGVps",
		 "TR_V_pseudogene" => "TRJps", "TR_J_pseudogene" => "TRJps",
		 "Mt_rRNA" => "MtrRNA", "Mt_tRNA" => "MttRNA", "miRNA" => "miRNA", "misc_RNA" => "miscRNA", "rRNA" => "rRNA", "snRNA" => "snRNA", "snoRNA" => "snoRNA",
		 "protein_coding" => "pc", "processed_transcript" => "protra", "sense_intronic" => "senintro", "sense_overlapping" => "senover", "antisense" => "as",
		 "pseudogene" => "ps", "processed_pseudogene" => "procps", "lincRNA" => "lincRNA", " 3prime_overlapping_ncrna " => "3povernc", "polymorphic_pseudogene" => "plpd");

my %transtypes = (
"pseudogene"				=> 'ps',	"processed_transcript"			=> 'protra',	"unprocessed_pseudogene"		=> 'uprops',
"transcribed_unprocessed_pseudogene"	=> 'truprops',	"lincRNA"				=> 'lincRNA',	"miRNA"					=> 'miRNA',
"protein_coding"			=> 'pc',	"processed_pseudogene"			=> 'procps',	"antisense"				=> 'as',
"snRNA"					=> 'snRNA',	"retained_intron"			=> 'retintr',	"nonsense_mediated_decay"		=> 'nmd',
"sense_intronic"			=> 'senintro',	"misc_RNA"				=> 'miscRNA',	"transcribed_processed_pseudogene"	=> 'trprops',
"snoRNA"				=> 'snoRNA',	"non_stop_decay"			=> 'nsd',	"rRNA"					=> 'rRNA',
"unitary_pseudogene"			=> 'unitps',	"3prime_overlapping_ncrna"		=> '3povernc',	"polymorphic_pseudogene"		=> 'plpd',
"sense_overlapping"			=> 'senover',	"IG_V_gene"				=> "IGV",	"IG_C_gene"				=> "IGC",
"IG_J_gene"				=> "IGJ",	"IG_V_pseudogene"			=> "IGVps",	"TR_C_gene"				=> "TRC",
"TR_J_gene"				=> "TRJ",	"TR_V_gene"				=> "TRV",	"TR_V_pseudogene"			=> "TRVps",
"IG_C_pseudogene"			=> "IGCps",	"TR_D_gene"				=> "TRD",	"TR_J_pseudogene"			=> "TRJps",
"translated_processed_pseudogene"	=> 'traprops',	"IG_J_pseudogene"			=> "IGJps",	"IG_D_gene"				=> "IGD",
"Mt_tRNA"				=> "MttRNA",	"Mt_rRNA"				=> "MtrRNA"
);

my %pos;
my (@transcripts, @line, @data, @temp);
my $line;

open(TRANS, "<$trfile");

while(<TRANS>) {			#Getting all the transcripts id's that we have put on the file
	my $line = $_;
	chomp($line);
	push(@transcripts, $line);
}

close(TRANS);

for (my $i=0; $i<scalar(@transcripts); $i++) {	#Looping through all the transcripts.

	#Declaring some variables that we will need
	my (@genpos, $chr, $chrfile, $trid, $gt, $trt, $trseq, $temp, $strand, $chrnum, $gid);
	my (@startpos0, @startpos1, @startpos2);

	open(ANNOT, "<$anfile");

	while(<ANNOT>) { #Opening and searching the annotation file. We look for the same transcript id and the exon parts (not the CDS, start codon, stop codon, UTR, ...)

		$line = $_;

		if ($line =~ m/^##/) {
			next;
		}

		@line = split("\t", $line);
		@data = split(";", $line[8]);
		@temp = split(" ", $data[1]);

		if ($temp[1] =~ m/"(ENSTR?[0-9]*)\.[0-9]*"/) {
			$trid = $1;
		} else {
			$trid = "";
		}

		if ($trid eq $transcripts[$i]) {
			if ($line[2] eq "exon") {

				($gid) = ($line =~ m/gene_id "(ENSG[0-9]*)\.[0-9]*"/);

				@temp = split(" ", $data[2]);
				($gt) = ($line =~ m/gene_type "(\w*)"/);
				$gt = $genetypes{$gt};			#With this we get a shorter version (normally) of the gene type

				($trt) = ($line =~ m/transcript_type "(\w*)"/);
				$trt = $transtypes{$trt};		#With this we get a shorter version (normally) of the transcript type

				
				#print "$gt\n";
				$chr = $line[0];
				$strand = $line[6];
				push(@genpos, $line[3]); #Starting position
				push(@genpos, $line[4]); #Ending position

				#print $chr."\n";
				#print $strand."\n";
				#print $line[2]."\n";
				#print $line[3]."\n";
				#print $line[4]."\n\n";
			}
		}
	}

	close(ANNOT);

	$chrfile = "$chrfolder/$chr.fa"; #Getting the chromosome file path

	my $chrseq = "";

	open(CHR, "<$chrfolder/$chr.fa");

	while(<CHR>) {

		$line = $_;	
		chomp($line);

		if ($line =~ m/^>chr[0-9]*.*/) {
			next;
		}

		$chrseq = $chrseq.$line; #Getting the sequence of the chromosome in a single line, so we can search with the genomic starting and ending positions that we get earlier.
	}

	close(CHR);

	my $count = 1; #Counter for the nucleotide of the transcript
	
	#my $totexons = scalar(@genpos)/2; #2 variables to count the number of exons (needed to look at the info of the prints below).
	#my $texns = scalar(@genpos) - 2;
		
	for (my $j=0; $j<$#genpos; $j+=2) {

		my $start = $genpos[$j];#Getting the starting and ending positions from the array
		my $end = $genpos[$j+1];
		my $trlen = $end-$start+1;
			
		#my $numexo = $totexons - ($texns - $j)/2;	#Some info (perhaps we can put this as an option).
		#print "Exon $numexo of $transcripts[$i]\n";
		#if ($strand eq "+") {
		#	print "Starts at: $start\n";
		#	print "Ends at: $end\n";
		#} elsif ($strand eq "-") {
		#	print "Starts at: $end\n";
		#	print "Ends at: $start\n";
		#}
		#	print "Trans length: $trlen\n\n";


		$temp = substr($chrseq, $start-1, $trlen);	#Getting the transcript sequence from the chromosome sequence.

		if ($strand eq "+") {

			#print "$temp\n";
			$trseq = $trseq.$temp;	#Concatenating the exon sequence to the rest of the sequence
					
			for (my $k=0; $k<=$trlen-1; $k++) {	#Counting the sequence. keys are sequence positions and values are genomic positions.
				$pos{$count} = $start+$k;
				#print "$count: $pos{$count}\n";
				$count++;
			}

		} elsif ($strand eq "-") {

			$temp = reverse($temp);		#Getting the complementary chain through reverse and translate.
			$temp =~ tr/ACGTacgt/TGCAtgca/;
			#print "$temp\n";

			$trseq = $trseq.$temp;
		
			for (my $k=$trlen-1; $k>=0; $k--) {	#The same as above, but we have to go backwards.
				$pos{$count} = $start+$k;
				#print "$count: $pos{$count}\n";
				$count++;
			}
		}
		#print "\n";
			
	}

	$trseq = uc($trseq);

	for(my $j=1; $j<=length($trseq); $j+=3) { #Counting start codons (ATG) in the 3 different frames (+0, +1, +2), and keeping the position of the A to find things later.
		for (my $k=0; $k<3; $k++) {

			if ((length($trseq)-$j-$k)<3) { #To avoid the end of the sequence.
				next;
			}

			my $cod = substr($trseq, $j+$k-1, 3);

			if ($cod eq "ATG") {
				if ($k==0) {
					push(@startpos0, $j+$k);
				} elsif ($k==1) {
					push(@startpos1, $j+$k);
				} elsif ($k==2) {
					push(@startpos2, $j+$k);
				}
			}
		}
	}

	($chrnum) = ($chr =~ m/chr([0-9XY]+)/);

	open(BED, ">>$BEDname.bed"); #open the BED file and append the data as needed.

	if (@startpos0) { #Starting with the first frame (+0)

		my %stops; #A hash to store the stop positions.

		for (my $j=0; $j<=$#startpos0; $j++) { #Looping through all the starting codons.

			my $array = [];

			my $subtrseq = substr($trseq, $startpos0[$j]-1); #Getting the sub-transcript (from the start codon to the end) to count stop codons

			for (my $k=1; $k<=length($subtrseq); $k+=3) {

				my $cod = substr($subtrseq, $k-1, 3);

				if ($cod eq "TAA"){
						push(@$array, $k+$startpos0[$j]-1); #Storing the position of th T.
				} elsif ($cod eq "TAG"){
						push(@$array, $k+$startpos0[$j]-1);
				} elsif ($cod eq "TGA"){
						push(@$array, $k+$startpos0[$j]-1);
				}
			}
		
			$stops{($j+1)} = $array;
		}

		my $startcount = 0; #More counters...

		foreach my $key (sort {$a <=> $b} keys %stops) {

			my $tempos = $startpos0[$startcount]; #To get an easier name.
			my $stopcount = 1;
			$startcount++;

			if ($strand eq "+") {

				print "$transcripts[$i]\t$gt\t$trt\t$strand\tATG0.$startcount\tATG\t$tempos,".($tempos+1).",".($tempos+2).
					"\t-\t$chr\_".($pos{$tempos})."_".($pos{$tempos+2})."_ATG\n"; #Parsing of the start codons in a transcript on the + strand. Below an example.
				
				#ENST00000473358.1	+	ATG0.1	ATG	16,17,18	-	chr1_29568_29570_ATG

				#Start codons are codified as follows: ATG.(frame{0,1,2}.(number)

			} elsif ($strand eq "-") {

				print "$transcripts[$i]\t$gt\t$trt\t$strand\tATG0.$startcount\tATG\t$tempos,".($tempos+1).",".($tempos+2).
					"\t-\t$chr\_".($pos{$tempos+2})."_".($pos{$tempos})."_CAT\n";#Parsing of the start codons in a transcript on the - strand. Below an example.

				#ENST00000441386.2	-	ATG0.1	ATG	94,95,96	-	chr3_5021526_5021528_CAT

				#Start codons are codified as follows: ATG.(frame{0,1,2}.(number)

			}

			foreach (@{$stops{$key}}) {

				if ($stopcount > 2 ) {
					last;
				}

				my $stopos = $_;
	       				
				my $cod = substr($trseq, $stopos-1, 3);

				my $len = $stopos+2-$tempos+1;

				if ($strand eq "+") {

					print "$transcripts[$i]\t$gt\t$trt\t$strand\tStop0.$startcount.$stopcount\t$cod\t$stopos,".($stopos+1).",".($stopos+2).
						"\t$len\t$chr\_".($pos{$stopos})."_".($pos{$stopos+2})."_$cod\n"; #Parsing for the stop codons in a transcript on the + strand. Below an ex.
					
					#ENST00000473358.1	+	Stop0.1.1	TGA	241,242,243	228	chr1_29793_29795_TGA

					#Stop codons are codified as follows: TNN.(frame{0,1,2}.(number).(number).
					#The frame and the first number is the same as the last start codon, so they are realted. ATG0.1 can have stop codons like TNN.0.1.X

					if ($stopcount == 1) { #Writting to a BED file so it can be intersected. Note that if it doesn't have a second STOP codon the final part (the length
								#of the transcript taking into account the second stop codon) will be 0.
						my $len2 = 0;						
						#my $len2 = length($trseq);
			
						if (defined(${$stops{$key}}[1]))  {
							my $stopos2 = ${$stops{$key}}[1];
							$len2 = $stopos2+2-$tempos+1;
						}
						
						my($p1, $p2, $p3) = split('',$cod);

				print BED "$chrnum\t".($pos{$stopos}-1).  "\t".($pos{$stopos}).  "\t$transcripts[$i]\_$gid\_$gt\_$trt\_$strand\_Stop0.$startcount.$stopcount.1\_$cod\_$p1\_$len\_$len2\n";
				print BED "$chrnum\t".($pos{$stopos+1}-1)."\t".($pos{$stopos+1})."\t$transcripts[$i]\_$gid\_$gt\_$trt\_$strand\_Stop0.$startcount.$stopcount.2\_$cod\_$p2\_$len\_$len2\n";
				print BED "$chrnum\t".($pos{$stopos+2}-1)."\t".($pos{$stopos+2})."\t$transcripts[$i]\_$gid\_$gt\_$trt\_$strand\_Stop0.$startcount.$stopcount.3\_$cod\_$p3\_$len\_$len2\n";
						#Parsing of the data to fit a BED file. Below an example.

						#1	12144	12145	ENST00000456328.2_ps_+_Stop0.1.1.1_TAA_T_57_285
						#1	12145	12146	ENST00000456328.2_ps_+_Stop0.1.1.2_TAA_A_57_285
						#1	12146	12147	ENST00000456328.2_ps_+_Stop0.1.1.3_TAA_A_57_285
		
					}
				
				} elsif ($strand eq "-") {

					my $rcod = reverse($cod);
					$rcod =~ tr/ACGTacgt/TGCAtgca/;

					print "$transcripts[$i]\t$gt\t$trt\t$strand\tStop0.$startcount.$stopcount\t$cod\t$stopos,".($stopos+1).",".($stopos+2).
						"\t$len\t$chr\_".($pos{$stopos+2})."_".($pos{$stopos})."_$rcod\n";

					if ($stopcount == 1) {
	
						my $len2 = 0;						
						#my $len2 = length($trseq);
			
						if (defined(${$stops{$key}}[1]))  {
							my $stopos2 = ${$stops{$key}}[1];
							$len2 = $stopos2+2-$tempos+1;
						}
	
						my($p1, $p2, $p3) = split('',$rcod);

				print BED "$chrnum\t".($pos{$stopos+2}-1)."\t".($pos{$stopos+2})."\t$transcripts[$i]\_$gid\_$gt\_$trt\_$strand\_Stop0.$startcount.$stopcount.1\_$cod\_$p1\_$len\_$len2\n";
				print BED "$chrnum\t".($pos{$stopos+1}-1)."\t".($pos{$stopos+1})."\t$transcripts[$i]\_$gid\_$gt\_$trt\_$strand\_Stop0.$startcount.$stopcount.2\_$cod\_$p2\_$len\_$len2\n";
				print BED "$chrnum\t".($pos{$stopos}-1).  "\t".($pos{$stopos}).  "\t$transcripts[$i]\_$gid\_$gt\_$trt\_$strand\_Stop0.$startcount.$stopcount.3\_$cod\_$p3\_$len\_$len2\n";
						#Parsing of the data to fit a BED file. Below an example.
	
						#1	14453	14454	ENST00000438504.2_ps_-_Stop0.1.1.1_TAA_T_234_0
						#1	14454	14455	ENST00000438504.2_ps_-_Stop0.1.1.2_TAA_T_234_0
						#1	14455	14456	ENST00000438504.2_ps_-_Stop0.1.1.3_TAA_A_234_0
		
					}
				
				}
				$stopcount++;
			}
		}

	} else { print "$transcripts[$i]\t$gid\t$gt\t$trt\t$strand\tdon't have a start codon on ORF0\n"; }

	if (@startpos1) { #Same explanations as @startpos0

		my %stops;

		for (my $j=0; $j<=$#startpos1; $j++) {

			my $array = [];

			my $subtrseq = substr($trseq, $startpos1[$j]-1);

			for (my $k=1; $k<=length($subtrseq); $k+=3) {

				my $cod = substr($subtrseq, $k-1, 3);

				if ($cod eq "TAA"){
						push(@$array, $k+$startpos1[$j]-1);
				} elsif ($cod eq "TAG"){
						push(@$array, $k+$startpos1[$j]-1);
				} elsif ($cod eq "TGA"){
						push(@$array, $k+$startpos1[$j]-1);
				}
			}
		
			$stops{($j+1)} = $array;
		}

		my $startcount = 0;

		foreach my $key (sort {$a <=> $b} keys %stops) {

			my $tempos = $startpos1[$startcount];
			my $stopcount = 1;
			$startcount++;

			if ($strand eq "+") {

				print "$transcripts[$i]\t$gid\t$gt\t$trt\t$strand\tATG1.$startcount\tATG\t$tempos,".($tempos+1).",".($tempos+2).
					"\t-\t$chr\_".($pos{$tempos})."_".($pos{$tempos+2})."_ATG\n";

			} elsif ($strand eq "-") {

				print "$transcripts[$i]\t$gid\t$gt\t$trt\t$strand\tATG1.$startcount\tATG\t$tempos,".($tempos+1).",".($tempos+2).
					"\t-\t$chr\_".($pos{$tempos+2})."_".($pos{$tempos})."_CAT\n";

			}
		
			foreach (@{$stops{$key}}) {

				if ($stopcount > 2 ) {
					last;
				}

				my $stopos = $_;
	       				
				my $cod = substr($trseq, $stopos-1, 3);

				my $len = $stopos+2-$tempos+1;

				if ($strand eq "+") {

					print "$transcripts[$i]\t$gid\t$gt\t$trt\t$strand\tStop1.$startcount.$stopcount\t$cod\t$stopos,".($stopos+1).",".($stopos+2).
						"\t$len\t$chr\_".($pos{$stopos})."_".($pos{$stopos+2})."_$cod\n";

					if ($stopcount == 1) {

						my $len2 = 0;						
						#my $len2 = length($trseq);
			
						if (defined(${$stops{$key}}[1]))  {
							my $stopos2 = ${$stops{$key}}[1];
							$len2 = $stopos2+2-$tempos+1;
						}
						
						my($p1, $p2, $p3) = split('',$cod);

				print BED "$chrnum\t".($pos{$stopos}-1).  "\t".($pos{$stopos}).  "\t$transcripts[$i]\_$gid\_$gt\_$trt\_$strand\_Stop1.$startcount.$stopcount.1\_$cod\_$p1\_$len\_$len2\n";
				print BED "$chrnum\t".($pos{$stopos+1}-1)."\t".($pos{$stopos+1})."\t$transcripts[$i]\_$gid\_$gt\_$trt\_$strand\_Stop1.$startcount.$stopcount.2\_$cod\_$p2\_$len\_$len2\n";
				print BED "$chrnum\t".($pos{$stopos+2}-1)."\t".($pos{$stopos+2})."\t$transcripts[$i]\_$gid\_$gt\_$trt\_$strand\_Stop1.$startcount.$stopcount.3\_$cod\_$p3\_$len\_$len2\n";
		
					}
				
				} elsif ($strand eq "-") {

					my $rcod = reverse($cod);
					$rcod =~ tr/ACGTacgt/TGCAtgca/;

					print "$transcripts[$i]\t$gid\t$gt\t$trt\t$strand\tStop1.$startcount.$stopcount\t$cod\t$stopos,".($stopos+1).",".($stopos+2).
						"\t$len\t$chr\_".($pos{$stopos+2})."_".($pos{$stopos})."_$rcod\n";

					if ($stopcount == 1) {
	
						my $len2 = 0;						
						#my $len2 = length($trseq);
			
						if (defined(${$stops{$key}}[1]))  {
							my $stopos2 = ${$stops{$key}}[1];
							$len2 = $stopos2+2-$tempos+1;
						}

						my($p1, $p2, $p3) = split('',$rcod);

				print BED "$chrnum\t".($pos{$stopos+2}-1)."\t".($pos{$stopos+2})."\t$transcripts[$i]\_$gid\_$gt\_$trt\_$strand\_Stop1.$startcount.$stopcount.1\_$cod\_$p1\_$len\_$len2\n";
				print BED "$chrnum\t".($pos{$stopos+1}-1)."\t".($pos{$stopos+1})."\t$transcripts[$i]\_$gid\_$gt\_$trt\_$strand\_Stop1.$startcount.$stopcount.2\_$cod\_$p2\_$len\_$len2\n";
				print BED "$chrnum\t".($pos{$stopos}-1).  "\t".($pos{$stopos}).  "\t$transcripts[$i]\_$gid\_$gt\_$trt\_$strand\_Stop1.$startcount.$stopcount.3\_$cod\_$p3\_$len\_$len2\n";
		
					}
				
				}
				$stopcount++;
			}
		}

	} else { print "$transcripts[$i]\t$gid\t$gt\t$trt\t$strand\tdon't have a start codon on ORF1\n"; }

	if (@startpos2) { #Same explanations as @startpos0

		my %stops;

		for (my $j=0; $j<=$#startpos2; $j++) {

			my $array = [];

			my $subtrseq = substr($trseq, $startpos2[$j]-1);

			for (my $k=1; $k<=length($subtrseq); $k+=3) {

				my $cod = substr($subtrseq, $k-1, 3);

				if ($cod eq "TAA"){
						push(@$array, $k+$startpos2[$j]-1);
				} elsif ($cod eq "TAG"){
						push(@$array, $k+$startpos2[$j]-1);
				} elsif ($cod eq "TGA"){
						push(@$array, $k+$startpos2[$j]-1);
				}
			}
			$stops{($j+1)} = $array;
		}

		my $startcount = 0;

		foreach my $key (sort {$a <=> $b} keys %stops) {

			my $tempos = $startpos2[$startcount];
			my $stopcount = 1;
			$startcount++;

			if ($strand eq "+") {

				print "$transcripts[$i]\t$gid\t$gt\t$trt\t$strand\tATG2.$startcount\tATG\t$tempos,".($tempos+1).",".($tempos+2).
					"\t-\t$chr\_".($pos{$tempos})."_".($pos{$tempos+2})."_ATG\n";

			} elsif ($strand eq "-") {

				print "$transcripts[$i]\t$gid\t$gt\t$trt\t$strand\tATG2.$startcount\tATG\t$tempos,".($tempos+1).",".($tempos+2).
					"\t-\t$chr\_".($pos{$tempos+2})."_".($pos{$tempos})."_CAT\n";

			}

		

			foreach (@{$stops{$key}}) {

				if ($stopcount > 2 ) {
					last;
				}

				my $stopos = $_;
				
	       				
				my $cod = substr($trseq, $stopos-1, 3);

				my $len = $stopos+2-$tempos+1;

				if ($strand eq "+") {

					print "$transcripts[$i]\t$gid\t$gt\t$trt\t$strand\tStop2.$startcount.$stopcount\t$cod\t$stopos,".($stopos+1).",".($stopos+2).
						"\t$len\t$chr\_".($pos{$stopos})."_".($pos{$stopos+2})."_$cod\n";

					if ($stopcount == 1) {

						my $len2 = 0;						
						#my $len2 = length($trseq);
			
						if (defined(${$stops{$key}}[1]))  {
							my $stopos2 = ${$stops{$key}}[1];
							$len2 = $stopos2+2-$tempos+1;
						}

						my($p1, $p2, $p3) = split('',$cod); # "Little" error on the bed file!!!!!

				print BED "$chrnum\t".($pos{$stopos}-1).  "\t".($pos{$stopos}).  "\t$transcripts[$i]\_$gid\_$gt\_$trt\_$strand\_Stop2.$startcount.$stopcount.1\_$cod\_$p1\_$len\_$len2\n";
				print BED "$chrnum\t".($pos{$stopos+1}-1)."\t".($pos{$stopos+1})."\t$transcripts[$i]\_$gid\_$gt\_$trt\_$strand\_Stop2.$startcount.$stopcount.2\_$cod\_$p2\_$len\_$len2\n";
				print BED "$chrnum\t".($pos{$stopos+2}-1)."\t".($pos{$stopos+2})."\t$transcripts[$i]\_$gid\_$gt\_$trt\_$strand\_Stop2.$startcount.$stopcount.3\_$cod\_$p3\_$len\_$len2\n";
		
					}
				
				} elsif ($strand eq "-") {

					my $rcod = reverse($cod);
					$rcod =~ tr/ACGTacgt/TGCAtgca/;

					print "$transcripts[$i]\t$gid\t$gt\t$trt\t$strand\tStop2.$startcount.$stopcount\t$cod\t$stopos,".($stopos+1).",".($stopos+2).
						"\t$len\t$chr\_".($pos{$stopos+2})."_".($pos{$stopos})."_$rcod\n";

					if ($stopcount == 1) {
	
						my $len2 = 0;						
						#my $len2 = length($trseq);
			
						if (defined(${$stops{$key}}[1]))  {
							my $stopos2 = ${$stops{$key}}[1];
							$len2 = $stopos2+2-$tempos+1;
						}

						my($p1, $p2, $p3) = split('',$rcod); # "Little" error on the bed file!!!!!

				print BED "$chrnum\t".($pos{$stopos+2}-1)."\t".($pos{$stopos+2})."\t$transcripts[$i]\_$gid\_$gt\_$trt\_$strand\_Stop2.$startcount.$stopcount.1\_$cod\_$p1\_$len\_$len2\n";
				print BED "$chrnum\t".($pos{$stopos+1}-1)."\t".($pos{$stopos+1})."\t$transcripts[$i]\_$gid\_$gt\_$trt\_$strand\_Stop2.$startcount.$stopcount.2\_$cod\_$p2\_$len\_$len2\n";
				print BED "$chrnum\t".($pos{$stopos}-1).  "\t".($pos{$stopos}).  "\t$transcripts[$i]\_$gid\_$gt\_$trt\_$strand\_Stop2.$startcount.$stopcount.3\_$cod\_$p3\_$len\_$len2\n";
		
					}
				
				}
				$stopcount++;
			}
		}
	} else { print "$transcripts[$i]\t$gid\t$gt\t$trt\t$strand\tdon't have a start codon on ORF2\n"; }

	#print "\n";
	#print BED "\n";

	close(BED);

}


