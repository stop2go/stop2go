stop2go is a program designed to check for polymorphisms on stop codons that change them to a coding state.

This program needs the following arguments to run:
	--input     || -i:   input file for the program to work. List of transcripts for the program to work. Ony ENSEMBL transcript IDs in the form ENST00000426466 are admitted.
	--annot     || -an:  annotation file. The annotation format admitted is GTF format.
	--folder    || -fo:  folder where FASTA files with the sequence of the chromosomes are stored.
	--freq      || -fr:  VCF file containing the frequences (only accepts 1 file)
	
stop2go also have the following optional commands:
	--output    || -o  : name for the output file. If this command is not used the name will be the same as the input file.
	--freqtag          : input your custom tags for the alternative frequencies on your VCF file (separate them by commas).
	--ancestral || -anc: file to search for the ancestral allele (only accepts 1 file).
	--anctag           : input your custom tags for the ancestral allele (separate them by commas).

By default stop2go assumes that the frequency tags are the following: AFR_AF, EUR_AF, ASN_AF, AMR_AF. They will be displayed in this order. If you use your own tags they will be displayed in the same order you input them.

By default stop2go assumes that the ancestral tags are the following: CAnc, GAnc, OAnc. They will be searched in this order, but only the first matching will be displayes in the results file.

You can add new information to your results file using the argument --add with your results file (e.g.: --add res_example.txt). With this you can avoid re-running the program if you have a large dataset.

If you want to run an example, please input this command in your terminal:

bash stop2go.sh --input ./example/example_list.txt --output example --annot ./example/example_annot.gtf --folder ./example/ --freq ./example/example.vcf

The program will output 3 different files:

output_list.txt			This file contains the collection of all ORFs found in the transcripts taking into account the 3 reading frames. Please be aware that a stop codon may have more than 1 associated start.
output.bed			BED file with all the stop codon positions to intersect with VCF files.
res_1kg_output.txt			This file contains information about those stops that change from a stop state to a coding state.

The output results file (res_output.txt) will have the following fields (separeted by tabulators):

Chromosome			19Y	
Stop position -1		1277	
Stop position			1278	
Trasncript ID			ENST00000426466	
Gene ID				ENSG00000163098
Gene type 			pc	
Transcript type 		pc	
Strand 				-	
Stop ID				Stop2.10.1.2	
Stop codon affected		TGA	
Base affected			C	
Original length			9	
Extended length			54	
Affected allele			A	
Alternative allele		G	
SNP id				rs34092035	
AFR frequency			0.13	
EUR frequency			0	
ASN frequency			0	
AMR frequency			0.02	
Maximum frequency difference	0.13


Other fields will be added, depending on the arguments used to run stop2go.


###================================================================================================================###

This are the correspondences between the gene type abrevations and the gene types found in the annotation:
IG_C_gene 			=> IGC 
IG_D_gene 			=> IGD 
IG_J_gene 			=> IGJ 
IG_V_gene 			=> IGV 		 
TR_C_gene 			=> TRC 
TR_J_gene 			=> TRJ 
TR_V_gene			=> TRV 
TR_D_gene 			=> TRD 
IG_C_pseudogene 		=> IGCps 
IG_J_pseudogene 		=> IGJps 
IG_V_pseudogene 		=> IGVps
TR_V_pseudogene 		=> TRJps 
TR_J_pseudogene 		=> TRJps
Mt_rRNA 			=> MtrRNA 
Mt_tRNA 			=> MttRNA 
miRNA				=> miRNA 
misc_RNA 			=> miscRNA 
rRNA 				=> rRNA 
snRNA 				=> snRNA 
snoRNA 				=> snoRNA
protein_coding 			=> pc 
processed_transcript 		=> protra 
sense_intronic 			=> senintro 
sense_overlapping 		=> senover 
antisense 			=> as
pseudogene 			=> ps 
processed_pseudogene 		=> procps
lincRNA 			=> lincRNA 
3prime_overlapping_ncrna  	=> 3povernc 
polymorphic_pseudogene 		=> plpd

This are the correspondences between the transcript type abrevation and the transcript types found on the annotation:
pseudogene				=> ps
processed_transcript			=> protra
unprocessed_pseudogene			=> uprops
transcribed_unprocessed_pseudogene	=> truprops
lincRNA					=> lincRNA
miRNA					=> miRNA
protein_coding				=> pc
processed_pseudogene			=> procps
antisense				=> as
snRNA					=> snRNA
retained_intron				=> retintr
nonsense_mediated_decay			=> nmd
sense_intronic				=> senintro
misc_RNA				=> miscRNA
transcribed_processed_pseudogene	=> trprops
snoRNA					=> snoRNA
non_stop_decay				=> nsd
rRNA					=> rRNA
unitary_pseudogene			=> unitps
3prime_overlapping_ncrna		=> 3povernc
polymorphic_pseudogene			=> plpd
sense_overlapping			=> senover
IG_V_gene				=> IGV
IG_C_gene				=> IGC
IG_J_gene				=> IGJ
IG_V_pseudogene				=> IGVps
TR_C_gene				=> TRC
TR_J_gene				=> TRJ
TR_V_gene				=> TRV
TR_V_pseudogene				=> TRVps
IG_C_pseudogene				=> IGCps
TR_D_gene				=> TRD
TR_J_pseudogene				=> TRJps
translated_processed_pseudogene		=> traprops
IG_J_pseudogene				=> IGJps
IG_D_gene				=> IGD
Mt_tRNA					=> MttRNA
Mt_rRNA					=> MtrRNA
