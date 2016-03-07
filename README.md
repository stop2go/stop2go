# stop2go: the way to go!

## Information

stop2go is a program designed to check for polymorphisms on stop codons that change them to a coding state. It is thought to work in a very flexible manner, with most of its options controlled by the user. It accepts generally used file formats:

- fasta format for genetic sequences
- vcf format for frequencies on the SNPs
- gtf format for annotation files

## Usage

To run stop2go you have to use stop2go.sh with the following arguments to run:

* --input (-i):   input file for the program to work. List of transcripts for the program to work. Ony ENSEMBL transcript IDs in the form ENST00000426466 are admitted.

* --annot (-an):  annotation file. The annotation format admitted is GTF format.

* --folder (-fo):  folder where FASTA files with the sequence of the chromosomes are stored.

* --freq (-fr):  VCF file containing the frequences (only accepts 1 file).By default stop2go assumes that the frequency tags are the following: AFR_AF, EUR_AF, ASN_AF, AMR_AF. They will be displayed in this order. If you use your own tags they will be displayed in the same order you input them.

Please, be aware that the fasta file of the chromosomes has to match the chromosome name on the gtf file (typically in the first field). For example, if the annotation file has "chr1" the chromosome fasta file has to be "chr1.fa". Otherwise the program will not work
	
stop2go also have the following optional commands:

* --output (-o): name for the output file. If this command is not used the name will be the same as the input file.
* --freqtag: input your custom tags for the alternative frequencies on your VCF file (separate them by commas).

* --ancestral (-anc): file to search for the ancestral allele (only accepts 1 file).By default stop2go assumes that the ancestral tags are the following: CAnc, GAnc, OAnc. They will be searched in this order, but only the first match will be displayed in the results file.
* --anctag: input your custom tags for the ancestral allele (separate them by commas).


* --qtl: name for the file where qtl information is stored.
* --qtl_spos: optional argument that indicates the position of the SNP ID inside the qtl file. By default stop2go assumes it will be in the first field of the file.
* --qtl_gpos: optional argument that indicates the position of the gene ID is the qtl affecting inside the qtl file. By default stop2go assumes it will be in the third field of the file.
* --qtl_qpos: ooptional argument that indicates the position of the strength and direction of the qtl inside the qtl file. By default stop2go assumes it will be in the tenth field of the file.


* --gwas: name for the file where GWAS information is stored.
* --gwas_dist: optional argument that indicates the distance (in kilobases) at which GWAS hits may be searched. By default stop2go assumes it will be 5.
* --gwas_spos: optional argument that indicates the position of the SNP ID inside the qtl file. By default stop2go assumes it will be in the first field of the file.
* --gwas_npos: optional argument that indicates the position of the SNP ID inside the qtl file. By default stop2go assumes it will be in the first field of the file.

You can add new information to your results file using the argument --add with your results file (e.g.: --add res_example.txt). With this you can avoid re-running the program if you have a large dataset and you wish to add some extra info.


## Example

If you want to run an example, please input this command in your terminal:

bash stop2go.sh --input ./example/example_list.txt --output example --annot ./example/example_annot.gtf --folder ./example/ --freq ./example/example.vcf

All the example files can be found in this project under the "example" folder if you wish to check them.

The program will output 3 different files:

output_list.txt: This file contains the collection of all ORFs found in the transcripts taking into account the 3 reading frames. Please be aware that a stop codon may have more than 1 associated start.

output.bed: BED file with all the stop codon positions to intersect with VCF files.

res_1kg_output.txt: This file contains information about those stops that change from a stop state to a coding state.

The output results file (res_output.txt) will have the following fields, separeted by tabulators (information between brackets corresponds to the first line of the example output):

1. Chromosome (19Y)
2. Stop position -1 (1277)
3. Stop position (1278)
4. Trasncript ID (ENST00000426466)
5. Gene ID (ENSG00000163098)
6. Gene type (pc)
7. Transcript type (pc)
8. Strand (-)
9. Stop ID (Stop2.10.1.2)
10. Stop codon affected (TGA)
11. Base affected (C)
12. Original length (9)
13. Extended length (54)
14. Affected allele (A)
15. Alternative allele (G)
16. SNP id (rs34092035)
17. AFR frequency (0.13)
18. EUR frequency (0)
19. ASN frequency (0)
20. AMR frequency (0.02)
21. Maximum frequency difference (0.13)

Be aware that fields may be added/changed, depending on the arguments used to run stop2go.
