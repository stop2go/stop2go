#!/bin/bash

if [ "$#" == "0" ]; then
		echo "stop2go is a program designed to check for polymorphisms on stop codons that change them to a coding state."
		echo ""
		echo "This program needs the following arguments to run:"
		echo "		--input     || -i:   input file for the program to work. List of transcripts for the program to work. Ony ENSEMBL transcript IDs in the form ENST00000426466 are admitted."
		echo "		--annot     || -an:  annotation file. The annotation format admitted is GTF format."
		echo "		--folder    || -fo:  folder where FASTA files with the sequence of the chromosomes are stored."
		echo "		--freq      || -fr:  VCF file containing the frequences."
		echo ""
		echo "stop2go also have the following optional commands:"
		echo "		--output    || -o  : name for the output file. If this command is not used the name will be the same as the input file."
		echo "		--freqtag          : input your custom tags for the alternative frequencies on your VCF file."
		echo "		--ancestral || -anc: file to search for the ancestral allele. "
		echo "		--anctag           : input your custom tags for the ancestral allele." 
		exit 1
	fi

ARRAY=("$@")
ELEMENTS=${#ARRAY[@]}

BEDTOOLS=$(command -v bedtools)
if [ "$BEDTOOLS" == "" ]; then 
	echo "stop2go needs bedtools to run. Please install bedtools and run stop2go."
	exit 1
fi

EIND=0
ETAG=""
ETYPE=""
OUTPUT=""
FREQTAG=""
ANCTAG=""
QTLS=""
QTLG=""
QTLQ=""
QIND=0
GWASD=0
GWASS=""
GWASN=""
GIND=0

for (( i=0;i<$ELEMENTS;i=i+2)); do
	CHECK=${ARRAY[${i}]}
	j=$(($i+1))

	if [ "$CHECK" = "--extract" ]; then
		EXTRACT=${ARRAY[${j}]}
	fi

	if [ "$CHECK" = "--extract_tag" ]; then
		ETAG=${ARRAY[${j}]}
		EIND=$(($EIND+1))
	fi

	if [ "$CHECK" = "--extract_type" ]; then
		ETYPE=${ARRAY[${j}]}
		EIND=$(($EIND+2))
	fi
	
	if [ "$CHECK" = "--input" ] || [ "$CHECK" = "-i" ]; then
		INPUT=${ARRAY[${j}]}
		MODE=1
	fi

	if [ "$CHECK" = "--add" ]; then
		INPUT=${ARRAY[${j}]}
		MODE=2
	fi

	if [ "$CHECK" = "--output" ] || [ "$CHECK" = "-o" ]; then
		OUTPUT=${ARRAY[${j}]}
	fi

	if [ "$CHECK" = "--annot" ] || [ "$CHECK" = "-an" ]; then
		ANNOT=${ARRAY[${j}]}
	fi

	if [ "$CHECK" = "--folder" ] || [ "$CHECK" = "-fo" ]; then
		FOLDER=${ARRAY[${j}]}
	fi

	if [ "$CHECK" = "--freq" ] || [ "$CHECK" = "-fr" ]; then
		FREQ=${ARRAY[${j}]}
	fi

	if [ "$CHECK" = "--freqtag" ]; then
		FREQTAG=${ARRAY[${j}]}
	fi

	if [ "$CHECK" = "--ancestral" ] || [ "$CHECK" = "-anc" ]; then
		ANCESTRAL=${ARRAY[${j}]}
	fi

	if [ "$CHECK" = "--anctag" ]; then
		ANCTAG=${ARRAY[${j}]}
	fi

	if [ "$CHECK" = "--qtl" ]; then
			QTL=${ARRAY[${j}]}
	fi

	if [ "$CHECK" = "--qtl_spos" ]; then
			QTLS=${ARRAY[${j}]}
			QIND=$(($QIND+1))
	fi

	if [ "$CHECK" = "--qtl_gpos" ]; then
			QTLG=${ARRAY[${j}]}
			QIND=$(($QIND+2))
	fi

	if [ "$CHECK" = "--qtl_qpos" ]; then
			QTLQ=${ARRAY[${j}]}
			QIND=$(($QIND+4))
	fi

	if [ "$CHECK" = "--gwas" ]; then
		GWAS=${ARRAY[${j}]}
	fi

	if [ "$CHECK" = "--gwas_dist" ]; then
		GWASD=${ARRAY[${j}]}
		GIND=$(($GIND+1))
	fi

	if [ "$CHECK" = "--gwas_spos" ]; then
		GWASS=${ARRAY[${j}]}
		GIND=$(($GIND+2))
	fi

	if [ "$CHECK" = "--gwas_npos" ]; then
		GWASN=${ARRAY[${j}]}
		GIND=$(($GIND+4))
	fi

	if [ "$CHECK" = "--help" ] || [ "$CHECK" = "-h" ]; then
		echo "stop2go is a program designed to check for polymorphisms on stop codons that change them to a coding state."
		echo ""
		echo "This program needs the following arguments to run:"
		echo "		--input     || -i:   input file for the program to work. List of transcripts for the program to work. Ony ENSEMBL transcript IDs in the form ENST00000426466 are admitted."
		echo "		--annot     || -an:  annotation file. The annotation format admitted is GTF format."
		echo "		--folder    || -fo:  folder where FASTA files with the sequence of the chromosomes are stored."
		echo "		--freq      || -fr:  VCF file containing the frequences."
		echo ""
		echo "stop2go also have the following optional commands:"
		echo "		--output    || -o  : name for the output file. If this command is not used the name will be the same as the input file."
		echo "		--freqtag          : input your custom tags for the alternative frequencies on your VCF file."
		echo "		--ancestral || -anc: file to search for the ancestral allele. "
		echo "		--anctag           : input your custom tags for the ancestral allele." 
		exit 1
	fi

done

if [ ! -z $EXTRACT ]; then 
	perl ./bin/extract.pl $EXTRACT $EIND $ETAG $ETYPE > id_list.txt
	sort id_list.txt | uniq > temp
	mv temp id_list.txt
	exit 1
fi

if [ "$MODE" = 1 ]; then
	perl ./bin/stop_searcher.pl $INPUT $ANNOT $FOLDER $OUTPUT > $OUTPUT\_list.txt

	$BEDTOOLS intersect -a $FREQ -b $OUTPUT.bed -wb > $OUTPUT\_int_1kg.txt
	cat $OUTPUT\_int_1kg.txt | grep SNP > $OUTPUT\_int_1kg_SNP.txt
	perl ./bin/1kgen_pars.pl $OUTPUT\_int_1kg_SNP.txt $FREQTAG > $OUTPUT\_1kg_SNP_freq.txt
	perl ./bin/stp2stp_out.pl $OUTPUT\_1kg_SNP_freq.txt > res_$OUTPUT.txt

	rm $OUTPUT\_int_1kg.txt
	rm $OUTPUT\_int_1kg_SNP.txt
	rm $OUTPUT\_1kg_SNP_freq.txt
fi

if [ ! -z $ANCESTRAL ]; then 
	if [ "$MODE" = 1 ]; then
		bedtools intersect -a $ANCESTRAL -b res_$OUTPUT.txt -wb > res_$OUTPUT\_anc.txt
		perl ./bin/ancest_pars.pl res_$OUTPUT\_anc.txt $ANCTAG > res_$OUTPUT.txt
		rm res_$OUTPUT\_anc.txt
	else
		bedtools intersect -a $ANCESTRAL -b $INPUT -wb > temp
		perl ./bin/ancest_pars.pl temp $ANCTAG > $INPUT
		rm temp
	fi
fi

if [ ! -z $QTL ]; then 
	if [ "$MODE" = 1 ]; then
		perl ./bin/eQTL_searcher.pl res_$OUTPUT.txt $QTL $QIND $QTLS $QTLG $QTLQ > temp
		mv temp res_$OUTPUT.txt
	else
		perl ./bin/eQTL_searcher.pl $INPUT $QTL $QIND $QTL $QTLS $QTLG $QTLQ > temp
		mv temp $INPUT
	fi
fi

if [ ! -z $GWAS ]; then 
	if [ "$MODE" = 1 ]; then
		perl ./bin/gwas_dist_searcher.pl res_$OUTPUT.txt $GWAS $GIND $GWASD $GWASS $GWASN > temp
		mv temp res_$OUTPUT.txt
	else
		perl ./bin/gwas_dist_searcher.pl $INPUT $GWAS $GIND $GWASD $GWASS $GWASN > temp
		mv temp $INPUT
	fi
fi
