#!/bin/bash
#SBATCH --partition=long
#SBATCH --time=UNLIMITED
#SBATCH --nodes=1
#SBATCH --mem=125gb
#SBATCH --ntasks=20
#SBATCH --output=/projects/sykesj/StdOut/R-%x.%j-Trinity.out
#SBATCH --error=/projects/sykesj/StdOut/R-%x.%j-Trinity.err


SPECIES=$1
LAYOUT=$2
WD=$3

#multi_qc () {
#	if [ $LAYOUT == 'PAIRED' ]
#	then
#		TRIMMED_LIBS=$(for file in $(ls "$WD"/analyses/$SPECIES/trimmomatic/PAIRED/male/*.fq "$WD"/analyses/$SPECIES/trimmomatic/PAIRED/female/*.fq) ; do readlink -f $file; done | paste -sd " " - )
#	elif [ $LAYOUT == 'SINGLE' ]
#	then
#		TRIMMED_LIBS=$(for file in $(ls "$WD"/analyses/$SPECIES/trimmomatic/SINGLE/male/*.fq "$WD"/analyses/$SPECIES/trimmomatic/SINGLE/female/*.fq \
#			"$WD"/analyses/$SPECIES/trimmomatic/PAIRED/male/*.fq "$WD"/analyses/$SPECIES/trimmomatic/PAIRED/female/*.fq) ; do readlink -f $file; done | paste -sd " " - )
#	fi#

#	/home/sykesj/software/FastQC/fastqc --outdir "$WD"/analyses/"$SPECIES"/fastqc2 "$TRIMMED_LIBS" \
#		; multiqc "$WD"/analyses/"$SPECIES"/fastqc/ -o "$WD"/analyses/"$SPECIES"/fastqc/ && rm -f "$WD"/analyses/"$SPECIES"/fastqc/SRR* \
#		; multiqc "$WD"/analyses/"$SPECIES"/fastqc2/ -o "$WD"/analyses/"$SPECIES"/fastqc2/ && rm -f "$WD"/analyses/"$SPECIES"/fastqc2/SRR*
#}

bash /home/sykesj/scripts/2020_gene_expression_study/QC2.sh "$SPECIES" "$LAYOUT"

trinity () {

	mkdir -p "$WD"/analyses/"$SPECIES"/trinity/trinity_"$SPECIES"_"$LAYOUT"


	if [ "$LAYOUT" == 'PAIRED' ]
	then

########### paired ###########

		LEFT=$(for file in $(ls "$WD"/analyses/"$SPECIES"/trimmomatic/"$LAYOUT"/male/*_1.fq "$WD"/analyses/"$SPECIES"/trimmomatic/"$LAYOUT"/female/*_1.fq); \
		do readlink -f "$file"; done | paste -sd "," - )
		echo $LEFT > "$WD"/analyses/"$SPECIES"/trinity/"$LAYOUT"_path.txt
		
		RIGHT=$(for file in $(ls "$WD"/analyses/"$SPECIES"/trimmomatic/"$LAYOUT"/male/*_2.fq "$WD"/analyses/"$SPECIES"/trimmomatic/"$LAYOUT"/female/*_2.fq); \
		do readlink -f "$file"; done | paste -sd "," - )
		echo $RIGHT >> "$WD"/analyses/"$SPECIES"/trinity/"$LAYOUT"_path.txt

	
		/home/sykesj/software/trinityrnaseq-v2.9.1/Trinity --SS_lib_type RF --seqType fq --left "$LEFT" --right "$RIGHT" --CPU 20 --max_memory 125G --output \
			"$WD"/analyses/"$SPECIES"/trinity/trinity_"$SPECIES"_"$LAYOUT" && \
			mv "$WD"/analyses/"$SPECIES"/trinity/trinity_"$SPECIES"_"$LAYOUT"/Trinity.fasta "$WD"/analyses/"$SPECIES"/trinity/"$LAYOUT"_assembly.fa && \
			rm -rf "$WD"/analyses/"$SPECIES"/trinity/trinity_"$SPECIES"_"$LAYOUT"


	elif [ "$LAYOUT" == 'SINGLE' ]
	then

########### single ############

		INPUT=$(for file in $(ls "$WD"/analyses/$SPECIES/trimmomatic/"$LAYOUT"/male/*.fq "$WD"/analyses/"$SPECIES"/trimmomatic/"$LAYOUT"/female/*.fq); \
		do readlink -f "$file"; done | paste -sd "," - )
		echo $INPUT > "$WD"/analyses/"$SPECIES"/trinity/"$LAYOUT"_path.txt

	
		/home/sykesj/software/trinityrnaseq-v2.9.1/Trinity --seqType fq --single "$INPUT" --CPU 20 --max_memory 125G --output \
			"$WD"/analyses/"$SPECIES"/trinity/trinity_"$SPECIES"_"$LAYOUT" && \
			mv "$WD"/analyses/"$SPECIES"/trinity/trinity_"$SPECIES"_"$LAYOUT"/Trinity.fasta "$WD"/analyses/"$SPECIES"/trinity/"$LAYOUT"_assembly.fa && \
			rm -rf "$WD"/analyses/"$SPECIES"/trinity/trinity_"$SPECIES"_"$LAYOUT"


	fi

}



filter_busco_blast_index () {

###### filter out < 1000 bp ##########

	python3 /home/sykesj/scripts/2020_gene_expression_study/1k_filter.py "$WD"/analyses/"$SPECIES"/trinity/"$LAYOUT"_assembly.fa \
		"$WD"/analyses/"$SPECIES"/trinity/"$SPECIES"_"$LAYOUT"_assembly_1k.fa && rm -f "$WD"/analyses/"$SPECIES"/trinity/"$LAYOUT"_assembly.fa

###### busco ##########

	python /home/sykesj/software/busco-master/src/busco/run_BUSCO.py -f --config /home/sykesj/software/busco-master/config/config.ini -i \
		"$WD"/analyses/"$SPECIES"/trinity/"$SPECIES"_"$LAYOUT"_assembly_1k.fa -o busco_"$LAYOUT"_"$SPECIES" \
		-l arthropoda_odb10 -m tran -c 20 \
		&& mv /projects/sykesj/BUSCO_tmp/busco_"$LAYOUT"_"$SPECIES"/short_summary.specific.arthropoda_odb10.busco_"$LAYOUT"_"$SPECIES".txt \
		"$WD"/analyses/"$SPECIES"/busco/BUSCO_out_"$SPECIES"_"$LAYOUT".txt \
		&& rm -rf /projects/sykesj/BUSCO_tmp/busco_"$LAYOUT"_"$SPECIES"

			
			
########## blast ##########

	export BLASTDB=:/home/sykesj/software/blastdb/nt/
	/home/sykesj/software/ncbi-blast-2.10.0+/bin/blastn -task megablast -query "$WD"/analyses/"$SPECIES"/trinity/"$SPECIES"_"$LAYOUT"_assembly_1k.fa -db nt -outfmt '6 qseqid staxids bitscore std' \
		-culling_limit 5 -num_threads 20 -evalue 1e-25 -out "$WD"/analyses/"$SPECIES"/blast/blastn_"$LAYOUT"_"$SPECIES".out \
		&& sort -k 13,13 -n "$WD"/analyses/"$SPECIES"/blast/blastn_"$LAYOUT"_"$SPECIES".out > "$WD"/analyses/"$SPECIES"/blast/"$SPECIES"_blastn_"$LAYOUT"_sorted.out \
		&& rm -f "$WD"/analyses/"$SPECIES"/blast/blastn_"$LAYOUT"_"$SPECIES".out

######### indexing #########

		mkdir -p "$WD"/analyses/"$SPECIES"/kallisto
		cd "$WD"/analyses/"$SPECIES"/kallisto/
		kallisto index -i "$LAYOUT"_"$SPECIES".idx "$WD"/analyses/"$SPECIES"/trinity/"$SPECIES"_"$LAYOUT"_assembly_1k.fa
		cd /home/sykesj

}


multiqc $SPECIES $LAYOUT
trinity $SPECIES $LAYOUT
filter_busco_blast_index $SPECIES $LAYOUT


