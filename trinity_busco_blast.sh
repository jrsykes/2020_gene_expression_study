#!/bin/bash
#SBATCH --partition=long
#SBATCH --time=5-00:00:00
#SBATCH --nodes=1
#SBATCH --mem=100gb
#SBATCH --ntasks=30
#SBATCH --output=/home/sykesj/scripts/StdOut/R-%x.%j.out
#SBATCH --error=/home/sykesj/scripts/StdOut/R-%x.%j.err


SPECIES=$1
LAYOUT=$2


multi_qc () {

	TRIMMED_LIBS=$(for file in $(ls /projects/sykesj/analyses/"$SPECIES"/trimmomatic/"$LAYOUT"/male/*.fq /projects/sykesj/analyses/"$SPECIES"/trimmomatic/"$LAYOUT"/female/*.fq) ; \
		do readlink -f "$file"; done | paste -sd " " - )
	
	/home/sykesj/software/FastQC/fastqc --outdir /projects/sykesj/analyses/$SPECIES/fastqc2 $TRIMMED_LIBS \
		&& multiqc /projects/sykesj/analyses/$SPECIES/fastqc/ -o /projects/sykesj/analyses/$SPECIES/fastqc/ && rm -f /projects/sykesj/analyses/$SPECIES/fastqc/SRR* \
		&& multiqc /projects/sykesj/analyses/$SPECIES/fastqc2/ -o /projects/sykesj/analyses/$SPECIES/fastqc2/ && rm -f /projects/sykesj/analyses/$SPECIES/fastqc2/SRR*
	
}


trinity_busco_blast () {

	mkdir -p /scratch/projects/sykesj/trinity_"$SPECIES"_"$LAYOUT"


	if [ "$LAYOUT" == 'PAIRED' ]
	then

########### paired ###########

		LEFT=$(for file in $(ls /projects/sykesj/analyses/"$SPECIES"/trimmomatic/"$LAYOUT"/male/*_1.fq /projects/sykesj/analyses/"$SPECIES"/trimmomatic/"$LAYOUT"/female/*_1.fq); \
		do readlink -f "$file"; done | paste -sd "," - )
		echo $LEFT > /projects/sykesj/analyses/"$SPECIES"/trinity/"$LAYOUT"_path.txt
		
		RIGHT=$(for file in $(ls /projects/sykesj/analyses/"$SPECIES"/trimmomatic/"$LAYOUT"/male/*_2.fq /projects/sykesj/analyses/"$SPECIES"/trimmomatic/"$LAYOUT"/female/*_2.fq); \
		do readlink -f "$file"; done | paste -sd "," - )
		echo $RIGHT >> /projects/sykesj/analyses/"$SPECIES"/trinity/"$LAYOUT"_path.txt

	
		/home/sykesj/software/trinityrnaseq-v2.9.1/Trinity --SS_lib_type RF --seqType fq --left "$LEFT" --right "$RIGHT" --CPU 30 --max_memory 100G --output \
			/scratch/projects/sykesj/trinity_"$SPECIES"_"$LAYOUT" && \
			mv /scratch/projects/sykesj/trinity_"$SPECIES"_"$LAYOUT"/Trinity.fasta /projects/sykesj/analyses/"$SPECIES"/trinity/"$LAYOUT"_assembly.fa && \
			rm -rf /scratch/projects/sykesj/trinity_"$SPECIES"_"$LAYOUT"*


	elif [ "$LAYOUT" == 'SINGLE' ]
	then

########### single ############

		INPUT=$(for file in $(ls /projects/sykesj/analyses/$SPECIES/trimmomatic/"$LAYOUT"/male/*.fq /projects/sykesj/analyses/"$SPECIES"/trimmomatic/"$LAYOUT"/female/*.fq); \
		do readlink -f "$file"; done | paste -sd "," - )
		echo $INPUT > /projects/sykesj/analyses/"$SPECIES"/trinity/"$LAYOUT"_path.txt

	
		/home/sykesj/software/trinityrnaseq-v2.9.1/Trinity --seqType fq --single "$INPUT" --CPU 30 --max_memory 100G --output /scratch/projects/sykesj/trinity_"$SPECIES"_"$LAYOUT" && \
			mv /scratch/projects/sykesj/trinity_"$SPECIES"_"$LAYOUT"/Trinity.fasta /projects/sykesj/analyses/"$SPECIES"/trinity/"$LAYOUT"_assembly.fa && \
			rm -rf /scratch/projects/sykesj/trinity_"$SPECIES"_"$LAYOUT"*


	fi

}



filter_busco_blast_index () {

###### filter out < 1000 bp ##########

	python3 /home/sykesj/scripts/2020_gene_expression_study/1k_filter.py /projects/sykesj/analyses/"$SPECIES"/trinity/"$LAYOUT"_assembly.fa \
		/projects/sykesj/analyses/"$SPECIES"/trinity/"$LAYOUT"_assembly_1k.fa && rm -f /projects/sykesj/analyses/"$SPECIES"/trinity/"$LAYOUT"_assembly.fa

###### busco ##########

	python /home/sykesj/software/busco-master/src/busco/run_BUSCO.py -f --config /home/sykesj/software/busco-master/config/config.ini -i \
		/projects/sykesj/analyses/"$SPECIES"/trinity/"$LAYOUT"_assembly_1k.fa -o busco_"$LAYOUT"_"$SPECIES" \
		-l arthropoda_odb10 -m tran -c 16 \
		&& mv /scratch/projects/sykesj/BUSCO_tmp/busco_"$LAYOUT"_"$SPECIES"/short_summary.specific.arthropoda_odb10.busco_"$LAYOUT"_"$SPECIES".txt \
		/projects/sykesj/analyses/"$SPECIES"/busco/BUSCO_out_"$SPECIES"_"$LAYOUT".txt \
		&& rm -rf /scratch/projects/sykesj/BUSCO_tmp/busco_"$LAYOUT"_"$SPECIES"

			
			
########## blast ##########

	export BLASTDB=:/home/sykesj/software/blastdb/nt/
	/home/sykesj/software/ncbi-blast-2.10.0+/bin/blastn -task megablast -query /projects/sykesj/analyses/"$SPECIES"/trinity/"$LAYOUT"_assembly_1k.fa -db nt -outfmt '6 qseqid staxids bitscore std' \
		-culling_limit 5 -num_threads 20 -evalue 1e-25 -out /scratch/projects/sykesj/blastn_"$LAYOUT"_"$SPECIES".out \
		&& mv /scratch/projects/sykesj/blastn_"$LAYOUT"_"$SPECIES".out /projects/sykesj/analyses/"$SPECIES"/blast/ \
		&& sort -k 13,13 -n /projects/sykesj/analyses/"$SPECIES"/blast/blastn_"$LAYOUT"_"$SPECIES".out > /projects/sykesj/analyses/"$SPECIES"/blast/"$SPECIES"_blastn_"$LAYOUT"_sorted.out \
		&& rm -f /projects/sykesj/analyses/"$SPECIES"/blast/blastn_"$LAYOUT"_"$SPECIES".out

######### indexing #########
		
		cd /projects/sykesj/analyses/"$SPECIES"/kallisto/
		kallisto index -i "$LAYOUT"_"$SPECIES".idx /projects/sykesj/analyses/"$SPECIES"/trinity/"$LAYOUT"_assembly_1k.fa
		cd /home/sykesj

}


#multiqc $SPECIES $LAYOUT
trinity_busco_blast $SPECIES $LAYOUT
filter_busco_blast_index $SPECIES $LAYOUT


