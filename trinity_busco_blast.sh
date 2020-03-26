#!/bin/bash
#SBATCH --partition=long
#SBATCH --time=5-00:00:00
#SBATCH --nodes=1
#SBATCH --mem=100gb
#SBATCH --ntasks=20
#SBATCH --output=/home/sykesj/scripts/StdOut/R-%x.%j.out
#SBATCH --error=/home/sykesj/scripts/StdOut/R-%x.%j.err


SPECIES=$1
LAYOUT=$2

trinity_busco_blast () {



	mkdir -p /scratch/projects/sykesj/trinity_$SPECIES_$LAYOUT


	if [ $LAYOUT == 'PAIRED' ]
	then

########### paired ###########

		LEFT=$(for file in $(ls /projects/sykesj/analyses/$SPECIES/trimmomatic/$LAYOUT/male/*_1.fq /projects/sykesj/analyses/$SPECIES/trimmomatic/$LAYOUT/female/*_1.fq); do readlink -f $file; done | paste -sd "," - )
		echo $LEFT > /projects/sykesj/analyses/$SPECIES/trinity/paired_path.txt
		RIGHT=$(for file in $(ls /projects/sykesj/analyses/$SPECIES/trimmomatic/$LAYOUT/male/*_2.fq /projects/sykesj/analyses/$SPECIES/trimmomatic/$LAYOUT/female/*_2.fq); do readlink -f $file; done | paste -sd "," - )
		echo $RIGHT >> /projects/sykesj/analyses/$SPECIES/trinity/paired_path.txt

	
		/home/sykesj/software/trinityrnaseq-v2.9.1/Trinity --SS_lib_type RF --seqType fq --left $LEFT --right $RIGHT --CPU 20 --max_memory 100G --output \
			/scratch/projects/sykesj/trinity_$SPECIES_$LAYOUT && rsync -a /scratch/projects/sykesj/trinity_$SPECIES_$LAYOUT/Trinity.fasta /projects/sykesj/analyses/$SPECIES/trinity/paired_assembly.fa

###### filter out < 1000 bp ##########

		python3 /home/sykesj/scripts/2020_gene_expression_study/1k_filter.py /projects/sykesj/analyses/$SPECIES/trinity/paired_assembly.fa \
			/projects/sykesj/analyses/$SPECIES/trinity/$LAYOUT_assembly_1k.fa && rm -f /projects/sykesj/analyses/$SPECIES/trinity/paired_assembly.fa

###### busco ##########

		python /home/sykesj/software/busco-master/src/busco/run_BUSCO.py -f --config /home/sykesj/software/busco-master/config/config.ini -i \
			/projects/sykesj/analyses/$SPECIES/trinity/paired_assembly_1k.fa -o busco_$LAYOUT_$SPECIES \
			-l arthropoda_odb10 -m tran -c 16 \
			&& mv /scratch/projects/sykesj/BUSCO_tmp/busco_paired_$SPECIES/short_summary.specific.arthropoda_odb10.busco_paired_$SPECIES.txt /projects/sykesj/analyses/$SPECIES/busco/BUSCO_out_$SPECIES_$LAYOUT.txt \
			&& rm -rf /scratch/projects/sykesj/BUSCO_tmp/busco_paired_$SPECIES
			
			
########## blast ##########

		mkdir /projects/sykesj/analyses/$SPECIES/blast 

		export BLASTDB=:/home/sykesj/software/blastdb/nt/
		/home/sykesj/software/ncbi-blast-2.10.0+/bin/blastn -task megablast -query /projects/sykesj/analyses/$SPECIES/trinity/paired_assembly_1k.fa -db nt -outfmt '6 qseqid staxids bitscore std' \
			-culling_limit 5 -num_threads 20 -evalue 1e-25 -out /scratch/projects/sykesj/blastn_paired_$SPECIES.out \
			&& rsync -a /scratch/projects/sykesj/blastn_paired_$SPECIES.out /projects/sykesj/analyses/$SPECIES/blast/ \
			&& sort -k 13,13 -n /projects/sykesj/analyses/$SPECIES/blast/blastn_paired_$SPECIES.out > /projects/sykesj/analyses/$SPECIES/blast/$SPECIES\_blastn_$LAYOUT_sorted.out \
			&& rm -f /projects/sykesj/analyses/$SPECIES/blast/blastn_paired_$SPECIES.out && rm -rf /scratch/projects/sykesj/blastn_paired_$SPECIES.out

######### indexing #########
		
		cd /projects/sykesj/analyses/$SPECIES/kallisto/
		kallisto index -i paired_$SPECIES.idx /projects/sykesj/analyses/$SPECIES/trinity/paired_assembly_1k.fa
		cd /home/sykesj


	elif [ $LAYOUT == 'SINGLE' ]
	then

########### single ############

		INPUT=$(for file in $(ls /projects/sykesj/analyses/$SPECIES/trimmomatic/$LAYOUT/male/*.fq /projects/sykesj/analyses/$SPECIES/trimmomatic/$LAYOUT/female/*.fq); do readlink -f $file; done | paste -sd "," - )
		echo $INPUT > /projects/sykesj/analyses/$SPECIES/trinity/single_path.txt

	
		/home/sykesj/software/trinityrnaseq-v2.9.1/Trinity --seqType fq --single $INPUT --CPU 20 --max_memory 100G --output /scratch/projects/sykesj/trinity_$SPECIES_$LAYOUT \
			&& rsync -a /scratch/projects/sykesj/trinity_$SPECIES_$LAYOUT/Trinity.fasta /projects/sykesj/analyses/$SPECIES/trinity/single_assembly.fa


###### filter < 1000 bp ##########

		python3 /home/sykesj/scripts/2020_gene_expression_study/1k_filter.py /projects/sykesj/analyses/$SPECIES/trinity/single_assembly.fa \
			/projects/sykesj/analyses/$SPECIES/trinity/$LAYOUT_assembly_1k.fa && rm -f /projects/sykesj/analyses/$SPECIES/trinity/single_assembly.fa

###### busco ##########

		python /home/sykesj/software/busco-master/src/busco/run_BUSCO.py -f --config /home/sykesj/software/busco-master/config/config.ini -i \
			/projects/sykesj/analyses/$SPECIES/trinity/single_assembly_1k.fa -o busco_$LAYOUT_$SPECIES \
			-l arthropoda_odb10 -m tran -c 16 \
			&& mv /scratch/projects/sykesj/BUSCO_tmp/busco_single_$SPECIES/short_summary.specific.arthropoda_odb10.busco_single_$SPECIES.txt /projects/sykesj/analyses/$SPECIES/busco/BUSCO_out_$SPECIES_$LAYOUT.txt \
			&& rm -rf /scratch/projects/sykesj/BUSCO_tmp/busco_single_$SPECIES



########## blast ##########

		mkdir /projects/sykesj/analyses/$SPECIES/blast

		export BLASTDB=:/home/sykesj/software/blastdb/nt/
		/home/sykesj/software/ncbi-blast-2.10.0+/bin/blastn -task megablast -query /projects/sykesj/analyses/$SPECIES/trinity/single_assembly_1k.fa -db nt -outfmt '6 qseqid staxids bitscore std' \
			-culling_limit 5 -num_threads 20 -evalue 1e-25 -out /scratch/projects/sykesj/blastn_single_$SPECIES.out \
			&& rsync -a /scratch/projects/sykesj/blastn_single_$SPECIES.out /projects/sykesj/analyses/$SPECIES/blast/ \
			&& sort -k 13,13 -n /projects/sykesj/analyses/$SPECIES/blast/blastn_single_$SPECIES.out > /projects/sykesj/analyses/$SPECIES/blast/$SPECIES\_blastn_$LAYOUT_sorted.out \
			&& rm -f /projects/sykesj/analyses/$SPECIES/blast/blastn_single_$SPECIES.out && rm -rf /scratch/projects/sykesj/blastn_single_$SPECIES.out


######### indexing #########

		cd /projects/sykesj/analyses/$SPECIES/kallisto/
		kallisto index -i single_$SPECIES.idx /projects/sykesj/analyses/$SPECIES/trinity/single_assembly_1k.fa
		cd /home/sykesj
	fi

	

}

rm -rf /scratch/projects/sykesj/BUSCO_tmp/busco_$LAYOUT_$SPECIES


multi_qc () {

	if [ $LAYOUT == 'PAIRED' ]
	then
		TRIMMED_LIBS=$(for file in $(ls /projects/sykesj/analyses/$SPECIES/trimmomatic/PAIRED/male/*.fq /projects/sykesj/analyses/$SPECIES/trimmomatic/PAIRED/female/*.fq) ; \
			do readlink -f $file; done | paste -sd " " - )
	
	elif [ $LAYOUT == 'SINGLE' ]
	then
		TRIMMED_LIBS=$(for file in $(ls /projects/sykesj/analyses/$SPECIES/trimmomatic/SINGLE/male/*.fq /projects/sykesj/analyses/$SPECIES/trimmomatic/SINGLE/female/*.fq \
			/projects/sykesj/analyses/$SPECIES/trimmomatic/PAIRED/male/*.fq /projects/sykesj/analyses/$SPECIES/trimmomatic/PAIRED/female/*.fq) ; do readlink -f $file; done | paste -sd " " - )
	fi
	
	/home/sykesj/software/FastQC/fastqc --outdir /projects/sykesj/analyses/$SPECIES/fastqc2 $TRIMMED_LIBS \
		&& multiqc /projects/sykesj/analyses/$SPECIES/fastqc/ -o /projects/sykesj/analyses/$SPECIES/fastqc/ && rm -f /projects/sykesj/analyses/$SPECIES/fastqc/SRR* \
		&& multiqc /projects/sykesj/analyses/$SPECIES/fastqc2/ -o /projects/sykesj/analyses/$SPECIES/fastqc2/ && rm -f /projects/sykesj/analyses/$SPECIES/fastqc2/SRR*
	
}

rm -rf /scratch/projects/sykesj/*$SPECIES_$LAYOUT*

multiqc $SPECIES $LAYOUT
trinity_busco_blast $SPECIES $LAYOUT

rm -rf /scratch/projects/sykesj/*$SPECIES_$LAYOUT*

