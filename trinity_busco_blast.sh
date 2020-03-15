#!/bin/bash
#SBATCH --partition=long
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --mem=100gb
#SBATCH --ntasks=20
#SBATCH --output=/home/sykesj/scripts/StdOut/R-%x.%j.out
#SBATCH --error=/home/sykesj/scripts/StdOut/R-%x.%j.err


SPECIES=$1
LAYOUT=$2

trinity_busco_blast () {



	mkdir /scratch/projects/sykesj/trinity_$SPECIES


	if [ $LAYOUT == 'PAIRED' ]
	then

########### paired ###########

		LEFT=$(for file in $(ls /projects/sykesj/analyses/$SPECIES/trimmomatic/male/*_1.fq /projects/sykesj/analyses/$SPECIES/trimmomatic/female/*_1.fq); do readlink -f $file; done | paste -sd "," - )
		echo $LEFT > /projects/sykesj/analyses/$SPECIES/trinity/path.txt
		RIGHT=$(for file in $(ls /projects/sykesj/analyses/$SPECIES/trimmomatic/male/*_2.fq /projects/sykesj/analyses/$SPECIES/trimmomatic/female/*_2.fq); do readlink -f $file; done | paste -sd "," - )
		echo $RIGHT >> /projects/sykesj/analyses/$SPECIES/trinity/path.txt

	
		/home/sykesj/software/trinityrnaseq-v2.9.1/Trinity --SS_lib_type RF --seqType fq --left $LEFT --right $RIGHT --CPU 20 --max_memory 100G --output /scratch/projects/sykesj/paired_trinity_$SPECIES && rsync -a /scratch/projects/sykesj/paired_trinity_$SPECIES/Trinity.fasta /projects/sykesj/analyses/$SPECIES/trinity/paired_assembly.fa

###### filter < 1000 bp ##########

		/home/sykesj/software/kentUtils/src/utils/faFilter/faFilter.c -minSize=1000 /projects/sykesj/analyses/$SPECIES/trinity/paired_assembly.fa /projects/sykesj/analyses/$SPECIES/trinity/paired_assembly_1k.fa && rm -f /projects/sykesj/analyses/$SPECIES/trinity/paired_assembly.fa

###### busco ##########

		python /home/sykesj/software/busco-master/src/busco/run_BUSCO.py -f --config /home/sykesj/software/busco-master/config/config.ini -i /projects/sykesj/analyses/$SPECIES/trinity/paired_assembly_1k.fa -o busco_paired_$SPECIES \
		-l arthropoda_odb10 -m tran -c 16 \
		&& mkdir /projects/sykesj/analyses/$SPECIES/busco/busco_paired_summaries \
		&& mv /projects/sykesj/analyses/temp_out/trinity/run_busco_paired_$SPECIES /projects/sykesj/analyses/$SPECIES/busco/busco_paired_summaries \
		&& python /home/sykesj/software/busco-master/scripts/generate_plot.py -wd /projects/sykesj/analyses/$SPECIES/busco/busco_paired_summaries/run_busco_paired_$SPECIES \
		&& mv /projects/sykesj/analyses/$SPECIES/busco/busco_paired_summaries/run_busco_paired_$SPECIES/busco_figure.png /projects/sykesj/analyses/$SPECIES/busco/busco_paired_$SPECIES.png \
		&& rm -rf /projects/sykesj/analyses/$SPECIES/busco/busco_paired_summaries \
		&& mv /home/sykesj/software/busco*.log /projects/sykesj/analyses/$SPECIES/busco/

########## blast ##########

		mkdir /projects/sykesj/analyses/$SPECIES/blast 


		blastn -task megablast -query /projects/sykesj/analyses/$SPECIES/trinity/paired_assembly_1k.fa -db /projects/sykesj/blast_db -outfmt '6 qseqid staxids bitscore std' \
			-culling_limit 5 -num_threads 20 -evalue 1e-25 -out /scratch/projects/sykesj/blastn_paired_$SPECIES.out \
			&& mv /scratch/projects/sykesj/blastn_paired_$SPECIES.out /projects/sykesj/analyses/$SPECIES/blast/ \
			&& sort -k 13,13 -n /projects/sykesj/analyses/$SPECIES/blast/blastn_paired_$SPECIES.out > /projects/sykesj/analyses/$SPECIES/blast/$SPECIES\_blastn_paired_sorted.out \
			&& rm -f /projects/sykesj/analyses/$SPECIES/blast/blastn_paired_$SPECIES.out && touch /projects/sykesj/analyses/$SPECIES/blast/paired_blast_$SPECIES.done

######### indexing #########

		kallisto index -i paired_$SPECIES.idx /projects/sykesj/analyses/$SPECIES/trinity/paired_assembly_1k.fa \
			&& mv /projects/sykesj/analyses/temp_out/trinity/paired_$SPECIES.idx /projects/sykesj/analyses/$SPECIES/kallisto/

	elif [ $LAYOUT == 'SINGLE' ]
	then

########### single ############

		INPUT=$(for file in $(ls /projects/sykesj/analyses/$SPECIES/trimmomatic/male/*\_s.fq /projects/sykesj/analyses/$SPECIES/trimmomatic/female/*\_.fq); do readlink -f $file; done | paste -sd "," - )
		echo $INPUT > /projects/sykesj/analyses/$SPECIES/trinity/path.txt

	
		/home/sykesj/software/trinityrnaseq-v2.9.1/Trinity --seqType fq --single $INPUT --CPU 32 --max_memory 100G --output /scratch/projects/sykesj/single_trinity_$SPECIES \
			&& rsync -a /scratch/projects/sykesj/single_trinity_$SPECIES/Trinity.fasta /projects/sykesj/analyses/$SPECIES/trinity/single_assembly.fa


###### filter < 1000 bp ##########

		/home/sykesj/software/kentUtils/src/utils/faFilter/faFilter.c -minSize=1000 /projects/sykesj/analyses/$SPECIES/trinity/single_assembly.fa \
			/projects/sykesj/analyses/$SPECIES/trinity/single_assembly_1k.fa && rm -f /projects/sykesj/analyses/$SPECIES/trinity/single_assembly.fa

###### busco ##########

		python /home/sykesj/software/busco-master/src/busco/run_BUSCO.py -f --config /home/sykesj/software/busco-master/config/config.ini -i /projects/sykesj/analyses/$SPECIES/trinity/single_assembly_1k.fa -o busco_single_$SPECIES \
			-l arthropoda_odb10 -m tran -c 16 \
			&& mkdir /projects/sykesj/analyses/$SPECIES/busco/busco_single_summaries \
			&& mv /projects/sykesj/analyses/temp_out/trinity/run_busco_single_$SPECIES /projects/sykesj/analyses/$SPECIES/busco/busco_single_summaries \
			&& python /home/sykesj/software/busco-master/scripts/generate_plot.py -wd /projects/sykesj/analyses/$SPECIES/busco/busco_single_summaries/run_busco_single_$SPECIES \
			&& mv /projects/sykesj/analyses/$SPECIES/busco/busco_single_summaries/run_busco_single_$SPECIES/busco_figure.png /projects/sykesj/analyses/$SPECIES/busco/busco_single_$SPECIES.png \
			&& rm -rf /projects/sykesj/analyses/$SPECIES/busco/busco_single_summaries \
			&& mv /home/sykesj/software/busco*.log /projects/sykesj/analyses/$SPECIES/busco/


########## blast ##########

		mkdir /projects/sykesj/analyses/$SPECIES/blast

	
		blastn -task megablast -query /projects/sykesj/analyses/$SPECIES/trinity/single_assembly_1k.fa -db /projects/sykesj/blast_db -outfmt '6 qseqid staxids bitscore std' \
			-culling_limit 5 -num_threads 20 -evalue 1e-25 -out /scratch/projects/sykesj/blastn_single_$SPECIES.out \
			&& mv /scratch/projects/sykesj/blastn_single_$SPECIES.out /projects/sykesj/analyses/$SPECIES/blast/ \
			&& sort -k 13,13 -n /projects/sykesj/analyses/$SPECIES/blast/blastn_single_$SPECIES.out > /projects/sykesj/analyses/$SPECIES/blast/$SPECIES\_blastn_single_sorted.out \
			&& rm -f /projects/sykesj/analyses/$SPECIES/blast/blastn_single_$SPECIES.out && touch /projects/sykesj/analyses/$SPECIES/blast/single_blast_$SPECIES.done

######### indexing #########

		kallisto index -i single_$SPECIES.idx /projects/sykesj/analyses/$SPECIES/trinity/single_assembly_1k.fa \
			&& mv /projects/sykesj/analyses/temp_out/trinity/single_$SPECIES.idx /projects/sykesj/analyses/$SPECIES/kallisto/

	fi

	rm -rf /scratch/projects/sykesj/*trinity_$SPECIES

}


multiqc /projects/sykesj/analyses/$SPECIES/fastqc/ -o /projects/sykesj/analyses/$SPECIES/fastqc/ && rm -f /projects/sykesj/analyses/$SPECIES/fastqc/SRR*
multiqc /projects/sykesj/analyses/$SPECIES/fastqc2/ -o /projects/sykesj/analyses/$SPECIES/fastqc2/ && rm -f /projects/sykesj/analyses/$SPECIES/fastqc2/SRR*


trinity_busco_blast $SPECIES $LAYOUT



