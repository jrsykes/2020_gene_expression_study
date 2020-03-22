#!/bin/bash
#SBATCH --partition=medium
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --mem=40gb
#SBATCH --ntasks=6
#SBATCH --output=/home/sykesj/scripts/StdOut/R-%x.%j.out
#SBATCH --error=/home/sykesj/scripts/StdOut/R-%x.%j.err


SPECIES=$1
SRR=$2
SEX=$3
LAYOUT=$4


kallisto_map () {


	mkdir /projects/sykesj/analyses/$SPECIES/kallisto/$SRR
	mkdir /scratch/projects/sykesj/map_$SRR

	if [ $LAYOUT == 'PAIRED' ]
	then

		kallisto quant -t 16 -i /projects/sykesj/analyses/$SPECIES/kallisto/paired_$SPECIES.idx -o /scratch/projects/sykesj/map_$SRR/$SRR \
			-b 100 /projects/sykesj/analyses/$SPECIES/trimmomatic/$SEX/$SRR\_1.fq /projects/sykesj/analyses/$SPECIES/trimmomatic/$SEX/$SRR\_2.fq


	elif [ $LAYOUT == 'SINGLE' ]
	then
		READ_LENGTH=$(awk 'BEGIN { t=0.0;sq=0.0; n=0;} ;NR%4==2 {n++;L=length($0);t+=L;sq+=L*L;}END{m=t/n;printf("%f\n",m);}'  /projects/sykesj/analyses/$SPECIES/trimmomatic/$SEX/$SRR\_s.fq)
		SD=$(awk 'BEGIN { t=0.0;sq=0.0; n=0;} ;NR%4==2 {n++;L=length($0);t+=L;sq+=L*L;}END{m=t/n;printf("%f\n",sq/n-m*m);}' /projects/sykesj/analyses/$SPECIES/trimmomatic/$SEX/$SRR\_s.fq)

		kallisto quant -t 16 -i /projects/sykesj/analyses/$SPECIES/kallisto/paired_$SPECIES.idx -o /scratch/projects/sykesj/map_$SRR/$SRR -b 100 \
			--single -l $READ_LENGTH -s $SD /projects/sykesj/analyses/$SPECIES/trimmomatic/$SEX/$SRR\_s.fq

	fi

	rsync -a /scratch/projects/sykesj/map_$SRR/$SRR /projects/sykesj/analyses/$SPECIES/kallisto && rm -rf /scratch/projects/sykesj/map_$SRR

####### setting up files for sleuth #########

	mkdir /projects/sykesj/analyses/$SPECIES/kallisto/kal_results
	touch /projects/sykesj/analyses/$SPECIES/kallisto/kal_results/hiseq_info.txt

	mkdir /projects/sykesj/analyses/$SPECIES/kallisto/kal_results/kal_files
	ln -s /projects/sykesj/analyses/$SPECIES/kallisto/$SRR /projects/sykesj/analyses/$SPECIES/kallisto/kal_results/kal_files/$SRR


}

kallisto_map $SPECIES $SRR $SEX $LAYOUT 