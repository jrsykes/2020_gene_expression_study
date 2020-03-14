#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=10gb
#SBATCH --ntasks=6
#SBATCH -o StdOut-%


SPECIES=$1
SRR=$2
SEX=$3
LAYOUT=$4



download_QC () {
	fastq-dump.2.10.3 --outdir /projects/sykesj/raw/$SPECIES/$SEX --defline-seq '@$sn[_$rn]/$ri' --split-files $SRR

	
	/home/sykesj/software/FastQC/fastqc --outdir /projects/sykesj/analyses/$SPECIES/fastqc /projects/sykesj/raw/$SPECIES/$SEX/$SRR\_1.fastq
	/home/sykesj/software/FastQC/fastqc --outdir /projects/sykesj/analyses/$SPECIES/fastqc /projects/sykesj/raw/$SPECIES/$SEX/$SRR\_2.fastq
	

	multiqc /projects/sykesj/analyses/$SPECIES/fastqc/ -o /projects/sykesj/analyses/$SPECIES/fastqc/

	rm -f /projects/sykesj/analyses/$SPECIES/fastqc/SRR*

}


make_dirs () {
	mkdir /projects/sykesj/raw/$SPECIES
	mkdir /projects/sykesj/raw/$SPECIES/male
	mkdir /projects/sykesj/raw/$SPECIES/female

	mkdir /projects/sykesj/analyses/$SPECIES
	mkdir /projects/sykesj/analyses/$SPECIES/fastqc
	mkdir /projects/sykesj/analyses/$SPECIES/fastqc2

	mkdir /projects/sykesj/analyses/$SPECIES/trinity
	mkdir /projects/sykesj/analyses/$SPECIES/busco
	mkdir /projects/sykesj/analyses/$SPECIES/kallisto

	mkdir /projects/sykesj/analyses/$SPECIES/trimmomatic
	mkdir /projects/sykesj/analyses/$SPECIES/trimmomatic/male
	mkdir /projects/sykesj/analyses/$SPECIES/trimmomatic/female
}


trim_QC () {

	#### SINGLE END MODE ####

	if [ $LAYOUT == "SINGLE" ]
	then
	java -jar /home/sykesj/software/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 \
		/projects/sykesj/raw/$SPECIES/$SEX/$SRR\_1.fastq /projects/sykesj/analyses/$SPECIES/trimmomatic/$SEX/$SRR\_s.fq \
		ILLUMINACLIP:/home/sykesj/software/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:12 \
		&& /home/sykesj/software/FastQC/fastqc --outdir /projects/sykesj/analyses/$SPECIES/fastqc2 /projects/sykesj/analyses/$SPECIES/trimmomatic/$SEX/$SRR\_s.fq \
		&& rm -f /projects/sykesj/raw/$SPECIES/$SEX/$SRR\_1.fastq


	#### PAIRED END MODE ####

	elif [ $LAYOUT == "PAIRED" ]
	then
	java -jar /home/sykesj/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 \
		/projects/sykesj/raw/$SPECIES/$SEX/$SRR\_1.fastq /projects/sykesj/raw/$SPECIES/$SEX/$SRR\_2.fastq \
		/projects/sykesj/analyses/$SPECIES/trimmomatic/$SEX/$SRR\_1.fq /projects/sykesj/analyses/$SPECIES/trimmomatic/$SEX/$SRR\_forward_unpaired.fq.gz \
		/projects/sykesj/analyses/$SPECIES/trimmomatic/$SEX/$SRR\_2.fq /projects/sykesj/analyses/$SPECIES/trimmomatic/$SEX/$SRR\_reverse_unpaired.fq.gz \
		ILLUMINACLIP:/home/sykesj/software/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:12 \
		&& /home/sykesj/software/FastQC/fastqc --outdir /projects/sykesj/analyses/$SPECIES/fastqc2 \
		/projects/sykesj/analyses/$SPECIES/trimmomatic/$SEX/$SRR\_1.fq /projects/sykesj/analyses/$SPECIES/trimmomatic/$SEX/$SRR\_2.fq \
		&& rm -f /projects/sykesj/raw/$SPECIES/$SEX/$SRR\_1.fastq && rm -f /projects/sykesj/raw/$SPECIES/$SEX/$SRR\_2.fastq

	

	rm -f /projects/sykesj/analyses/$SPECIES/trimmomatic/$SEX/$SRR\_forward_unpaired.fq.gz ; rm -f /projects/sykesj/analyses/$SPECIES/trimmomatic/$SEX/$SRR\_reverse_unpaired.fq.gz

	sed 's/_F\/1/_1/g' /projects/sykesj/analyses/$SPECIES/trimmomatic/$SEX/$SRR\_1.fq | sed 's/_f\/1/_1/g' | sed 's/_forward\/1/_1/g' | sed 's/_Forward\/1/_1/g' 
	sed 's/_R\/2/_2/g' /projects/sykesj/analyses/$SPECIES/trimmomatic/$SEX/$SRR\_2.fq | sed 's/_r\/2/_2/g' | sed 's/_reverse\/2/_2/g' | sed 's/_Reverse\/2/_2/g' 

	else echo 'Error in mode selection at command line'
	fi

	multiqc /projects/sykesj/analyses/$SPECIES/fastqc2/ -o /projects/sykesj/analyses/$SPECIES/fastqc2/
	rm -f /projects/sykesj/analyses/$SPECIES/fastqc2/SRR*
}


clear_SPECIES_dat () {
	
	rm -rf /projects/sykesj/raw/$SPECIES
		
	rm -rf /projects/sykesj/analyses/$SPECIES

	rm -rf /scratch/projects/sykesj/*$SPECIES*

}

## one = SPECIES, two = SRR, three = SEX, four = paired or single LAYOUT

clear_SPECIES_dat $SPECIES
make_dirs $SPECIES
download_QC $SPECIES $SRR $SEX $LAYOUT
#trim_QC $SPECIES $SRR $SEX $LAYOUT



