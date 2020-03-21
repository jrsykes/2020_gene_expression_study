#!/bin/bash
#SBATCH --partition=medium
#SBATCH --time=0-01:00:00
#SBATCH --nodes=1
#SBATCH --mem=10gb
#SBATCH --ntasks=6
#SBATCH --output=/home/sykesj/scripts/StdOut/R-%x.%j.out
#SBATCH --error=/home/sykesj/scripts/StdOut/R-%x.%j.err


SPECIES=$1
SRR=$2
SEX=$3
LAYOUT=$4



download_QC () {
	fastq-dump.2.10.3 --outdir /projects/sykesj/raw/$SPECIES/$SEX --defline-seq '@$sn[_$rn]/$ri' --split-files $SRR
	
	/home/sykesj/software/FastQC/fastqc --outdir /projects/sykesj/analyses/$SPECIES/fastqc /projects/sykesj/raw/$SPECIES/$SEX/$SRR\_1.fastq
	/home/sykesj/software/FastQC/fastqc --outdir /projects/sykesj/analyses/$SPECIES/fastqc /projects/sykesj/raw/$SPECIES/$SEX/$SRR\_2.fastq
	

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


	mkdir /scratch/projects/sykesj/BUSCO_tmp
	mkdir /scratch/projects/sykesj/kalisto_tmp
}


trim_QC () {

	#### SINGLE END MODE ####

	if [ $LAYOUT == "SINGLE" ]
	then
	java -jar /home/sykesj/software/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 \
		/projects/sykesj/raw/$SPECIES/$SEX/$SRR\_1.fastq /projects/sykesj/analyses/$SPECIES/trimmomatic/$SEX/$SRR\_s.fq \
		ILLUMINACLIP:/home/sykesj/software/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:12 \
		&& rm -f /projects/sykesj/raw/$SPECIES/$SEX/$SRR*.fastq


	#### PAIRED END MODE ####

	elif [ $LAYOUT == "PAIRED" ]
	then
	java -jar /home/sykesj/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 \
		/projects/sykesj/raw/$SPECIES/$SEX/$SRR\_1.fastq /projects/sykesj/raw/$SPECIES/$SEX/$SRR\_2.fastq \
		/projects/sykesj/analyses/$SPECIES/trimmomatic/$SEX/$SRR\_1.fq /projects/sykesj/analyses/$SPECIES/trimmomatic/$SEX/$SRR\_forward_unpaired.fq.gz \
		/projects/sykesj/analyses/$SPECIES/trimmomatic/$SEX/$SRR\_2.fq /projects/sykesj/analyses/$SPECIES/trimmomatic/$SEX/$SRR\_reverse_unpaired.fq.gz \
		ILLUMINACLIP:/home/sykesj/software/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:12 \
		&& rm -f /projects/sykesj/raw/$SPECIES/$SEX/$SRR*.fastq 

	

	rm -f /projects/sykesj/analyses/$SPECIES/trimmomatic/$SEX/$SRR\_forward_unpaired.fq.gz ; rm -f /projects/sykesj/analyses/$SPECIES/trimmomatic/$SEX/$SRR\_reverse_unpaired.fq.gz

	sed 's/_F\/1/_1/g' /projects/sykesj/analyses/$SPECIES/trimmomatic/$SEX/$SRR\_1.fq | sed 's/_f\/1/_1/g' | sed 's/_forward\/1/_1/g' | sed 's/_Forward\/1/_1/g' 
	sed 's/_R\/2/_2/g' /projects/sykesj/analyses/$SPECIES/trimmomatic/$SEX/$SRR\_2.fq | sed 's/_r\/2/_2/g' | sed 's/_reverse\/2/_2/g' | sed 's/_Reverse\/2/_2/g' 

	else echo 'Error in mode selection at command line'
	fi

	
}




## one = SPECIES, two = SRR, three = SEX, four = paired or single LAYOUT

make_dirs $SPECIES
download_QC $SPECIES $SRR $SEX $LAYOUT
trim_QC $SPECIES $SRR $SEX $LAYOUT

TRIMMED_LIBS=$(for file in $(ls /projects/sykesj/analyses/$SPECIES/trimmomatic/male/*.fq /projects/sykesj/analyses/$SPECIES/trimmomatic/female/*.fq); do readlink -f $file; done | paste -sd "," - )
		echo $TRIMMED_LIBS > /projects/sykesj/analyses/$SPECIES/trinity/fastQCpath.txt

/home/sykesj/software/FastQC/fastqc --outdir $TRIMMED_LIBS

