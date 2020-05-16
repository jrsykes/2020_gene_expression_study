#!/bin/bash
#SBATCH --partition=medium
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --mem=100gb
#SBATCH --ntasks=6
#SBATCH --output=/scratch/projects/sykesj/StdOut/R-%x.%j-download.out
#SBATCH --error=/scratch/projects/sykesj/StdOut/R-%x.%j-download.err


SPECIES=$1
SRR=$2
SEX=$3
LAYOUT=$4
WD=$5


download_QC () {
	fastq-dump.2.10.3 --outdir "$WD"/raw/$SPECIES/$SEX --defline-seq '@$sn[_$rn]/$ri' --split-files $SRR
	
	/home/sykesj/software/FastQC/fastqc --outdir "$WD"/analyses/$SPECIES/fastqc "$WD"/raw/$SPECIES/$SEX/$SRR\_1.fastq
	/home/sykesj/software/FastQC/fastqc --outdir "$WD"/analyses/$SPECIES/fastqc "$WD"/raw/$SPECIES/$SEX/$SRR\_2.fastq
	

}


make_dirs () {
	mkdir -p "$WD"/raw/$SPECIES/male
	mkdir -p "$WD"/raw/$SPECIES/female

	mkdir -p "$WD"/analyses/$SPECIES
	mkdir -p "$WD"/analyses/$SPECIES/fastqc
	mkdir -p "$WD"/analyses/$SPECIES/fastqc2

	mkdir -p "$WD"/analyses/$SPECIES/trinity
	mkdir -p "$WD"/analyses/$SPECIES/busco
	mkdir "$WD"/analyses/$SPECIES/blast
	
	mkdir "$WD"/analyses/$SPECIES/trimmomatic
	
	mkdir /scratch/"$WD"/BUSCO_tmp
	mkdir /scratch/"$WD"/kalisto_tmp


}


trim_QC () {

	#### SINGLE END MODE ####

	if [ $LAYOUT == "SINGLE" ]
	then

	mkdir "$WD"/analyses/$SPECIES/trimmomatic/$LAYOUT	
	mkdir "$WD"/analyses/$SPECIES/trimmomatic/$LAYOUT/male
	mkdir "$WD"/analyses/$SPECIES/trimmomatic/$LAYOUT/female

	java -jar /home/sykesj/software/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 \
		"$WD"/raw/$SPECIES/$SEX/$SRR\_1.fastq "$WD"/analyses/$SPECIES/trimmomatic/$LAYOUT/$SEX/$SRR\_s.fq \
		ILLUMINACLIP:/home/sykesj/software/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:12 \
		&& rm -f "$WD"/raw/$SPECIES/$SEX/$SRR*.fastq


	#### PAIRED END MODE ####

	elif [ $LAYOUT == "PAIRED" ]
	then

	mkdir "$WD"/analyses/$SPECIES/trimmomatic/$LAYOUT	
	mkdir "$WD"/analyses/$SPECIES/trimmomatic/$LAYOUT/male
	mkdir "$WD"/analyses/$SPECIES/trimmomatic/$LAYOUT/female

	java -jar /home/sykesj/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 \
		"$WD"/raw/$SPECIES/$SEX/$SRR\_1.fastq "$WD"/raw/$SPECIES/$SEX/$SRR\_2.fastq \
		"$WD"/analyses/$SPECIES/trimmomatic/$LAYOUT/$SEX/$SRR\_1.fq "$WD"/analyses/$SPECIES/trimmomatic/$LAYOUT/$SEX/$SRR\_forward_unpaired.fq.gz \
		"$WD"/analyses/$SPECIES/trimmomatic/$LAYOUT/$SEX/$SRR\_2.fq "$WD"/analyses/$SPECIES/trimmomatic/$LAYOUT/$SEX/$SRR\_reverse_unpaired.fq.gz \
		ILLUMINACLIP:/home/sykesj/software/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:12 \
		&& rm -f "$WD"/raw/$SPECIES/$SEX/$SRR*.fastq 

	

	rm -f "$WD"/analyses/$SPECIES/trimmomatic/$LAYOUT/$SEX/$SRR\_forward_unpaired.fq.gz ; rm -f "$WD"/analyses/$SPECIES/trimmomatic/$LAYOUT/$SEX/$SRR\_reverse_unpaired.fq.gz

	sed 's/_F\/1/_1/g' "$WD"/analyses/$SPECIES/trimmomatic/"$LAYOUT"/$SEX/$SRR\_1.fq | sed 's/_f\/1/_1/g' | sed 's/_forward\/1/_1/g' | sed 's/_Forward\/1/_1/g' 
	sed 's/_R\/2/_2/g' "$WD"/analyses/$SPECIES/trimmomatic/"$LAYOUT"/$SEX/$SRR\_2.fq | sed 's/_r\/2/_2/g' | sed 's/_reverse\/2/_2/g' | sed 's/_Reverse\/2/_2/g' 

	else echo 'Error in mode selection at command line'
	fi

	
}




## one = SPECIES, two = SRR, three = SEX, four = paired or single LAYOUT

make_dirs $SPECIES
download_QC $SPECIES $SRR $SEX $LAYOUT
trim_QC $SPECIES $SRR $SEX $LAYOUT


