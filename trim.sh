#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=10gb
#SBATCH --ntasks=4
#SBATCH -o StdOut-%

species=$1
SRR=$2
sex=$3
layout=$4


#### SINGLE END MODE ####

if [ $layout == "single" ]
then
java -jar /home/sykesj/software/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 /projects/sykesj/raw/$species/$sex/$SRR\_1.fastq /projects/sykesj/analyses/$species/trimmomatic/$sex/$SRR\_s.fq ILLUMINACLIP:/home/sykesj/software/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:12 && /home/sykesj/software/FastQC/fastqc --outdir /projects/sykesj/analyses/$species/fastqc2 /projects/sykesj/analyses/$species/trimmomatic/$sex/$SRR\_s.fq && rm -f /projects/sykesj/raw/$species/$sex/$SRR\_1.fastq

#### PAIRED END MODE ####

elif [ $layout == "paired" ]
then
java -jar /home/sykesj/software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 /projects/sykesj/raw/$species/$sex/$SRR\_1.fastq /projects/sykesj/raw/$species/$sex/$SRR\_2.fastq /projects/sykesj/analyses/$species/trimmomatic/$sex/$SRR\_1.fq /projects/sykesj/analyses/$species/trimmomatic/$sex/$SRR\_forward_unpaired.fq.gz /projects/sykesj/analyses/$species/trimmomatic/$sex/$SRR\_2.fq /projects/sykesj/analyses/$species/trimmomatic/$sex/$SRR\_reverse_unpaired.fq.gz ILLUMINACLIP:/home/sykesj/software/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 HEADCROP:12 && /home/sykesj/software/FastQC/fastqc --outdir /projects/sykesj/analyses/$species/fastqc2 /projects/sykesj/analyses/$species/trimmomatic/$sex/$SRR\_1.fq /projects/sykesj/analyses/$species/trimmomatic/$sex/$SRR\_2.fq && rm -f /projects/sykesj/raw/$species/$sex/$SRR\_1.fastq && rm -f /projects/sykesj/raw/$species/$sex/$SRR\_2.fastq

rm -f /projects/sykesj/analyses/$species/trimmomatic/$sex/$SRR\_forward_unpaired.fq.gz ; rm -f /projects/sykesj/analyses/$species/trimmomatic/$sex/$SRR\_reverse_unpaired.fq.gz

sed 's/_F\/1/_1/g' /projects/sykesj/analyses/$species/trimmomatic/$sex/$SRR\_1.fq | sed 's/_f\/1/_1/g' | sed 's/_forward\/1/_1/g' | sed 's/_Forward\/1/_1/g' # > /data/projects/lross_ssa/analyses/$species/trimmomatic/$sex/$SRR\_1.fq && rm -f /data/projects/lross_ssa/analyses/$species/trimmomatic/$sex/$SRR\_11.fq
sed 's/_R\/2/_2/g' /projects/sykesj/analyses/$species/trimmomatic/$sex/$SRR\_2.fq | sed 's/_r\/2/_2/g' | sed 's/_reverse\/2/_2/g' | sed 's/_Reverse\/2/_2/g' # > /data/projects/lross_ssa/analyses/$species/trimmomatic/$sex/$SRR\_2.fq && rm -f /data/projects/lross_ssa/analyses/$species/trimmomatic/$sex/$SRR\_22.fq

else echo 'Error in mode selection at command line'
fi



