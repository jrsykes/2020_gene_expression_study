#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=10gb
#SBATCH --ntasks=1
#SBATCH -o StdOut-%


species=$1
SRR=$2
sex=$3
mode=$4
fastq-dump.2.10.3 --outdir /home/sykesj/raw/$species/$sex --defline-seq '@$sn[_$rn]/$ri' --split-files $SRR

if [ $mode == 'paired' ]
then
/home/sykesj/software/FastQC/fastqc --outdir /home/sykesj/analyses/$species/fastqc /home/sykesj/raw/$species/$sex/$SRR\_1.fastq
/home/sykesj/software/FastQC/fastqc --outdir /home/sykesj/analyses/$species/fastqc /home/sykesj/raw/$species/$sex/$SRR\_2.fastq

elif 
[ $mode == 'single' ]
then
/home/sykesj/software/FastQC/fastqc --outdir /home/sykesj/analyses/$species/fastqc /home/sykesj/raw/$species/$sex/$SRR\_1.fastq
fi


fastq-dump.2.10.3 --outdir /home/sykesj/raw/veroa_destructor/female --defline-seq '@$sn[_$rn]/$ri' --split-files SRR8100122
