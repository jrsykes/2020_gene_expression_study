#!/bin/bash
#$ -S /bin/bash
#$ -q main.q
#$ -pe smp 2
#$ -l h_vmem=10G
#$ -wd /data/projects/lross_ssa/analyses/temp_out/download


species=$1
SRR=$2
sex=$3
mode=$4
/exports/software/sratoolkit/sratoolkit.2.8.2-1-ubuntu64/bin/fastq-dump --outdir /data/projects/lross_ssa/raw/$species/$sex --defline-seq '@$sn[_$rn]/$ri' --split-files $SRR

if [ $mode == 'paired' ]
then
/exports/software/fastqc/fastqc_v0.11.5/FastQC/fastqc --outdir /data/projects/lross_ssa/analyses/$species/fastqc /data/projects/lross_ssa/raw/$species/$sex/$SRR\_1.fastq
/exports/software/fastqc/fastqc_v0.11.5/FastQC/fastqc --outdir /data/projects/lross_ssa/analyses/$species/fastqc /data/projects/lross_ssa/raw/$species/$sex/$SRR\_2.fastq

elif 
[ $mode == 'single' ]
then
/exports/software/fastqc/fastqc_v0.11.5/FastQC/fastqc --outdir /data/projects/lross_ssa/analyses/$species/fastqc /data/projects/lross_ssa/raw/$species/$sex/$SRR\_1.fastq
fi

#mkdir /data/projects/lross_ssa/raw/$species/genome
#wget https://www.ncbi.nlm.nih.gov/genome/14032?genome_assembly_id=59542 -O /data/projects/lross_ssa/raw/$species/genome/$species.fa
#wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/695/645/GCA_000695645.1_Pven_1.0/GCA_000695645.1_Pven_1.0_genomic.gbff.gz -O /data/projects/lross_ssa/raw/$species/genome/$species.gbff
#wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/695/645/GCA_000695645.1_Pven_1.0/GCA_000695645.1_Pven_1.0_genomic.gff.gz -O /data/projects/lross_ssa/raw/$species/genome/$species.gff
#/exports/software/sratoolkit/sratoolkit.2.8.2-1-ubuntu64/bin/fastq-dump --outdir /data/projects/lross_ssa/raw/pachypsylla_venusta/genome GCA_000695645.1

###### The code below is to be used to download libraries useing the wget command. This is much faster than dowoading with the fastq-dump command used above but trinity wil only run on libraries downloaded useing the above fastq-dump command and so unless you are using another assembler such as SOAP the following code should be ignored.

#echo -n $SRR | wc -c > $SRR\n.txt
#read SRRn < $SRR\n.txt ; rm -f $SRR\n.txt
#
#if [ $mode == 'paired' ]
#	then
#if [ $SRRn -eq 9 ]
#	then
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${SRR:0:6}/$SRR/$SRR\_1.fastq.gz -O /data/projects/lross_ssa/raw/$species/$sex/$SRR\_1.fastq.gz && gunzip /data/projects/lross_ssa/raw/$species/$sex/$SRR\_1.fastq.gz && /exports/software/fastqc/fastqc_v0.11.5/FastQC/fastqc --outdir /data/projects/lross_ssa/analyses/$species/fastqc /data/projects/lross_ssa/raw/$species/$sex/$SRR\_1.fastq
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${SRR:0:6}/$SRR/$SRR\_2.fastq.gz -O /data/projects/lross_ssa/raw/$species/$sex/$SRR\_2.fastq.gz && gunzip /data/projects/lross_ssa/raw/$species/$sex/$SRR\_2.fastq.gz && /exports/software/fastqc/fastqc_v0.11.5/FastQC/fastqc --outdir /data/projects/lross_ssa/analyses/$species/fastqc /data/projects/lross_ssa/raw/$species/$sex/$SRR\_2.fastq	
#	elif [ $SRRn -eq 10 ]
#	then
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${SRR:0:6}/00${SRR:9:11}/$SRR/$SRR\_1.fastq.gz -O /data/projects/lross_ssa/raw/$species/$sex/$SRR\_1.fastq.gz && gunzip /data/projects/lross_ssa/raw/$species/$sex/$SRR\_1.fastq.gz && /exports/software/fastqc/fastqc_v0.11.5/FastQC/fastqc --outdir /data/projects/lross_ssa/analyses/$species/fastqc /data/projects/lross_ssa/raw/$species/$sex/$SRR\_1.fastq
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${SRR:0:6}/00${SRR:9:11}/$SRR/$SRR\_2.fastq.gz -O /data/projects/lross_ssa/raw/$species/$sex/$SRR\_2.fastq.gz && gunzip /data/projects/lross_ssa/raw/$species/$sex/$SRR\_2.fastq.gz && /exports/software/fastqc/fastqc_v0.11.5/FastQC/fastqc --outdir /data/projects/lross_ssa/analyses/$species/fastqc /data/projects/lross_ssa/raw/$species/$sex/$SRR\_2.fastq	
#	elif [ $SRRn -eq 11 ]
#	then
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${SRR:0:6}/0${SRR:9:12}/$SRR/$SRR\_1.fastq.gz -O /data/projects/lross_ssa/raw/$species/$sex/$SRR\_1.fastq.gz && gunzip /data/projects/lross_ssa/raw/$species/$sex/$SRR\_1.fastq.gz && /exports/software/fastqc/fastqc_v0.11.5/FastQC/fastqc --outdir /data/projects/lross_ssa/analyses/$species/fastqc /data/projects/lross_ssa/raw/$species/$sex/$SRR\_1.fastq
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${SRR:0:6}/0${SRR:9:12}/$SRR/$SRR\_2.fastq.gz -O /data/projects/lross_ssa/raw/$species/$sex/$SRR\_2.fastq.gz && gunzip /data/projects/lross_ssa/raw/$species/$sex/$SRR\_2.fastq.gz && /exports/software/fastqc/fastqc_v0.11.5/FastQC/fastqc --outdir /data/projects/lross_ssa/analyses/$species/fastqc /data/projects/lross_ssa/raw/$species/$sex/$SRR\_2.fastq
#	elif [ $SRRn -eq 12 ]
#	then
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${SRR:0:6}/${SRR:9:13}/$SRR/$SRR\_1.fastq.gz -O /data/projects/lross_ssa/raw/$species/$sex/$SRR\_1.fastq.gz && gunzip /data/projects/lross_ssa/raw/$species/$sex/$SRR\_1.fastq.gz && /exports/software/fastqc/fastqc_v0.11.5/FastQC/fastqc --outdir /data/projects/lross_ssa/analyses/$species/fastqc /data/projects/lross_ssa/raw/$species/$sex/$SRR\_1.fastq
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${SRR:0:6}/${SRR:9:13}/$SRR/$SRR\_2.fastq.gz -O /data/projects/lross_ssa/raw/$species/$sex/$SRR\_2.fastq.gz && gunzip /data/projects/lross_ssa/raw/$species/$sex/$SRR\_2.fastq.gz && /exports/software/fastqc/fastqc_v0.11.5/FastQC/fastqc --outdir /data/projects/lross_ssa/analyses/$species/fastqc /data/projects/lross_ssa/raw/$species/$sex/$SRR\_2.fastq
#	else
#echo 'Error: command line arguments incorrect'
#	fi
#elif [ $mode == 'single' ]
#	then
#	if [ $SRRn -eq 9 ]
#then
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${SRR:0:6}/$SRR/$SRR.fastq.gz -O /data/projects/lross_ssa/raw/$species/$sex/$SRR.fastq.gz && gunzip /data/projects/lross_ssa/raw/$species/$sex/$SRR.fastq.gz && /exports/software/fastqc/fastqc_v0.11.5/FastQC/fastqc --outdir /data/projects/lross_ssa/analyses/$species/fastqc /data/projects/lross_ssa/raw/$species/$sex/$SRR.fastq
#elif [ $SRRn -eq 10 ]
#then
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${SRR:0:6}/00${SRR:9:11}/$SRR/$SRR.fastq.gz -O /data/projects/lross_ssa/raw/$species/$sex/$SRR.fastq.gz && gunzip /data/projects/lross_ssa/raw/$species/$sex/$SRR.fastq.gz && /exports/software/fastqc/fastqc_v0.11.5/FastQC/fastqc --outdir /data/projects/lross_ssa/analyses/$species/fastqc /data/projects/lross_ssa/raw/$species/$sex/$SRR.fastq
#elif [ $SRRn -eq 11 ]
#then
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${SRR:0:6}/0${SRR:9:12}/$SRR/$SRR.fastq.gz -O /data/projects/lross_ssa/raw/$species/$sex/$SRR.fastq.gz && gunzip /data/projects/lross_ssa/raw/$species/$sex/$SRR.fastq.gz && /exports/software/fastqc/fastqc_v0.11.5/FastQC/fastqc --outdir /data/projects/lross_ssa/analyses/$species/fastqc /data/projects/lross_ssa/raw/$species/$sex/$SRR.fastq
#elif [ $SRRn -eq 12 ]
#then
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${SRR:0:6}/${SRR:9:13}/$SRR/$SRR.fastq.gz -O /data/projects/lross_ssa/raw/$species/$sex/$SRR.fastq.gz && gunzip /data/projects/lross_ssa/raw/$species/$sex/$SRR.fastq.gz && /exports/software/fastqc/fastqc_v0.11.5/FastQC/fastqc --outdir /data/projects/lross_ssa/analyses/$species/fastqc /data/projects/lross_ssa/raw/$species/$sex/$SRR.fastq
#else
#echo 'Error: command line arguments incorrect'
#fi
#fi

