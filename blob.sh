#!/bin/bash
#SBATCH --partition=short
#SBATCH --time=0-01:00:00
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --output=/home/sykesj/scripts/StdOut/R-%x.%j.out
#SBATCH --error=/home/sykesj/scripts/StdOut/R-%x.%j.err


SPECIES=$1
SRR=$2
SEX=$3
LAYOUT=$4

SPECIES=testp
SRR=SRR567165
SEX=female
LAYOUT=PAIRED

mkdir -p /projects/sykesj/analyses/$SPECIES/blobtools/$SRR
cd /projects/sykesj/analyses/$SPECIES/blobtools/$SRR

# Getting direcories for SRA libs
#LIBS=$(for dir in $(ls /projects/sykesj/analyses/$SPECIES/kallisto/kal_results/kal_files/); do echo $dir; done | sed -e 's/^/-c /g' | sed 's/:/.cov/g' | paste -sd " " - )

# extracting coverage from kallisto's abundance TSV
for dir in /projects/sykesj/analyses/$SPECIES/kallisto/$SRR; do echo $dir; cut -f1,5 $dir/abundance.tsv | grep -v target > $dir.cov; done

#Ensure that the arguments trinity/*, blast/* and -c *.cov in the following line are correct.
# Creating a blobDB  
/home/sykesj/software/blobtools-blobtools_v1.1.1/./blobtools create -i /projects/sykesj/analyses/$SPECIES/trinity/$LAYOUT"_assembly_1k.fa" -t /projects/sykesj/analyses/$SPECIES/blast/$SPECIES"_blastn_"$LAYOUT"_sorted.out" -o blobplot -c /projects/sykesj/analyses/$SPECIES/kallisto/$SRR".cov"

# Extracting a "view" table
/home/sykesj/software/blobtools-blobtools_v1.1.1/./blobtools view -i blobplot.blobDB.json --rank all --hits

# Making blobplots (not really necessary)
#/exports/software/blobtools/blobtools plot -i blobplot.blobDB.json && rm -f *[0-9].png


#Open blobplot.blobDB.table.txt useing less and find the numbers next to superkingdom.t and phylum.t
#Use these numbers to edit the following six lines of code accordingly.

# Getting a distribution of kingdoms/phyla
# Kingdom
#grep -v '^#' blobplot.blobDB.table.txt | cut -f9 | sort | uniq -c | less ### superkingdom.t.
# Phylum
#grep -v '^#' blobplot.blobDB.table.txt | cut -f13 | sort | uniq -c | less ### phylum.t.

# Look at those that were annotated as Viruses 
#awk '$18=="Viruses"' blobplot.blobDB.table.txt | less
#awk '$13=="Streptophyta"' blobplot.blobDB.table.txt | less


# Filtering abundance before SLEUTH
awk '$10=="Viruses"' blobplot.blobDB.table.txt | cut -f1 > viruses.contig_ids.txt
awk '$14=="Streptophyta"' blobplot.blobDB.table.txt | cut -f1 > streptophyta.contig_ids.txt


grep -v -wFf viruses.contig_ids.txt /projects/sykesj/analyses/$SPECIES/kallisto/$SRR/abundance.tsv | grep -v -wFf streptophyta.contig_ids.txt \
	> /projects/sykesj/analyses/$SPECIES/kallisto/$SRR/abundance.filtered.tsv

