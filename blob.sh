#!/bin/bash
#SBATCH --partition=short
#SBATCH --time=0-01:00:00
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --ntasks=1
#SBATCH --output=/projects/sykesj/StdOut/R-%x.%j-blob.out
#SBATCH --error=/projects/sykesj/StdOut/R-%x.%j-blob.err


SPECIES=$1
SRR=$2
LAYOUT=$3
WD=$4


mkdir -p "$WD"/analyses/$SPECIES/blobtools/$SRR
cd "$WD"/analyses/$SPECIES/blobtools/$SRR

# Getting direcories for SRA libs
#LIBS=$(for dir in $(ls "$WD"/analyses/$SPECIES/kallisto/kal_results/kal_files/); do echo $dir; done | sed -e 's/^/-c /g' | sed 's/:/.cov/g' | paste -sd " " - )

# extracting coverage from kallisto's abundance TSV
for dir in "$WD"/analyses/$SPECIES/kallisto/"$SRR"; do echo $dir; cut -f1,5 $dir/abundance.tsv | grep -v target > $dir.cov; done

#Ensure that the arguments trinity/*, blast/* and -c *.cov in the following line are correct.
# Creating a blobDB  
/home/sykesj/software/blobtools-blobtools_v1.1.1/./blobtools create -i "$WD"/analyses/"$SPECIES"/trinity/"$SPECIES"_"$LAYOUT"_assembly_1k.fa \
	-t "$WD"/analyses/"$SPECIES"/blast/"$SPECIES"_blastn_"$LAYOUT"_sorted.out -o blobplot -c "$WD"/analyses/"$SPECIES"/kallisto/"$SRR".cov

# Extracting a "view" table
/home/sykesj/software/blobtools-blobtools_v1.1.1/./blobtools view -i blobplot.blobDB.json --rank all --hits

# Custom Python script to extract virus and streptophyte contigs

python3  /home/sykesj/scripts/2020_gene_expression_study/BlobFilter.py blobplot.blobDB.table.txt contig_ids.txt

# Filter out contaminant contigs

grep -v -wFf contig_ids.txt "$WD"/analyses/"$SPECIES"/kallisto/"$SRR"/abundance.tsv > "$WD"/analyses/"$SPECIES"/kallisto/"$SRR"_abundance.filtered.tsv