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

STRING=$(grep superkingdom.t. blobplot.blobDB.table.txt)
SUB='superkingdom.t.'

for word in $STRING; do
	if [[ "$word" == *"$SUB"* ]];
	then
	KINGDOM=$word
	KINGDOM_NUM=$(sed 's/[^0-9]*//g' <<< $word)
	fi
done

STRING=$(grep phylum.t. blobplot.blobDB.table.txt)
SUB='phylum.t.'

for word in $STRING; do
	if [[ "$word" == *"$SUB"* ]];
	then
	PHYLUM=$word
	PHYLUM_NUM=$(sed 's/[^0-9]*//g' <<< $word)
	fi
done

# Getting a distribution of kingdoms/phyla
# Kingdom
grep -v '^#' blobplot.blobDB.table.txt | cut -f"$KINGDOM" | sort | uniq -c | less ### superkingdom.t.
# Phylum
grep -v '^#' blobplot.blobDB.table.txt | cut -f"$PHYLUM" | sort | uniq -c | less ### phylum.t.

########################################################################################3

cat t.awk
NR==1 {
    for (i=1; i<=NF; i++) {
        ix[$i] = i
    }
}
NR>1 {
    print $ix[c1], $ix[c2]
}
awk -f t.awk "$KINGDOM" "$PHYLUM"  blobplot.blobDB.table.txt

########################################################################################

# Look at those that were annotated as Viruses 
awk '$6=="Viruses"' blobplot.blobDB.table.txt | less
awk '$13=="Streptophyta"' blobplot.blobDB.table.txt | less


# Filtering abundance before SLEUTH
awk '$(ech=="Viruses"' blobplot.blobDB.table.txt | cut -f1 > viruses.contig_ids.txt
awk '$PHYLUM=="Streptophyta"' blobplot.blobDB.table.txt | cut -f1 > streptophyta.contig_ids.txt


grep -v -wFf viruses.contig_ids.txt /projects/sykesj/analyses/"$SPECIES"/kallisto/"$SRR"/abundance.tsv | grep -v -wFf streptophyta.contig_ids.txt \
	> /projects/sykesj/analyses/"$SPECIES"/kallisto/"$SRR"/abundance.filtered.tsv

grep -v -wFf viruses.contig_ids.txt /projects/sykesj/analyses/"$SPECIES"/kallisto/"$SRR"/abundance.tsv > /projects/sykesj/analyses/"$SPECIES"/kallisto/"$SRR"/abundance.filtered.tsv