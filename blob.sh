# BLOBTOOLS

#https://blobtools.readme.io/

# Activate blob tools virtenv
conda activate blobtools

SPECIES=$1
LAYOUT=$2

cd /projects/sykesj/analyses/$SPECIES

# Getting direcories for SRA libs
LIBS=$(for dir in $(ls /projects/sykesj/analyses/$SPECIES/kallisto/kal_results/kal_files/); do echo $dir; done | sed -e 's/^/-c /g' | sed 's/:/.cov/g' | paste -sd " " - )

# extracting coverage from kallisto's abundance TSV
for dir in kallisto/SRR*; do echo $dir; cut -f1,5 $dir/abundance.tsv | grep -v target > $dir.cov; done

#Ensure that the arguments trinity/*, blast/* and -c *.cov in the following line are correct.
# Creating a blobDB  
./blobtools create -i trinity/$LAYOUT_assembly_1k.fa -t blast/$SPECIES_blastn_$LAYOUT_sorted.out -o blobplot $LIBS

# Extracting a "view" table
/exports/software/blobtools/blobtools view -i blobplot.blobDB.json --rank all --hits

# Making blobplots (not really necessary)
#/exports/software/blobtools/blobtools plot -i blobplot.blobDB.json && rm -f *[0-9].png

mkdir blobtools
mv blob* blobtools/
cd blobtools

# deactivate blob tools virtenv
conda deactivate

#Open blobplot.blobDB.table.txt useing less and find the numbers next to superkingdom.t and phylum.t
#Use these numbers to edit the following six lines of code accordingly.

# Getting a distribution of kingdoms/phyla
# Kingdom
grep -v '^#' blobplot.blobDB.table.txt | cut -f9 | sort | uniq -c | less ### superkingdom.t.
# Phylum
grep -v '^#' blobplot.blobDB.table.txt | cut -f13 | sort | uniq -c | less ### phylum.t.

# Look at those that were annotated as Viruses 
#awk '$18=="Viruses"' blobplot.blobDB.table.txt | less
#awk '$13=="Streptophyta"' blobplot.blobDB.table.txt | less


# Filtering abundance before SLEUTH
awk '$10=="Viruses"' blobplot.blobDB.table.txt | cut -f1 > viruses.contig_ids.txt
awk '$14=="Streptophyta"' blobplot.blobDB.table.txt | cut -f1 > streptophyta.contig_ids.txt



abundance_files=()
mapfile -t abundance_files < <(for file in $(ls /projects/sykesj/analyses/$SPECIES/kallisto/SRR*/*.tsv); do readlink -f $file; done)
echo ${abundance_files[0]}

for i in "${abundance_files[@]}"
do
	LIB_SRR=$($i | sed 's/\/projects\/sykesj\/analyses\/testp\/kallisto\///g' | sed 's/\/abundance.tsv//g')
	grep -v -wFf viruses.contig_ids.txt $i | grep -v -wFf streptophyta.contig_ids.txt > /projects/sykesj/analyses/$SPECIES/kallisto/$LIB_SRR/abundance.filtered.tsv
done




#abundance_files=()
#mapfile -t abundance_files < <(for file in $(ls /projects/sykesj/analyses/$SPECIES/kallisto/SRR*/*.tsv); do readlink -f $file; done)
#echo ${abundance_files[0]}
#for i in "${abundance_files[@]}"
#do
#	LIB_SRR=$($i | sed 's/\/projects\/sykesj\/analyses\/testp\/kallisto\///g' | sed 's/\/abundance.tsv//g')
#	grep -v -wFf viruses.contig_ids.txt $i | grep -v -wFf streptophyta.contig_ids.txt > /projects/sykesj/analyses/$SPECIES/kallisto/$LIB_SRR/abundance.filtered.tsv
#done