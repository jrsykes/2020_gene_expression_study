# BLOBTOOLS

#https://blobtools.readme.io/

# Activate blob tools virtenv
conda activate blobtools

SPECIES=$1
SRR=$2
SEX=$3
LAYOUT=$4

cd /projects/sykesj/analyses/$SPECIES

# extracting coverage from kallisto's abundance TSV
for dir in kallisto/SRR*; do echo $dir; cut -f1,5 $dir/abundance.tsv | grep -v target > $dir.cov; done

#Ensure that the arguments trinity/*, blast/* and -c *.cov in the following line are correct.
# Creating a blobDB  
/exports/software/blobtools/blobtools create -i trinity/$LAYOUT_assembly_1k.fa -t blast/$SPECIES_blastn_$LAYOUT_sorted.out -o blobplot -c kallisto/SRR1.cov -c kallisto/SRR3.cov -c kallisto/SRR4.cov -c kallisto/SRR2.cov #-c kallisto/SRR1239356.cov -c kallisto/SRR1239354.cov #-c kallisto/SRR4043739.cov -c kallisto/SRR4043738.cov -c kallisto/SRR1582619.cov -c kallisto/SRR1582618.cov -c kallisto/SRR1582616.cov -c kallisto/SRR1582617.cov #-c kallisto/SRR3934330.cov -c kallisto/SRR3934331.cov -c kallisto/SRR3934338.cov -c kallisto/SRR3934339.cov #-c kallisto/SRR1566026.cov -c kallisto/SRR1566025.cov 

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

a=SRR1
b=SRR3
c=SRR4
d=SRR2
#e=SRR1239356
#f=SRR1239354
#g=SRR4043739
#h=SRR4043738
#i=SRR1582619
#j=SRR1582618
#k=SRR1582616
#l=SRR1582617
#m=SRR3934330
#n=SRR3934331
#o=SRR3934338
#p=SRR3934339
#q=SRR1566026
#r=SRR1566025

grep -v -wFf viruses.contig_ids.txt ../kallisto/$a/abundance.tsv | grep -v -wFf streptophyta.contig_ids.txt > ../kallisto/$a/abundance.filtered.tsv
grep -v -wFf viruses.contig_ids.txt ../kallisto/$b/abundance.tsv | grep -v -wFf streptophyta.contig_ids.txt > ../kallisto/$b/abundance.filtered.tsv
grep -v -wFf viruses.contig_ids.txt ../kallisto/$c/abundance.tsv | grep -v -wFf streptophyta.contig_ids.txt > ../kallisto/$c/abundance.filtered.tsv
grep -v -wFf viruses.contig_ids.txt ../kallisto/$d/abundance.tsv | grep -v -wFf streptophyta.contig_ids.txt > ../kallisto/$d/abundance.filtered.tsv
#grep -v -wFf viruses.contig_ids.txt ../kallisto/$e/abundance.tsv | grep -v -wFf streptophyta.contig_ids.txt > ../kallisto/$e/abundance.filtered.tsv
#grep -v -wFf viruses.contig_ids.txt ../kallisto/$f/abundance.tsv | grep -v -wFf streptophyta.contig_ids.txt > ../kallisto/$f/abundance.filtered.tsv
#grep -v -wFf viruses.contig_ids.txt ../kallisto/$g/abundance.tsv | grep -v -wFf streptophyta.contig_ids.txt > ../kallisto/$g/abundance.filtered.tsv
#grep -v -wFf viruses.contig_ids.txt ../kallisto/$h/abundance.tsv | grep -v -wFf streptophyta.contig_ids.txt > ../kallisto/$h/abundance.filtered.tsv
#grep -v -wFf viruses.contig_ids.txt ../kallisto/$i/abundance.tsv | grep -v -wFf streptophyta.contig_ids.txt > ../kallisto/$i/abundance.filtered.tsv
#grep -v -wFf viruses.contig_ids.txt ../kallisto/$j/abundance.tsv | grep -v -wFf streptophyta.contig_ids.txt > ../kallisto/$j/abundance.filtered.tsv
#grep -v -wFf viruses.contig_ids.txt ../kallisto/$k/abundance.tsv | grep -v -wFf streptophyta.contig_ids.txt > ../kallisto/$k/abundance.filtered.tsv
#grep -v -wFf viruses.contig_ids.txt ../kallisto/$l/abundance.tsv | grep -v -wFf streptophyta.contig_ids.txt > ../kallisto/$l/abundance.filtered.tsv
#grep -v -wFf viruses.contig_ids.txt ../kallisto/$m/abundance.tsv | grep -v -wFf streptophyta.contig_ids.txt > ../kallisto/$m/abundance.filtered.tsv
#grep -v -wFf viruses.contig_ids.txt ../kallisto/$n/abundance.tsv | grep -v -wFf streptophyta.contig_ids.txt > ../kallisto/$n/abundance.filtered.tsv
#grep -v -wFf viruses.contig_ids.txt ../kallisto/$o/abundance.tsv | grep -v -wFf streptophyta.contig_ids.txt > ../kallisto/$o/abundance.filtered.tsv
#grep -v -wFf viruses.contig_ids.txt ../kallisto/$p/abundance.tsv | grep -v -wFf streptophyta.contig_ids.txt > ../kallisto/$p/abundance.filtered.tsv
#grep -v -wFf viruses.contig_ids.txt ../kallisto/$q/abundance.tsv | grep -v -wFf streptophyta.contig_ids.txt > ../kallisto/$q/abundance.filtered.tsv
#grep -v -wFf viruses.contig_ids.txt ../kallisto/$r/abundance.tsv | grep -v -wFf streptophyta.contig_ids.txt > ../kallisto/$r/abundance.filtered.tsv



