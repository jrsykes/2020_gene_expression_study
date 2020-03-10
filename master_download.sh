#This script is will download the SRA libraries, run fastqc on those libraries and create the directory infrastructure for the rest of the analysis.
#To run it, first amend the relevant lines beginning with qsub. Ensure that only two qsub lines are unhashed per run.
#Once this is complete run the script as follows: 'bash master_download.sh <species_name>' . 
#Once these libraries have downloaded, re-run the scrip, unhashing the next two qsub lines. 
#Repeat this process until all libraries for this this have been downloaded.
#IMPORTANT. Only download two libraries at a time otherwise you will slow the operation of the cluster for other users.
#IMPORTNT. Ensure that species name is all in lower case and the space between the genus and the species is replaced with an underscore.
#Once downloading is complete for all libraries of a given species, run master_trim.sh

species=$1

if [ -z "$species" ]
then
echo 'Species name not present'
exit 
fi

mkdir /home/sykesj/raw/$species
mkdir /home/sykesj/raw/$species/male
mkdir /home/sykesj/raw/$species/female

mkdir /home/sykesj/analyses/$species
mkdir /home/sykesj/analyses/$species/fastqc
mkdir /home/sykesj/analyses/$species/fastqc2

mkdir /home/sykesj/analyses/$species/trinity
mkdir /home/sykesj/analyses/$species/busco
mkdir /home/sykesj/analyses/$species/kallisto

mkdir /home/sykesj/analyses/$species/trimmomatic
mkdir /home/sykesj/analyses/$species/trimmomatic/male
mkdir /home/sykesj/analyses/$species/trimmomatic/female


## one = species, two = SRR, three = sex, four = paired or single end mode

script=/home/sykesj/scripts/download_qc.sh

#sbatch $script veroa_destructor SRR5377264	female PAIRED
#sbatch $script veroa_destructor SRR5377265	male PAIRED
sbatch $script veroa_destructor SRR5377267	female PAIRED
sbatch $script veroa_destructor SRR5377268	female PAIRED
sleep 2h
sbatch $script veroa_destructor SRR5760818	male SINGLE
sbatch $script veroa_destructor SRR5760828	male SINGLE
sleep 2h
sbatch $script veroa_destructor SRR5760838	male SINGLE
sbatch $script veroa_destructor SRR5760848	male SINGLE
sleep 2h
sbatch $script veroa_destructor SRR8864012	female PAIRED
sbatch $script veroa_destructor SRR8100123	female PAIRED
sleep 2h
sbatch $script veroa_destructor SRR8100124	female PAIRED
sbatch $script veroa_destructor SRR8100122	female PAIRED








