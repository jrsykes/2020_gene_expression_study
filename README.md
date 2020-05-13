# Hoonah


### About hoonah 


This pipeline is made up of three distinct stages. 1. Quantifying transcriptome data with the hoonah bioinformatics peipline, 2. Cleaning and organising that data with cali.py and 3. Analysing that data with a Baysian phylogenetic mixed model using principal componenet analysis for dimensionality reeduction and then MCMCglmm.

### Input
1.A csv file, of any length, containg RNA-seq library acession numbers and library meta-data in the following format:
	species(genus_species),NCBI library acession number,condition 1,library layout,condition 2
	e.g.
	diachasma_alloeum,SRR2041626,female,PAIRED,haplodiploid 
2. An ultrametric clodogram, in nexus format, containg all species to be analysed. (Only necessary for MCMCglmm analysis)

### Output
1. 	Two multiQC files for each species (before and after trimming of reads). 					<WD>/analyses/<species>/fastqc & <WD>/analyses/<species>/fastqc2
2.	BUSCO results for the completeness of the de novo assembled transcriptome.					<WD>/analyses/<species>/busco
3.	Lists of all viral and streptophye reads removed form the transcriptomes. 					<WD>/analyses/<species>/blobtools/<SRR>/contig_ids.txt
		If you are curious what reads were found, these can be searched againsed 
		the blast output where you will find a blast acsession number. 							<WD>/analyses/<species>/blast/<species>_blastn_<layout>_sorted.out
		You can that search via NCBI. 
4. 	Files containing Kalisto tpm counts for all filtered reads from every input SRA library. 	<WD>/analyses/<species>/kallisto/<SRR>_abundance.filtered.tsv 
	At this point, you may want my to depart from this analysis as use other tools such as sleauth or Gfold to analyse this data, though these methods have problems with arbitrary definitions of biased gene expression and do not allow for within-gene expression varyance. These methods also treat isoforms as different genes.
5.	A single file from cali.py which contains metadata and transcript counts for all librarties in the analysis. These transcripts are named by a NCBI acession IDs and all isoforms have been pooled.
6.	Principal component plots showing the effect of two conditions on whole transcriptome gene expression. e.g. Sex determiantion system and sex.
7.	MCMCglmm output showing the relationship between the specfied conditions and whole transcriptome gene expression, controlled for non-independance of species.  


### Running hoonah 


1. Install all dependancies such as fastq-bump, fastQC, trimmomatic, trinity etc.. This may take a while.

2. Create a working directory and in that directory create two subdirectories names 'raw' and 'analyses' Also create an empty working directory by the same name on your scratch drive. e.g. /scratch/WD

3. Create a csv data file with the following format:

species(genus_species),NCBI library acession number,sex,library layout,condition
e.g.
diachasma_alloeum,SRR2041626,female,PAIRED,haplo 

This file should contain all libraries that you with to analyses

4. This pipeline was writen to run on with SLURM queing system and makes heavy use of the commands 'sbatch' and 'squeue'. If, for example, your cluser uses an SGE queing system where jobs are submited with qsub instead of sbatch, you will need to edit two lines in the hoonah.py file as follows, ensureing that a smapce is left between the end of the command on the second apostrophy:

squeue = 'qstat '
sbatch = 'qsub '


If this patch fails, you will have to edit all instances of 'squeue' and 'sbatch' in the hoonah.py file acordingly... not the end of the world!


5. Edit the beginning of all .sh scripts and cali.py to match your cluster. e.g.:

#!/bin/bash
#SBATCH --partition=medium
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --mem=100gb
#SBATCH --ntasks=6
#SBATCH --output=/home/<USER>/scripts/StdOut/R-%x.%j-download.out
#SBATCH --error=/home/<USER>/scripts/StdOut/R-%x.%j-download.err

6. Run hoonah.py as follows: pyhton3 hoonah.py <PATH TO DATA FILE> <WORKING DIRECTORY> <PATH TO hoonah DIR> <USER NAME>

e.g.

'<pyhton3 hoonah.py ~/dat/data.csv ~/hoonah_WD/ ~/software/hoonah fryphilipj>'

7. Once the above analysis is complete, which may atke several weeks depending on your resources and the number of libraries that you include, it is time to run cali.py
	1. Create a working directory (call it what you like) and within it, two subdirectories names 'in' and 'out'.
	2. In the 'in' subdirectory place all of the blast output files and filtered kalisto acundance files created by the hoonah.py pipeline. Symbolic liks work fine.
		These files can be found in: <hoonah WD>/analyses/<species>/kallisto/
									 <hoonsh WD>/analyses/<species>/blast/
		Where both paired end and single end libraries have been given to hoonah.py for a single species, be carful to give cali.py the correct blast file.
		In the cases, look in <hoonah WD>/analyses/<species>/kallisto/ where you will see either 'mapped_to_SINGLE_idx' or 'mapped_to_SINGLE_idx'
			This is your indication as weather the paired or single end libraries performed better in trascriptome assembly and were thus mapped to and so which blast file should be given to cali.py
	3. Edit the beginig of cali.py to best suit your cluster.
		e.g.:
			#!miniconda3/bin/python3
			#SBATCH --partition=long
			#SBATCH --time=UNLIMITED
			#SBATCH --nodes=1
			#SBATCH --mem=200gb
			#SBATCH --ntasks=2
			#SBATCH --output=/home/sykesj/scripts/StdOut/R-%x.%j-cali.out
			#SBATCH --error=/home/sykesj/scripts/StdOut/R-%x.%j-cali.err

	4. Run cali.py with the following command:	'<sbatch cali.py "PATH TO DATA FILE"/data.csv "path to"/CaliWD "n CPUs">'
			**IMPORTANT. _cali.py is writen to use multiprocessing. Do not give cali.py more CPUs than there are libraries in you data.csv file. It won't run._**
		Depending you resources and number of libraries, cali.py will take several hours/days to run

	5. Run the MCMCglmm R script. This will need to be done interactivly as you will need to adjust graphs, review statistic, select number of principal componenets, adjust priors for MCMCglmm etc. throughout. This should be relativly self explanitory.

8. Enjoy! Any questions or bug reports? Please contact me as jamie_r_sykes@outlook.com 


FAQs

Q. 	What to do if the pipline is terminated before it is finished?
A. 	Find out which species it was working on as the time, delete all data for these species from <WD/analyses>, <WD/raw> and </scratch/WD> and relaunch. 
	This will work for hoonah.py and cali.py. Both programs will pick up where they left off. 


