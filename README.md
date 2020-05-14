# Hoonah


### About hoonah 

Given a collection of RNA-seq libraries, this pipeline will analyse the effect of one or two conditions and their interaction on whole transcriptome data across a diversity of species simultaenously, while controlling for the non-independance of those species. For each species, at least one library of both the first condition (_e.g._ sex) must be supplied as well as a roughly even number of libraries of the second condition (_e.g._ sex determination system).
A combination of paired and single end layout libraries can be input. In this this case, a single end- and paired end- transcriptome will be assembled and all libraries for that species will be mapped to the more complete of those two transcriptomes.

A schematic of the pipeline can be seen at: https://github.com/jrsykes/hoonah/blob/master/pipeline_flow_chart.png

This pipeline was written to answer the question:
Can predicted patterns of sexually antagonistic selection be discerned from the transcriptomes of haplodiploid and diplodiploid arthropods? Paper not yet published.

This pipline can be run on oganisms other than arthropods, however the query data set for BUSCO will need to be changed acordingly and if your study concerns virual or streptophyte RNA, you will need to address this in the blob tools filtering step as these reads will e filtered out. 

This pipeline is made up of three distinct stages. 
1. Quantification of transcriptome data with the hoonah bioinformatics pipeline 
2. Cleaning and organising that data with cali.py
3. Analysing that data with a Bayesian phylogenetic mixed model in MCMCglmm, using principal component analysis for dimensionality reduction.

### Input
1. A csv file, of any length, containing RNA-seq library accession numbers and library meta-data in the following format:
	
		species(genus_species),NCBI library accession number,condition 1,library layout,condition 2
	
	_e.g._
	
		diachasma_alloeum,SRR2041626,female,PAIRED,haplodiploid 
2. An ultrametric cladogram, in nexus format, containing all species to be analysed. (Only necessary for MCMCglmm analysis)

### Output
1. 	Two multiQC files for each species (before and after trimming of reads). 	Found in:	"WD"/analyses/"species"/fastqc & "WD"/analyses/"species"/fastqc2
2.	BUSCO results for the completeness of the de novo assembled transcriptome. 	Found in:	"WD"/analyses/"species"/busco
3.	Lists of all viral and streptophyte reads removed form the transcriptomes.	Found in:	"WD"/analyses/"species"/blobtools/"SRR"/contig_ids.txt
		
	If you are curious what reads were found, these can be searched againsed the blast output where you will find a blast accsession number.	Found in:	"WD"/analyses/"species"/blast/"species"_blastn_"layout"_sorted.out
4. 	Files containing Kalisto tpm counts for all filtered reads from every input SRA library.	Found in:	"WD"/analyses/"species"/kallisto/"SRR"_abundance.filtered.tsv 
	_At this point, you may want my to depart from this analysis and use other tools such as sleuth or Gfold to analyse this data, though these methods have problems with arbitrary definitions of biased gene expression and do not allow for within-gene expression variance. These methods also treat isoforms as different genes._
5.	A single file from cali.py which contains metadata and transcript counts for all libraries in the analysis. These transcripts are named by NCBI accession IDs and all isoforms have been pooled.
6.	Principal component plots showing the effect of two conditions on whole transcriptome gene expression across all species. e.g. Sex determination system and sex.
7.	MCMCglmm output showing the relationship(or lack there of) between the specified conditions and whole transcriptome gene expression, controlled for non-independance of species.  


### Running hoonah 


1. Clone this repository and install all dependencies such as fastq-dump, fastQC, trimmomatic, trinity, kalisto, python3, blobtools, blast and R (I think that's everything).. This may take a while.

2. Create a working directory and, in that directory, create two subdirectories named 'raw' and 'analyses'. Also create an empty working directory by the same name on your scratch drive. _e.g._ /scratch/WD

3. Create a csv data file with the following format:

		species(genus_species),NCBI library accession number,sex,library layout,condition

_e.g._

		diachasma_alloeum,SRR2041626,female,PAIRED,haplo 

This file should contain all libraries that you with to analyses

4. This pipeline was written to run with a SLURM queuing system and makes heavy use of the commands 'sbatch' and 'squeue'. If, for example, your cluster uses an SGE queuing system, where jobs are submitted with qsub instead of sbatch, you will need to edit two lines in the hoonah.py file as follows, ensuring that a space is left between the end of the command on the second apostrophy:

		squeue = 'qstat '
		sbatch = 'qsub '


If this patch fails, you will have to edit all instances of 'squeue' and 'sbatch' in the hoonah.py file accordingly... not the end of the world!


5. Edit the beginning of all .sh scripts and cali.py to match your cluster. e.g.:

		#!/bin/bash
		#SBATCH --partition=medium
		#SBATCH --time=1-00:00:00
		#SBATCH --nodes=1
		#SBATCH --mem=100gb
		#SBATCH --ntasks=6
		#SBATCH --output=/home/<USER>/scripts/StdOut/R-%x.%j-download.out
		#SBATCH --error=/home/<USER>/scripts/StdOut/R-%x.%j-download.err

6. Run hoonah.py as follows:

    	pyhton3 hoonah.py DATA_FILE WORKING_DIRECTORY hoonah_DIR USER_NAME
		pyhton3 hoonah.py ~/dat/data.csv ~/hoonah_WD ~/software/hoonah fryphilipj

If you have large data set, you can run as many instances of hoonah as you like and they will not interfere with each other. Multiple instances will never download more that three libraries simultaneously or run more than two transcriptome assemblers (usually one) at a time. There are two good reasons for this. A) You won't use up all of the band width or compute resources on a shared cluster. _Being a good neighbour!_ and B) When analysing both paired and single end sequences for one species, this means that hoonah can choose the best transcriptome before progressing.

7. Once the above analysis is complete, which may take several weeks depending on your resources and the number of libraries that you include, it is time to run cali.py
	1. Create a working directory (call it what you like) and within it, two subdirectories named 'in' and 'out'.
	2. In the 'in' subdirectory place all of the blast output files and filtered kallisto abundance files created by the hoonah.py pipeline. Symbolic links work fine.
		These files can be found in: 

		"hoonah WD"/analyses/"species"/kallisto/
		"hoonah WD"/analyses/"species"/blast/
		
	Where both paired end and single end libraries have been given to hoonah.py for a single species, be careful to give cali.py the correct blast file.
	In this cases, look in hoonah_WD/analyses/species/kallisto/ where you will see either 'mapped_to_SINGLE_idx' or 'mapped_to_SINGLE_idx'
			This is your indication as weather the paired or single end libraries performed better in transcriptome assembly and were thus mapped to and so which blast file should be given to cali.py
	3. Edit the beginning of cali.py to best suit your cluster.
		_e.g.:_
			
			#!miniconda3/bin/python3
			#SBATCH --partition=long
			#SBATCH --time=UNLIMITED
			#SBATCH --nodes=1
			#SBATCH --mem=200gb
			#SBATCH --ntasks=2
			#SBATCH --output=/home/<USER>/scripts/StdOut/R-%x.%j-cali.out
			#SBATCH --error=/home/<USER>/scripts/StdOut/R-%x.%j-cali.err

	4. Run cali.py with the following command:	
				
		sbatch cali.py data.csv "path to"/CaliWD "n CPUs"
		
	cali.py makes use of the same data.csv file input to hoonah.py
	**IMPORTANT. _cali.py is written to use multiprocessing. Do not give cali.py more CPUs than there are libraries in you data.csv file. It won't run._**
	Depending you resources and number of libraries, cali.py will take several hours/days to run.

8. Run the MCMCglmm.R script. This will need to be done interactively as you will need to adjust graphs, review statistics, select number of principal components, adjust priors for MCMCglmm and more throughout. This should be relatively self explanatory.

9. Enjoy! Any questions or bug reports? Please contact me at jamie_r_sykes@outlook.com 


FAQs

Q. What to do if the pipline is terminated before it is finished?

A. Find out which species it was working on at the time, delete all data for those species from WD/analyses, WD/raw and /scratch/WD and relaunch. 
	This will work for hoonah.py and cali.py. Both programs will pick up where they left off. 


