

import os
import pandas as pd

species = pd.read_csv("/home/sykesj/raw/refined_species.csv")
fastq_dump = '/home/jamie/Documents/2020_gene_expression_study/sratoolkit.2.10.3-centos_linux64/bin/fastq-dump.2.10.3'



for index, row in species_dat.iterrows():
	species = row[0]
	SRR = row[4]
	sex = row[3]
	mode = row[5]
	command = fastq_dump + ' --outdir /home/jamie/Documents/2020_gene_expression_study/' + species + '/' + sex + '--defline-seq \'@$sn[_$rn]/$ri\' --split-files' + SRR
	print (command)
	#os.system(fastq_dump + '--outdir /home/jamie/Documents/2020_gene_expression_study/' + species + '/' + sex + '--defline-seq \'@$sn[_$rn]/$ri\' --split-files' + SRR)

#/exports/software/sratoolkit/sratoolkit.2.8.2-1-ubuntu64/bin/fastq-dump --outdir /data/projects/lross_ssa/raw/$species/$sex --defline-seq '@$sn[_$rn]/$ri' --split-files $SRR


