import pandas as pd
import os
import time
import subprocess

species = input("species_name: ")
#species = 'amblyomma_americanum'
dat = pd.read_csv("/home/sykesj/dat/SRA_list_refined.csv", header=None)
#dat = pd.read_csv("/home/jamie/Documents/2020_gene_expression_study/SRA_list_refined.csv", header=None)


###################
#Clear species data
###################
#command = 'rm -rf /projects/sykesj/raw/' + species + '; rm -rf /projects/sykesj/analyses/' + species + '; rm -rf /scratch/projects/sykesj/*' + species + '*'
#subprocess.Popen([command], shell=True)
####################

for index, row in dat.iterrows():
	try:
		if row[0] == species:
			species = row[0]
			SRR = row[1]
			sex = row[2]
			layout = row[3]

			command = 'sbatch /home/sykesj/scripts/2020_gene_expression_study/new_download.sh ' + species + ' ' + SRR + ' ' + sex + ' ' + layout
			#subprocess.Popen([command], shell=True)

			time.sleep(20)
			check = str(subprocess.check_output('squeue', shell=True))
			
			while 'sykesj' in check:
				time.sleep(60)
				check = str(subprocess.check_output('squeue', shell=True))
	except:
		pass

#########################################################################################################


df_paired = pd.DataFrame()
df_single = pd.DataFrame()


for index, row in dat.iterrows():
	if row[0] == species:
		if row[3] == 'PAIRED':
			df_paired = df_paired.append(row[0:4], ignore_index=True)
		if row[3] == 'SINGLE':
			df_single = df_single.append(row[0:4], ignore_index=True)


if df_paired.empty == False:
	command = 'sbatch /home/sykesj/scripts/2020_gene_expression_study/trinity_busco_blast.sh ' + species + ' PAIRED'
	subprocess.Popen([command], shell=True)

if df_single.empty == False:
	command = 'sbatch /home/sykesj/scripts/2020_gene_expression_study/trinity_busco_blast.sh ' + species + ' SINGLE'
	subprocess.Popen([command], shell=True)


time.sleep(20)
check = str(subprocess.check_output('squeue', shell=True))

exit()
		
while 'sykesj' in check:
	time.sleep(60)
	check = str(subprocess.check_output('squeue', shell=True))


for index, row in dat.iterrows():
	try:
		if row[0] == species:
			species = row[0]
			SRR = row[1]
			sex = row[2]
			layout = row[3]

			command = 'sbatch /home/sykesj/scripts/2020_gene_expression_study/map.sh ' + species + ' ' + SRR + ' ' + sex + ' ' + layout
			subprocess.Popen([command], shell=True)
	except:
		pass