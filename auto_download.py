import pandas as pd
import os
import time
import subprocess

species = input("species_name: ")

class bcolors:
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'

print(f"{bcolors.OKBLUE}######################")
print(f"{bcolors.OKGREEN}Clearing species data")
print(f"{bcolors.OKBLUE}######################")

command = 'rm -rf /projects/sykesj/raw/*' + species + '*; rm -rf /projects/sykesj/analyses/*' + species + '*; rm -rf /scratch/projects/sykesj/*' + species + '*'
subprocess.Popen([command], shell=True)


dat = pd.read_csv("/home/sykesj/dat/SRA_list_refined.csv", header=None)


print(f"{bcolors.OKBLUE}################################")
print(f"{bcolors.OKGREEN}Download SRR files, QC and trim")
print(f"{bcolors.OKBLUE}################################")


for index, row in dat.iterrows():
	try:
		if row[0] == species:
			species = row[0]
			SRR = row[1]
			sex = row[2]
			layout = row[3]

			command = 'sbatch /home/sykesj/scripts/2020_gene_expression_study/new_download.sh ' + species + ' ' + SRR + ' ' + sex + ' ' + layout
			subprocess.Popen([command], shell=True)

			time.sleep(20)
			check = int(subprocess.check_output('squeue --user=sykesj | wc -l', shell=True))

			while check > 2:
				time.sleep(20)
				check = int(subprocess.check_output('squeue --user=sykesj | wc -l', shell=True))
	except:
		pass


check2 = str(subprocess.check_output('squeue --user=sykesj', shell=True))
while 'new_' in check2:
				time.sleep(20)
				check2 = str(subprocess.check_output('squeue --user=sykesj', shell=True))


print(f"{bcolors.OKBLUE}################################")
print(f"{bcolors.OKGREEN}Runing Trinity, BUSCO and Blast")
print(f"{bcolors.OKBLUE}################################")

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

check2 = str(subprocess.check_output('squeue --user=sykesj', shell=True))
while 'trinity' in check2:
				time.sleep(20)
				check2 = str(subprocess.check_output('squeue --user=sykesj', shell=True))


print(f"{bcolors.OKBLUE}##############################################")
print(f"{bcolors.OKGREEN}Maping SRA libraries to de novo transcriptome")
print(f"{bcolors.OKBLUE}##############################################")


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


print(f"{bcolors.OKBLUE}##################")
print(f"{bcolors.OKGREEN}Pipeline complete")
print(f"{bcolors.OKBLUE}##################")