import pandas as pd
import os
import time
import subprocess

species = input("species_name: ")

dat = pd.read_csv("/home/sykesj/dat/SRA_list_refined.csv", header=None)

for index, row in dat.iterrows():
	try:
		if row[0] == species:
			species = row[0]
			SRR = row[1]
			sex = row[2]
			layout = row[3]

			#command = 'sbatch /home/sykesj/scripts/2020_gene_expression_study/new_download.sh ' + species + ' ' + SRR + ' ' + sex + ' ' + layout
			command = 'echo /home/sykesj/scripts/2020_gene_expression_study/new_download.sh ' + species + ' ' + SRR + ' ' + sex + ' ' + layout

			subprocess.Popen([command], shell=True)
			time.sleep(20)
			check = str(subprocess.check_output('squeue', shell=True))
			
			while 'sykesj' in check:
				time.sleep(60)
				check = str(subprocess.check_output('squeue', shell=True))
	except:
		pass
exit()

command = 'sbatch /home/sykesj/scripts/2020_gene_expression_study/trinity_busco_blast.sh ' + species + ' ' + layout
subprocess.Popen([command], shell=True)

time.sleep(20)
check = str(subprocess.check_output('squeue', shell=True))
			
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