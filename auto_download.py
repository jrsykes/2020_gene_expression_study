import pandas as pd
import os
import time
import subprocess

class bcolors:
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'


dat = pd.read_csv("/home/sykesj/dat/SRA_list_refined.csv", header=None)

species_list = []
for index, row in dat.iterrows():
	species_list.append(row[0])

for i in (list(dict.fromkeys(species_list))):
	SpeciesCheck = str(subprocess.check_output('ls /projects/sykesj/raw/', shell=True))
	if i not in SpeciesCheck:
		species = i

#		print(f"{bcolors.OKBLUE}######################")
#		print(f"{bcolors.OKGREEN}Clearing species data and preparing files")
#		print(f"{bcolors.OKBLUE}######################")

		#command = 'rm -rf /home/sykesj/busco_*.log ; rm -rf /projects/sykesj/raw/*' + species + '* ; rm -rf /projects/sykesj/analyses/*' + species + '* ; rm -rf /scratch/projects/sykesj/*' + species + '*'
		#subprocess.Popen([command], shell=True)


		print(f"{bcolors.OKBLUE}################################")
		print(f"{bcolors.OKGREEN}Download SRR files, QC and trim")
		print(f"{bcolors.OKBLUE}################################")


		for index, row in dat.iterrows():
			try:
				if row[0] == species:
					SRR = row[1]
					sex = row[2]
					layout = row[3]
					command = 'sbatch /home/sykesj/scripts/2020_gene_expression_study/new_download.sh ' + species + ' ' + SRR + ' ' + sex + ' ' + layout
					subprocess.Popen([command], shell=True)
					time.sleep(5)
					check = int(subprocess.check_output('squeue --user=sykesj | grep new_ | wc -l', shell=True))
					while check > 1:
						time.sleep(10)
						check = int(subprocess.check_output('squeue --user=sykesj | grep new_ | wc -l', shell=True))
			except:
				pass
		check2 = str(subprocess.check_output('squeue --user=sykesj', shell=True))
		while 'new_' in check2:
						time.sleep(10)
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

##############################################################################################################################
# Trinity paired

		check3 = int(subprocess.check_output('squeue --user=sykesj | grep trinity | wc -l', shell=True))
		while check3 > 2:
			time.sleep(10)
			check3 = int(subprocess.check_output('squeue --user=sykesj | grep trinity | wc -l', shell=True))


		if df_paired.empty == False:
			command = 'sbatch /home/sykesj/scripts/2020_gene_expression_study/trinity_busco_blast.sh ' + species + ' PAIRED'
			subprocess.Popen([command], shell=True)

##############################################################################################################################
# Trinity single 

		check4 = int(subprocess.check_output('squeue --user=sykesj | grep trinity | wc -l', shell=True))
		while check4 > 2:
			time.sleep(10)
			check4 = int(subprocess.check_output('squeue --user=sykesj | grep trinity | wc -l', shell=True))


		if df_single.empty == False:
			command = 'sbatch /home/sykesj/scripts/2020_gene_expression_study/trinity_busco_blast.sh ' + species + ' SINGLE'
			subprocess.Popen([command], shell=True)

##############################################################################################################################


		time.sleep(10)

		check5 = str(subprocess.check_output('squeue --user=sykesj', shell=True))
		while 'trinity' in check5:
						time.sleep(10)
						check5 = str(subprocess.check_output('squeue --user=sykesj', shell=True))


		print(f"{bcolors.OKBLUE}##############################################")
		print(f"{bcolors.OKGREEN}Maping SRA libraries to de novo transcriptome")
		print(f"{bcolors.OKBLUE}##############################################")

		if df_paired.empty == False and df_single.empty == False:
			DUEL_LAYOUT = 'YES'
		else:
			DUEL_LAYOUT = 'NO'


		for index, row in dat.iterrows():
			try:
				if row[0] == species:
					SRR = row[1]
					sex = row[2]
					layout = row[3]		

					command = 'sbatch /home/sykesj/scripts/2020_gene_expression_study/map.sh ' + species + ' ' + SRR + ' ' + sex + ' ' + layout + ' ' + DUEL_LAYOUT
					subprocess.Popen([command], shell=True)
			except:
				pass

		time.sleep(5)

		check6 = str(subprocess.check_output('squeue --user=sykesj', shell=True))
		while 'map' in check6:
						time.sleep(10)
						check6 = str(subprocess.check_output('squeue --user=sykesj', shell=True))


		print(f"{bcolors.OKBLUE}################################")
		print(f"{bcolors.OKGREEN}Filtering contaminant sequences")
		print(f"{bcolors.OKBLUE}################################")



		for index, row in dat.iterrows():
			try:
				if row[0] == species:
					SRR = row[1]
					sex = row[2]
					layout = row[3]		

					command = 'sbatch /home/sykesj/scripts/2020_gene_expression_study/blob.sh ' + species + ' ' + SRR + ' ' + layout
					subprocess.Popen([command], shell=True)
			except:
				pass




		print(f"{bcolors.OKBLUE}##################")
		print(f"{bcolors.OKGREEN}Pipeline complete")
		print(f"{bcolors.OKBLUE}##################")
