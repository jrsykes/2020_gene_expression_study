import pandas as pd
import os
import time
import subprocess
import sys
import re

class bcolors:
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    OKWARNING = '\033[93m'

#################################
squeue = 'squeue '
sbatch = 'sbatch '
################################

dat_file = sys.argv[1]
WD = sys.argv[2]
scripts = sys.argv[3]
user = sys.argv[4]


dat = pd.read_csv(dat_file, header=None)

species_list = []
for index, row in dat.iterrows():
	species_list.append(row[0])

for i in (list(dict.fromkeys(species_list))):
	SpeciesCheck = str(subprocess.check_output('ls ' +  WD + '/raw/', shell=True))
	if i not in SpeciesCheck:
		species = i


		print(f"{bcolors.OKBLUE}")
		print(f"{bcolors.OKGREEN}Download SRR files, QC and trim")
		print(f"{bcolors.OKBLUE}")


		for index, row in dat.iterrows():
			try:
				if row[0] == species:
					SRR = row[1]
					sex = row[2]
					layout = row[3]
					command = sbatch + scripts + '/download_trim_qc.sh ' + species + ' ' + SRR + ' ' + sex + ' ' + layout + ' ' + WD
					subprocess.Popen([command], shell=True)
					time.sleep(5)
					check = int(subprocess.check_output(squeue + ' --user=' + user + ' | grep download | wc -l', shell=True))
					while check > 1:
						time.sleep(10)
						check = int(subprocess.check_output(squeue + ' --user=' + user + ' | grep download | wc -l', shell=True))
			except:
				pass
		check2 = str(subprocess.check_output(squeue + ' --user=' + user, shell=True))
		while 'download' in check2:
						time.sleep(10)
						check2 = str(subprocess.check_output(squeue + ' --user=' + user, shell=True))


		print(f"{bcolors.OKBLUE}")
		print(f"{bcolors.OKGREEN}Runing Trinity, BUSCO and Blast")
		print(f"{bcolors.OKBLUE}")

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

		check3 = int(subprocess.check_output(squeue + ' --user=' + user + ' | grep trinity | wc -l', shell=True))
		while check3 > 2:
			time.sleep(10)
			check3 = int(subprocess.check_output(squeue + ' --user=' + user + ' | grep trinity | wc -l', shell=True))


		if df_paired.empty == False:
			command = sbatch + scripts + '/trinity_busco_blast.sh ' + species + ' PAIRED ' + WD
			subprocess.Popen([command], shell=True)

##############################################################################################################################
# Trinity single 

		check4 = int(subprocess.check_output(squeue + ' --user=' + user + ' | grep trinity | wc -l', shell=True))
		while check4 > 2:
			time.sleep(10)
			check4 = int(subprocess.check_output(squeue + ' --user=' + user + ' | grep trinity | wc -l', shell=True))


		if df_single.empty == False:
			command = sbatch + scripts + '/trinity_busco_blast.sh ' + species + ' SINGLE ' + WD
			subprocess.Popen([command], shell=True)

##############################################################################################################################
# Pause pipeline for both transcriptomes to be assembled and compared if both paired and single end layout libraries are used
		
		time.sleep(30)
		print(f"{bcolors.OKWARNING}Paused for completion of transcriptome assembly(s)")
		print(f"{bcolors.OKBLUE}")
		trinity_check = 'yes'
		while trinity_check == 'yes':
			
			command = squeue + '--user=' + user + ' -h | awk \'{print $1}\''
			check5 = str(subprocess.check_output(command, shell=True))
			result = re.sub('[^0-9]',' ', check5)
			for i in result.split():
				try:
					file = '/projects/sykesj/StdOut/R-%x.' + i + '-Trinity.out'
					with open(file, 'r') as f:
						if species not in f.read():
							trinity_check = 'no'
				except:
					pass
			
			time.sleep(1800)

##########################################################

		print(f"{bcolors.OKBLUE}")
		print(f"{bcolors.OKGREEN}Maping SRA libraries to de novo transcriptome")
		print(f"{bcolors.OKBLUE}")


		for index, row in dat.iterrows():
			try:
				if row[0] == species:
					SRR = row[1]
					sex = row[2]
					layout = row[3]		

					command = sbatch + scripts + '/map.sh ' + species + ' ' + SRR + ' ' + sex + ' ' + layout + ' ' + DUEL_LAYOUT + ' ' + WD
					subprocess.Popen([command], shell=True)
			except:
				pass

		time.sleep(5)

		check6 = str(subprocess.check_output(squeue + ' --user=' + user, shell=True))
		while 'map' in check6:
						time.sleep(10)
						check6 = str(subprocess.check_output(squeue + ' --user=' + user, shell=True))


		print(f"{bcolors.OKBLUE}")
		print(f"{bcolors.OKGREEN}Filtering contaminant sequences")
		print(f"{bcolors.OKBLUE}")



		for index, row in dat.iterrows():
			try:
				if row[0] == species:
					SRR = row[1]
					sex = row[2]
					layout = row[3]		

					command = sbatch + scripts + '/blob.sh ' + species + ' ' + SRR + ' ' + layout + ' ' + WD
					subprocess.Popen([command], shell=True)
			except:
				pass




		print(f"{bcolors.OKBLUE}")
		print(f"{bcolors.OKGREEN}Pipeline complete for " + species)
		print(f"{bcolors.OKBLUE}")
