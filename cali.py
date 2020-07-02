#!miniconda3/bin/python3
#SBATCH --partition=long
#SBATCH --time=UNLIMITED
#SBATCH --nodes=1
#SBATCH --mem=100gb
#SBATCH --ntasks=40
#SBATCH --output=/projects/sykesj/StdOut/R-%x.%j-cali.out
#SBATCH --error=/projects/sykesj/StdOut/R-%x.%j-cali.err


import pandas as pd
import subprocess
import sys
import os
import multiprocessing as mp
####################################
# Parsing input files



SRAFile = sys.argv[1]


dat = pd.read_csv(SRAFile, header=None, dtype = str)



blast_list = []
abundace_list= []

CaliWDin = str(sys.argv[2]) + '/in/'
CaliWDout = str(sys.argv[2]) + '/out/'


n_processes = int(sys.argv[3])




#Create list of blast files

for i in os.listdir(CaliWDin):
	if 'blastn' in i:
		blast_list.append(CaliWDin + i)
	elif 'abundance.filtered.tsv' in i:
		abundace_list.append(CaliWDin + i)

##################################


#######################################################
# Create list of all contig blast IDs with no duplicates

print('\nBuilding list of blast IDs \n')


blast_out = CaliWDout + 'blast_id_list.csv'



def blast_id_list_builder():
	blast_id_list = []
	counter = 0
	for item in blast_list:
		df = pd.read_csv(item, sep='\t', usecols = [0,4,13], dtype={0: str, 4: str, 13: float}, header=None)
		trin_id_list = list(dict.fromkeys(df[0].tolist()))
		for i in trin_id_list:
			search_out = df[df[0] == i].sort_values(by=[13], ascending = False)
			if search_out.iloc[0,2] > 0:
				blast_id_list.append(search_out.iloc[0,1])
		counter += 1
		print('List', round(counter/len(blast_list)*100),'% complete')
	blast_id_list = list(dict.fromkeys(blast_id_list))
	blast_id_df = pd.DataFrame(blast_id_list)
	blast_id_df.to_csv(blast_out, index=False)
	


command = 'ls ' + CaliWDout
check = str(subprocess.check_output(command, shell=True))



if 'blast_id_list.csv' not in check:
	blast_id_list_builder()
	blast_id_list = pd.read_csv(blast_out, header=None, dtype = str).values.tolist()
else:
	blast_id_list = pd.read_csv(blast_out, header=None, dtype = str).values.tolist()

blast_id_list.remove(blast_id_list[0])
blast_id_list = [item for sublist in blast_id_list for item in sublist]
dot_lst = ['sex', 'SexDeterm', 'species']


for i in dot_lst:
	blast_id_list.insert(0, i)


print('\nBlast ID list built \n')

####################################################################


# Create a file for each SRA library with trinity ID, blast ID and TPM or each contig
###############################





####################################################
# For SRA in SRA list
#Merge Trinity IDs, BLAST IDs and Kalisoto TPMs for each contig



def compiler(chunk):
	"""Merge Trinity IDs, BLAST IDs and Kalisoto TPMs for each contig in a library"""
	for index, row in chunk.iterrows():
		SRR = row[1]
		command = 'ls ' + CaliWDout
		check = str(subprocess.check_output(command, shell=True))
		if SRR not in check:
			try:
				blast = CaliWDin + row[0] + '_blastn_PAIRED_sorted.out'
				blast_df = pd.read_csv(blast, sep='\t', header=None, usecols = [0, 4, 13], dtype={0: str, 4: str, 13: float})
			except:
				blast = CaliWDin + row[0] + '_blastn_SINGLE_sorted.out'
				blast_df = pd.read_csv(blast, sep='\t', header=None, usecols = [0, 4, 13], dtype={0: str, 4: str, 13: float})
			abundance_file = CaliWDin + SRR + '_abundance.filtered.tsv'
			abundance_df = pd.read_csv(abundance_file, sep='\t', usecols = ['target_id', 'tpm'], header=0, dtype={"target_id": str, "tpm": float})
						
			df_to_write = pd.DataFrame()
			for index, row in abundance_df.iterrows():
					trinity_id = row[0]
					tpm = row[1]
					try:
						search_out = blast_df[blast_df[0] == trinity_id].sort_values(by=[13], ascending = False)
						if search_out.iloc[0,2] > 0:
							blast_id = search_out.iloc[0,1]
							row = pd.Series([str(trinity_id), str(blast_id), str(tpm)])
							df_append = pd.DataFrame([row])
							df_to_write = pd.concat([df_to_write, df_append], ignore_index = True)
					except:
						pass
			out_file = CaliWDout + SRR + '_CaliOut.csv'
			df_to_write.to_csv(out_file, index=False)


   





###############################################################################################################################################

###########################################################
# Create data frame with all contigs related to blast IDs and tpm. Then write to CSV


def ID_tpm_combiner(chunk):
	
	out_df = pd.DataFrame()

	for index, row in chunk.iterrows():
		
		command = 'ls ' + CaliWDout
		check = str(subprocess.check_output(command, shell=True))
		search_str = row[1] + '_ID_tpm_combiner'
		if search_str not in check:

			species = row[0]
			SRR = row[1]
			sex = row[2]
			SexDeterm = row[4]
			tpm_list = []
			tpm_list.append(species)
			tpm_list.append(SexDeterm)
			tpm_list.append(sex)
			cali_abundance_file = CaliWDout + SRR + '_CaliOut.csv'
			cali_abundance_df = pd.read_csv(cali_abundance_file, header=None, dtype={0: str, 1: str, 2: float})
			for i in blast_id_list[3:]:
				search_out = cali_abundance_df[cali_abundance_df[1] == str(i)]
				tpm = search_out[2].sum()
				tpm_list.append(tpm)
			file = CaliWDout + SRR + '_ID_tpm_combiner'
			with open(file, 'w') as f:
				for i in tpm_list:
					f.write(str(i))
					f.write('\n')



chunk_size = int(dat.shape[0]/n_processes)
chunks = [dat.iloc[dat.index[i:i + chunk_size]] for i in range(0, dat.shape[0], chunk_size)]

print('Compiling Trinity IDs, Blast IDs & TMP counts \n')

pool = mp.Pool(processes=n_processes)
result = pool.map(compiler, chunks)
pool.close() 

print('Complilation complete \n')

################################################################################################33

print('Mapping IDs to respective TPMs \n')

pool = mp.Pool(processes=n_processes)
result = pool.map(ID_tpm_combiner, chunks)
pool.close() 

print('Mapping complete \n')

print('Producing final output file \n')

combiner_file_list = []
for i in os.listdir(CaliWDout):
	if 'ID_tpm_combiner' in i:
		combiner_file_list.append(CaliWDout + i)


final_out_df = pd.DataFrame()
final_out_df['SRR'] = blast_id_list



for index, row in dat.iterrows():
	file = CaliWDout + row[1] + '_ID_tpm_combiner'
	df = pd.read_csv(file, dtype = str, names=['x'])
	values = df.x.tolist()
	col_name = row[1]
	final_out_df[col_name] = values

final_out_file = CaliWDout + 'cali_out.csv'
final_out_df.to_csv(final_out_file, index=False)

print('Program complete \n')
