#!miniconda3/bin/python3
#SBATCH --partition=long
#SBATCH --time=5-00:00:00
#SBATCH --nodes=1
#SBATCH --mem=300gb
#SBATCH --ntasks=25
#SBATCH --output=/home/sykesj/scripts/StdOut/R-%x.%j.out
#SBATCH --error=/home/sykesj/scripts/StdOut/R-%x.%j.err

n_processes = 25


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
out_file = CaliWDout + 'CaliOut.csv'


def blast_id_list_builder():
	blast_id_list = []
	for item in blast_list:
		df = pd.read_csv(item, sep='\t', usecols = [4], dtype = str)
		for index, row in df.iterrows():
			if row[0] not in blast_id_list:
				blast_id_list.append(row[0])
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



blast_id_list = [item for sublist in blast_id_list for item in sublist]
dot_lst = ['.', '.', '.']
for i in dot_lst:
	blast_id_list.insert(1, i)

blast_id_df = pd.DataFrame(blast_id_list[1:])


print('Blast ID list built \n')

####################################################################


# Create a file for each SRA library with trinity ID, blast ID and TPM or each contig
###############################





####################################################
# For SRA in SRA list
#Merge Trinity IDs, BLAST IDs and Kalisoto TPMs for each contig



def compiler(chunk):
	"""Merge Trinity IDs, BLAST IDs and Kalisoto TPMs for each contig"""
	for index, row in chunk.iterrows():
		SRR = row[1]
		command = 'ls ' + CaliWDout
		check = str(subprocess.check_output(command, shell=True))
		if SRR not in check:
			outFile = CaliWDout + SRR + '_CaliOut.csv'
			try:
				blast = CaliWDin + row[0] + '_blastn_PAIRED_sorted.out'
				blast_df = pd.read_csv(blast, sep='\t', header=None, usecols = [0, 4], dtype=str)
			except:
				blast = CaliWDin + row[0] + '_blastn_SINGLE_sorted.out'
				blast_df = pd.read_csv(blast, sep='\t', header=None,  usecols = [0, 4], dtype=str)
			abundance_file = CaliWDin + SRR + '_abundance.filtered.tsv'
			abundance_df = pd.read_csv(abundance_file, sep='\t', usecols = ['target_id', 'tpm'], header=0, dtype={"target_id": str, "tpm": float})
			for index, row in abundance_df.iterrows():
					trinity_id = row[0]
					tpm = row[1]
					try:
						search_out = list(blast_df[blast_df[0].str.match(trinity_id)].iloc[0])
						blast_id = search_out[1]
						with open(outFile, 'a+') as f:
							f.write(str(trinity_id) + ',' + str(blast_id) + ',' + str(tpm) + '\n')
					except:
						pass


   





###############################################################################################################################################

###########################################################
# Create data frame with all contigs related to blast IDs and tpm. Then write to CSV


def ID_tpm_combiner(dat):
	
	out_df = pd.DataFrame()

	for index, row in dat.iterrows():
	#for index, row in chunk.iterrows():
		species = row[0]
		SRR = row[1]
		sex = row[2]
		SexDeterm = row[4]
		tpm_list = []
		#tpm_list.append(SRR)
		tpm_list.append(species)
		tpm_list.append(SexDeterm)
		tpm_list.append(sex)
		cali_abundance_file = CaliWDout + SRR + '_CaliOut.csv'
		cali_abundance_df = pd.read_csv(cali_abundance_file, header=None, dtype={0: str, 1: str, 2: float})
		for i in blast_id_list[4:]:
			search_out = cali_abundance_df[cali_abundance_df[1].str.match(str(i))]
			tpm = search_out[2].sum()
			tpm_list.append(tpm)
		
		out_df[SRR] = tpm_list
	
	return out_df.values.tolist()


def collect_results(result):
    """Uses apply_async's callback to setup up a separate Queue for each process"""
    results.extend(result)	
	



chunk_size = int(dat.shape[0]/n_processes)
chunks = [dat.iloc[dat.index[i:i + chunk_size]] for i in range(0, dat.shape[0], chunk_size)]

print('Compiling Trinity IDs, Blast IDs & TMP counts \n')

pool = mp.Pool(processes=n_processes)
result = pool.map(compiler, chunks)
pool.close() 

print('Complilation complete \n')




print('Producing final output file \n')

#pool = mp.Pool(processes=n_processes)
#result = pool.map(ID_tpm_combiner, chunks)
#pool.close()  

################################################################################################33

results = []

if __name__ == "__main__":


	# Repeats the compute intensive operation on 10 data frames concurrently
	pool = mp.Pool(processes=n_processes)
	for i in range(dat.shape[1]): 
		pool.apply_async(ID_tpm_combiner, args=(dat, ), callback=collect_results)
	pool.close()
	pool.join()

	# Converts list of lists to a data frame
	out_df = pd.DataFrame(results)
	
	final_df = pd.concat([blast_id_df, out_df], axis = 1, sort = False)

	final_df.to_csv(out_file, index=False)

print('Program complete \n')