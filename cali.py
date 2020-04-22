import pandas as pd
import subprocess
import sys
import os
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

FinalOutfile_dir = str(sys.argv[2]) + '/out/'


##################################


#######################################################
# Create list of all contig blast IDs with no duplicates

print('\n Building list of blast IDs \n')

blast_id_list = []

for item in blast_list:
	df = pd.read_csv(item, sep='\t', usecols = [4], dtype = str)
	for index, row in df.iterrows():
		if row[0] not in blast_id_list:
			blast_id_list.append(row[0])

blast_id_list = list(dict.fromkeys(blast_id_list))

print('\n Blast ID list built \n')

####################################################################


# Create a file for each SRA library with trinity ID, blast ID and TPM or each contig
###############################

print('\n Compiling Trinity IDs, Blast IDs & TMP counts \n')


####################################################
# For SRA in SRA list
#Merge Trinity IDs, BLAST IDs and Kalisoto TPMs for each contig


blast_id_list = []

for index, row in dat.iterrows():
	SRR = row[1]
	command = 'ls ' + CaliWDout
	check = str(subprocess.check_output(command, shell=True))
	if SRR not in check:
		try:
			blast = CaliWDin + row[0] + '_blastn_PAIRED_sorted.out'
			blast_df = pd.read_csv(blast, sep='\t', header=None, usecols = [0, 4], dtype=str)
		except:
			blast = CaliWDin + row[0] + '_blastn_SINGLE_sorted.out'
			blast_df = pd.read_csv(blast, sep='\t', header=None,  usecols = [0, 4], dtype=str)
		
		outFile = FinalOutfile_dir + SRR + '_CaliOut.csv'
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
									
					if blast_id not in blast_id_list:
						blast_id_list.append(blast_id)
				except:
					pass





print('Complilation complete \n')
	
###############################################################################################################################################

###########################################################
# Create data frame with all contigs related to blast IDs and tpm. Then write to CSV

print('Producing final output file \n')


final_df = pd.DataFrame()

final_df['Blast_ID'] = blast_id_list

for index, row in dat.iterrows():
	SRR = row[1]
	sex = str(row[2])
	SexDeterm = str(row[4])
	species = str(row[0])
	tpm_list = []
	tpm_list.append(species)
	tpm_list.append(SexDeterm)

	cali_abundance_file = FinalOutfile_dir + SRR + '_CaliOut.csv'
	cali_abundance_df = pd.read_csv(cali_abundance_file, header=None, dtype={0: str, 1: str, 2: float})

	for i in blast_id_list:
		search_out = cali_abundance_df[cali_abundance_df[1].str.match(i)]
		tpm = search_out[2].sum()


		tpm_list.append(tpm)

	final_df[SRR] = tpm_list
	
out = FinalOutfile_dir + 'CaliOut.csv'
final_df.to_csv(out, index=False)

print('Program complete \n')
