import pandas as pd
import subprocess
import sys



####################################
# Parsing input files


dat = pd.read_csv("SRA_list_refined.csv", header=None)



blast_list = []
abundace_list= []



for filename in sys.argv[1:]:
	if 'blastn' in filename:
		blast_list.append(filename)
	elif 'abundance.filtered.tsv' in filename:
		abundace_list.append(filename)

Outfile = sys.argv[2]

##################################


#######################################################
# Create list of all contig blast IDs with no duplicates

#print('\n Building list of blast IDs \n')#

#blast_id_list = []#

#for item in blast_list:
#	df = pd.read_csv(item, sep='\t')
#	for index, row in df.iterrows():
#		if row[4] not in blast_id_list:
#			blast_id_list.append(row[4])#

#blast_id_list = list(dict.fromkeys(blast_id_list))#

#print('\n Blast ID list built \n')

####################################################################33


# Create a file for each SRA library with trinity ID, blast ID and TPM or each contig
###############################

print('\n Compiling Trinity IDs, Blast IDs & TMP counts \n')


####################################################
# For SRA in SRA list
#Merge Trinity IDs, BLAST IDs and Kalisoto TPMs for each contig


blast_id_list = []

for index, row in dat.iterrows():
	SRR = row[1]
	check = str(subprocess.check_output('ls', shell=True))
	if SRR not in check:
		try:
			blast = 'dat/' + row[0] + '_blastn_PAIRED_sorted.out'
			blast_df = pd.read_csv(blast, sep='\t', header=None, low_memory=False)
		except:
			blast = 'dat/' + row[0] + '_blastn_SINGLE_sorted.out'
			blast_df = pd.read_csv(blast, sep='\t', header=None, low_memory=False)
		
		outFile = SRR + '_CaliOut.csv'
		abundance_file = 'dat/' + SRR + '_abundance.filtered.tsv'
		abundance_df = pd.read_csv(abundance_file, sep='\t', header=0)
		
		for index, row in abundance_df.iterrows():
				trinity_id = row[0]
				tpm = row[4]
				try:
					search_out = list(blast_df[blast_df[0].str.match(trinity_id)].iloc[0])
					blast_id = search_out[4]
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
	condition = str(row[2]) + str(row[4])
	species = str(row[0])
	tpm_list = []
	tpm_list.append(species)
	tpm_list.append(condition)

	cali_abundance_file = SRR + '_CaliOut.csv'
	cali_abundance_df = pd.read_csv(cali_abundance_file, header=None)

	for i in blast_id_list:
		search_out = cali_abundance_df[cali_abundance_df[1].str.match(i)]
		tpm = search_out[2].sum()


		tpm_list.append(tpm)

	final_df[SRR] = tpm_list
	

final_df.to_csv(Outfile, index=False)

print('Program complete \n')
