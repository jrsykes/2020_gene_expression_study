import pandas as pd
import sys

SRAFile = sys.argv[1]
outDIR = sys.argv[2]

dat = pd.read_csv(SRAFile, header=None, dtype = str)

fail_list=[]

for index, row in dat.iterrows():
	species=row[0]
	file = '~/Documents/2020_gene_expression_study/multiqc_data/' + species + '/fastqc2/multiqc_data/multiqc_fastqc.txt'
	fastqc_file = pd.read_csv(file, dtype = 'str', header = 0, sep = '\t')
	for index, row in fastqc_file.iterrows():
		fail_counter = 0
		for i in row:
			if i == 'fail':
				fail_counter += 1
		if fail_counter > 1:
			fail_list.append(species + ',' + row['Sample'] + ',' + str(fail_counter))
			

df = pd.DataFrame(fail_list)

file_name = outDIR + 'QC_failed_libraries.txt'
df.to_csv(file_name, index = False)


