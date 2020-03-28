import pandas as pd
import re
import sys

inFile = sys.argv[1]
outFile = sys.argv[2]

with open(inFile, 'r') as f:
	header = f.readlines()[10:11]

newheader = re.sub('[^a-zA-Z0-9\n\.]', ' ', str(header)).replace(' s ', ' ').strip()[:-2].split()

newheader_list =[]

for i in newheader:
	if 'name' in i or 'tsuperkingdom.t.' in i or 'tphylum.t.' in i:
		newheader_list.append(i)

df = pd.DataFrame(columns=newheader_list)


with open(inFile, 'r') as f2:
	body = f2.readlines()[11:]


body_list = []
for i in body:
	line = i.replace('\t', ' ').split()
	
	for f, b in zip(newheader, line):
		if 'name' in f or 'tsuperkingdom.t.' in f or 'tphylum.t.' in f:
			body_list.append(b)
	df_length = len(df)
	df.loc[df_length] = body_list
	body_list =[]


contaminant_contig_list = []
for index, row in df.iterrows():
	if row['tsuperkingdom.t.6'] == 'Viruses' or row['tphylum.t.10'] == 'Streptophyta':
		contaminant_contig_list.append(row['name'])	

with open(outFile, 'a') as f:
	for i in contaminant_contig_list:
		f.write(i)
