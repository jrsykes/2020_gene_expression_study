 #

import requests
from bs4 import BeautifulSoup

page = requests.get('https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1799852')
soup = BeautifulSoup(page.content, 'html.parser')
raw_dat = soup.find(id='sra-viewer-app')

dat = list((raw_dat.find_all(class_='Read')))


for item in dat:
	if 'σ' in str(item):
		print (item)

#import re
#intro = "<>I'm Tom."
#intro = re.sub(r'.*I', 'I', intro)
#print(intro)