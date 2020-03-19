import subprocess
import re
import pandas as pd
from io import StringIO


check = str(subprocess.check_output('squeue --user=sykesj', shell=True), 'utf-8')


data = StringIO(check)

df = pd.read_csv(data)

for row in df.iterrows():
	print (row)
