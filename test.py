import subprocess

queue = str(subprocess.check_output('squeue', shell=True))
check = 0

for line in queue:
	if 'sykesj' in line:
		check += 1

print (check)