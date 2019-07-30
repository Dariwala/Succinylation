import os

protein_file = open("Cd hit 40% result\\1564510281.fas.1")
lines = protein_file.readlines()
proteins = open("proteins_in_dataset.txt",'w')

for i in range(len(lines)):
	if i%4 == 0:
		print(lines[i][1:-1])
		proteins.write(lines[i][1:-1]+"\n")