import os

files = os.listdir("Saifur sir dataset\\FASTAs\\")
fastas = open("fastas.txt",'w')

for file in files:
	fasta = open("Saifur sir dataset\\FASTAs\\" + file)
	lines = fasta.readlines()
	fastas.write(lines[0])
	fastas.write(lines[1] + "\n")
	#break
	