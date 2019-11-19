import os
import numpy as np

def valid_indices(index,length):
	indices = []
	if index < 15:
		for i in range(1,16):
			indices.append(index+i)
		indices.append(index)
		for i in range(1,16):
			indices.append(index+i)
	elif index >= length - 15:
		for i in range(15,0,-1):
			indices.append(index-i)
		indices.append(index)
		for i in range(15,0,-1):
			indices.append(index-i)
	else:
		for i in range(15,0,-1):
			indices.append(index-i)
		indices.append(index)
		for i in range(1,16):
			indices.append(index+i)
	return indices

def count(str,substr,indices):
	co = 0
	for i in range(0,len(indices)-len(substr)):
		is_equal = True
		for j in range(i,i+len(substr)):
			if substr[j-i] != str[indices[j]]:
				is_equal = False
		if is_equal:
			co += 1
	return co

def compute_pssm(prot,indices):
	#print(len(prot),len(indices))
	acids = "ARNDCQEGHILKMFPSTWYV"
	B = np.zeros((8,1), dtype = float)
	
	B[0][0] = count(prot,"CFILMVW",indices)/25
	B[1][0] = count(prot,"NQSTY",indices)/27
	B[2][0] = count(prot,"DEKR",indices)/28
	B[3][0] = count(prot,"AG",indices)/30
	B[4][0] = count(prot,"HP",indices)/30
	B[5][0] = count(prot,"DE",indices)/30
	B[6][0] = count(prot,"IL",indices)/30
	B[7][0] = count(prot,"NQ",indices)/30
	'''B[8][0] = count(prot,"NQSTY",indices)/27
	B[9][0] = count(prot,"NQSTY",indices)/27
	B[2][0] = count(prot,"NQSTY",indices)/27
	B[2][0] = count(prot,"NQSTY",indices)/27
	B[2][0] = count(prot,"NQSTY",indices)/27
	B[2][0] = count(prot,"NQSTY",indices)/27
	B[2][0] = count(prot,"NQSTY",indices)/27
	B[2][0] = count(prot,"NQSTY",indices)/27'''
	return B


protein_file = open("proteins_in_dataset.txt")

protein_file = protein_file.readlines()

protein = []

for line in protein_file:
	protein.append(line[:-1])
	
succ_file = open("Saifur sir dataset//Succinylation.txt")
succ_file = succ_file.readlines()

succ_indices = {}

for line in succ_file[1:]:
	l = line.split("\t")
	prot = l[0]
	if prot in protein:
		if prot not in succ_indices:
			succ_indices[prot] = []
		succ_indices[prot].append(l[2])

output = open("final_output_special_bigram.txt",'w')

for prot in protein:
	print(prot)
	fasta = open("Saifur sir dataset//FASTAs//"+prot+".fasta")
	lines = fasta.readlines()
	fasta = lines[1]
	
	for i in range(len(fasta)):
		if fasta[i] == 'K':
			indices = valid_indices(i,len(fasta))
			#pssm_values = fetch_pssm_value(prot,indices)
			#spd3_values = fetch_spd3_value(prot,indices)
			
			B = compute_pssm(fasta,indices)
			#Bp = compute_spd3(spd3_values)
			
			B = B.flatten()
			#Bp = Bp.flatten()
			#total = np.concatenate([B,Bp])
			
			label = 0
			
			if str(i+1) in succ_indices[prot]:
				label = 1
			
			for j in range(len(B)):
				output.write(str(B[j]) + "\t")
			output.write(str(label)+"\n")
			
output.close()