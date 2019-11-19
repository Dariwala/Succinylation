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

def fetch_pssm_value(prot,indices):
	pssm_file = open("Saifur sir dataset//PSSM//"+prot+".fasta.pssm")
	pssm_file = pssm_file.readlines()
	
	pssm = np.zeros((31,20), dtype = float)
	
	for i in range(len(indices)):
		index = indices[i]
		line = pssm_file[index+3]
		line = line.split()
		for j in range(20):
			pssm[i][j] = float(line[22+j])/100 # 2 for first 20 values, 22 for last 22 values
			
		
	return pssm

def compute_pssm(pssm):
	B = np.zeros((20,20), dtype = float)
	
	for i in range(20):
		for j in range(20):
			for k in range(30):
				B[i][j] += pssm[k][i] * pssm[k+1][j]
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

output = open("final_output_only_pssm.txt",'w')

for prot in protein:
	fasta = open("Saifur sir dataset//FASTAs//"+prot+".fasta")
	lines = fasta.readlines()
	fasta = lines[1]
	
	for i in range(len(fasta)):
		if fasta[i] == 'K':
			indices = valid_indices(i,len(fasta))
			pssm_values = fetch_pssm_value(prot,indices)
			#spd3_values = fetch_spd3_value(prot,indices)
			
			B = compute_pssm(pssm_values)
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