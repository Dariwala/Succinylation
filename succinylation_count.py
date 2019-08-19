proteins = open("proteins_in_dataset.txt")

lines = proteins.readlines()

protein = []
for line in lines:
	protein.append(line[:-1])

succs = open("Saifur sir dataset//Succinylation.txt")

succs = succs.readlines()
prot_seq={}
succ_count = {}

for succ in succs[1:]:	
	succ_splitted = succ.split("\t")
	prot = succ_splitted[0]
	seq = succ_splitted[1]
	k_count = seq.count("K")
	prot_seq[prot] = k_count
	
	if prot not in succ_count:
		succ_count[prot] = 0
	succ_count[prot] += 1

succ_pos = 0
succ_neg = 0

for prot in succ_count:
	if prot in protein:
		succ_pos += succ_count[prot]
		succ_neg += prot_seq[prot] - succ_count[prot]

print(succ_pos,succ_neg)
	