data = open("dataset.txt",'r')

proteins = data.readlines()
unique_protein = {}
cplm = open("plmd.txt",'w')

for protein in proteins:
    protein = protein.split("\t")
    if protein[1] not in unique_protein:
        unique_protein[protein[1]] = 1
        fasta_file = open("Fastas//"+protein[1]+".txt",'w')
        fasta_file.write(">"+protein[1]+"\n"+protein[4]+"\n")
        cplm.write(">"+protein[1]+"\n"+protein[4]+"\n")
        fasta_file.close()
cplm.close()
