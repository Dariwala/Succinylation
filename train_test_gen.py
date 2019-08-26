import numpy as np
from sklearn.model_selection import train_test_split

dataset = open("dataset.txt")

dataset = dataset.readlines()

data = np.zeros((len(dataset),len(dataset[0].split("\t"))),dtype=float)

i = -1

for line in dataset:
	i += 1
	line = line.split("\t")
	j = -1
	for l in line:
		j += 1
		data[i][j] = float(l)

X = data[:,:-1]
y = data[:,-1]


X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=42)

output = open("train.txt","w")

for i in range(X_train):
	row = X_train[i]
	for col in row:
		output.write(str(col) + "\t")
	output.write(y_train[i]+"\n")

output = open("test.txt","w")

for i in range(X_test):
	row = X_test[i]
	for col in row:
		output.write(str(col) + "\t")
	output.write(y_test[i]+"\n")
	