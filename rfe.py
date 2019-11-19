from sklearn.svm import SVC
from sklearn.feature_selection import RFE
import numpy as np
from sklearn.model_selection import train_test_split

dataset = open("dataset_pssm_spd3_01.txt","r")

dataset = dataset.readlines()

data = np.zeros((len(dataset),len(dataset[0].split("\t")[:-1])))

i = -1
for line in dataset:
	i += 1
	l = line.split("\t")
	j = -1
	for val in l[:-1]:
		j += 1
		data[i][j] = float(val)

X, y = data[:,:-1],data[:,-1]
#y = np.array([y])
#y = np.reshape(y,(y.shape[1],y.shape[0]))

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33,stratify=y, random_state=42)

clf = SVC(gamma='auto',kernel='linear')
selector = RFE(clf, 100, step=1)
selector = selector.fit(X_train, y_train)

y_pred = selector.estimator_.predict(X_test.compress(selector.support_,axis=1))

curr_pos = curr_neg = inc_pos = inc_neg = 0

for i in range(len(y_test)):
	if y_test[i] == 1:
		if y_pred[i] == 1:
			curr_pos += 1
		else:
			inc_neg += 1
	else:
		if y_pred[i] == 1:
			inc_pos += 1
		else:
			curr_neg += 1

accuracy = (curr_pos + curr_neg)/(curr_pos + curr_neg + inc_neg + inc_pos)
sensitivity = (curr_pos)/(curr_pos + inc_neg)
specificity = (curr_neg)/(curr_neg + inc_pos)
mcc = (curr_pos * curr_neg - inc_neg*inc_pos)/(np.sqrt((curr_pos + inc_pos)*(curr_pos + inc_neg)*(curr_neg+inc_pos)*(curr_neg+inc_neg)))
print(accuracy,specificity,sensitivity,mcc)