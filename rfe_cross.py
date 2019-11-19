from sklearn.svm import SVC
from sklearn.feature_selection import RFE
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold
import sys

no_of_folds = -1

dataset = open("dataset_pssm_spd3.txt","r")

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
for pop in [5,80,100,150,200,250,300,350]:
	for no_of_folds in [6,8,10]:
		for r in [3,19,11,1001]:
			accuracy = sensitivity = specificity = mcc = 0.0
			kf = KFold(n_splits=no_of_folds,shuffle=True,random_state = r)
			for train,test in kf.split(data):
				X_train, X_test, y_train, y_test = X[train], X[test], y[train], y[test]

				clf = SVC(gamma='auto',kernel='linear')
				selector = RFE(clf, pop, step=1)
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

				accuracy += (curr_pos + curr_neg)/(curr_pos + curr_neg + inc_neg + inc_pos)
				sensitivity += (curr_pos)/(curr_pos + inc_neg)
				specificity += (curr_neg)/(curr_neg + inc_pos)
				mcc += (curr_pos * curr_neg - inc_neg*inc_pos)/(np.sqrt((curr_pos + inc_pos)*(curr_pos + inc_neg)*(curr_neg+inc_pos)*(curr_neg+inc_neg)))
			print(pop,r,no_of_folds)
			print(accuracy/no_of_folds,specificity/no_of_folds,sensitivity/no_of_folds,mcc/no_of_folds)