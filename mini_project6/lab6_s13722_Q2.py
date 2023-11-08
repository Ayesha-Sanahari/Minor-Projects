"""
 Author : Ayesha Sanahari
 Date : 22/Jan/2021
 Use Python scikit-learn package to implement K-nearest neighbors(KNN) algorithm to learn from the popular iris data set.

"""
import numpy as np
from sklearn.datasets import load_iris
from sklearn import preprocessing
from sklearn.neighbors import KNeighborsClassifier

# Loading dataset into variables
iris = load_iris()
irisData = iris['data']
irisTargets = iris['target']
irisSpp = iris['target_names']

#scaling iris data
scalar = preprocessing.StandardScaler().fit(irisData)
scaledIris = scalar.transform(irisData)

#plant 1 data and scaling
plant1=np.array([[4.6, 3.0, 1.5, 0.2]])
scaledpalnt1= scalar.transform(plant1)

#both plants data and scaling
plants = np.array([[4.6, 3.0, 1.5, 0.2], [6.2, 3.0, 4.1, 1.2]])
scaledpalnts= scalar.transform(plants)

#selecting sepal collumn
sepal=scaledIris[:, :2]
testDataplant1 = scaledpalnt1[:, :2]
testDataPlants = scaledpalnts[:, :2]


#for sepal data
#KNeighborsClassifier for 2 neighbors
classifier = KNeighborsClassifier(n_neighbors=2).fit(sepal, irisTargets)

#print stats
for test in testDataplant1:
    print("Predicted species:")
    for i in classifier.predict(testDataplant1):
        print(irisSpp[i], end='\n')
    print("Nearest Neightbours:")
    for n in classifier.kneighbors([test], 2)[1][0]:
        print("index:", n, end='\t')
        print("coordinates:", sepal[n], end='\t')
        tmp =  classifier.predict([sepal[n]])[0]
        print('label:', tmp, end='\t')
        print('species:', irisSpp[tmp])

##KNeighborsClassifier for 5 neighbors
print("for 5 neighbours")
classifier = KNeighborsClassifier(n_neighbors=5).fit(sepal, irisTargets)
for test in testDataPlants:
    print("Predicted species:")
    for i in classifier.predict(testDataPlants):
        print(irisSpp[i], end='\n')
    print("Nearest Neightbours:")
    for n in classifier.kneighbors([test], 5)[1][0]:
        print("index:", n, end='\t')
        print("coordinates:", sepal[n], end='\t')
        tmp =  classifier.predict([sepal[n]])[0]
        print('label:', tmp, end='\t')
        print('species:', irisSpp[tmp])

#predictions
print("Predictions:")
for pred_species, pred_probability in zip(classifier.predict(testDataPlants),
                                          classifier.predict_proba(testDataPlants)):
    print(pred_species, end='\t')
    print(irisSpp[pred_species], end='\t')
    print(f'(probability: {pred_probability})')


#for petal data
print("\n\n for petal data")
#select petal data column
testDataplant1 = scaledpalnt1[:, 2:]
testDataPlants = scaledpalnts[:, 2:]
petal = scaledIris[:, 2:]
print("with 2 n_neighbors ")

#KNeighborsClassifier for 2 neighbors
classifier = KNeighborsClassifier(n_neighbors=2).fit(petal, irisTargets)

#print result
for test in testDataplant1:
    print("Predicted species:")
    for i in classifier.predict(testDataplant1):
        print(irisSpp[i], end='\n')
    print("Nearest Neightbours:")
    for n in classifier.kneighbors([test], 2)[1][0]:
        print("index:", n, end='\t')
        print("coordinates:", petal[n], end='\t')
        tmp =  classifier.predict([petal[n]])[0]
        print('label:', tmp, end='\t')
        print('species:', irisSpp[tmp])

print("for 5 neighbours")
#KNeighborsClassifier for 2 neighbors and both plants
classifier = KNeighborsClassifier(n_neighbors=5).fit(petal, irisTargets)

#print result
for test in testDataPlants:
    print("Predicted species:")
    for i in classifier.predict(testDataPlants):
        print(irisSpp[i], end='\n')
    print("Nearest Neightbours:")
    for n in classifier.kneighbors([test], 5)[1][0]:
        print("index:", n, end='\t')
        print("coordinates:", petal[n], end='\t')
        tmp =  classifier.predict([petal[n]])[0]
        print('label:', tmp, end='\t')
        print('species:', irisSpp[tmp])

#predictions
print("Predictions:")
testDataplant1 = scaledpalnt1[:, 2:]
for pred_species, pred_probability in zip(classifier.predict(testDataPlants),
                                          classifier.predict_proba(testDataPlants)):
    print(pred_species, end='\t')
    print(irisSpp[pred_species], end='\t')
    print(f'(probability: {pred_probability})')

#using both petal and sepal data
print("\n\n for sepal and petal data")
testDataPlants = scaledpalnts[:, :]
data = scaledIris[:, :]
testDataplant1 = scaledpalnt1[:, :]
print("with 2 n_neighbors ")

#KNeighborsClassifier for all data
classifier = KNeighborsClassifier(n_neighbors=2).fit(data, irisTargets)

#print result
for test in testDataplant1:
    print("Predicted species:")
    for i in classifier.predict(testDataplant1):
        print(irisSpp[i], end='\n')
    print("Nearest Neightbours:")
    for n in classifier.kneighbors([test], 2)[1][0]:
        print("index:", n, end='\t')
        print("coordinates:", data[n], end='\t')
        tmp =  classifier.predict([data[n]])[0]
        print('label:', tmp, end='\t')
        print('species:', irisSpp[tmp])

#KNeighborsClassifier for all data with 5 n_neighbors
print("for 5 neighbours")
classifier = KNeighborsClassifier(n_neighbors=5).fit(data, irisTargets)
for test in testDataPlants:
    print("Predicted species:")
    for i in classifier.predict(testDataPlants):
        print(irisSpp[i], end='\n')
    print("Nearest Neightbours:")
    for n in classifier.kneighbors([test], 5)[1][0]:
        print("index:", n, end='\t')
        print("coordinates:", data[n], end='\t')
        tmp =  classifier.predict([data[n]])[0]
        print('label:', tmp, end='\t')
        print('species:', irisSpp[tmp])


print("Predictions:")
for pred_species, pred_probability in zip(classifier.predict(testDataPlants),
                                          classifier.predict_proba(testDataPlants)):
    print(pred_species, end='\t')
    print(irisSpp[pred_species], end='\t')
    print(f'(probability: {pred_probability})')