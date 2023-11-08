"""
 Author : Ayesha Sanahari
 Date : 22/Jan/2021
 use Python scikit-learn package to implement K-means algorithms to learn from the popular iris data set.
 Input: iris.csv data file
 Output: scatter plot of iris data with predictions
"""
import numpy as np
from sklearn.cluster import KMeans
import pandas as pd
import matplotlib.pyplot as plt

#reading data from pandas
iris=pd.read_csv('iris.csv')
#define label column
df2=iris['v_short']
#print(df2)

# Add sepal data columns to Dataframe
df=pd.DataFrame(iris,columns=['sepal.length','sepal.width'])
df3=pd.DataFrame(iris,columns=['v_short'])
#print(df)

kmeans= KMeans(n_clusters=3).fit(df)
# getting centroids
centroids=kmeans.cluster_centers_
print("kmeans cluster_centers:", centroids)

# specise unknown data sample
plant=np.array([[4.6, 3.0, 1.5, 0.2],[6.2, 3.0, 4.1, 1.2]])
X = plant[:, [0, 1]]

# Predicting the lables of new inserted data
predicted=kmeans.predict(X)
print("predicted lables for sepals data:", predicted)
labels=kmeans.labels_
#print("labels:", labels)

#scatter plot of sepal length ,sepal width ,unknown saple data with cetroids
plt.scatter(df['sepal.length'],df['sepal.width'],c=kmeans.labels_.astype(float),s=50,alpha=0.5)
plt.scatter(X[:,0],X[:,1],c='green',s=50)
plt.scatter(centroids[:,0],centroids[:,1],c='red',s=50)
plt.title("Scatter plot for Sepal measurement")
plt.xlabel('Sepal length')
plt.ylabel('Sepal width')


#annotation
for df2, x, y,z in zip(df2, df.iloc[:, 0], df.iloc[:, 1],df3.iloc[:,0]):
    if z=='s':
        plt.annotate(
            df2,
            xy=(x, y), xytext=(-50, 20),
            textcoords='offset points', ha='right', va='bottom',
            bbox=dict(boxstyle='circle,pad=0.3',fc="white", alpha=0.2),
            arrowprops=dict(arrowstyle = '->', connectionstyle="arc,rad=40",alpha=0.2))
    if z == 've':
        plt.annotate(
            df2,
            xy=(x, y), xytext=(50, -50),
            textcoords='offset points', ha='right', va='bottom',
            bbox=dict(boxstyle='circle,pad=0.3', fc="pink", alpha=0.2),
            arrowprops=dict(arrowstyle='->', connectionstyle="arc,rad=40", alpha=0.2))
    if z == 'vi':
        plt.annotate(
            df2,
            xy=(x, y), xytext=(50, 50),
            textcoords='offset points', ha='right', va='bottom',
            bbox=dict(boxstyle='circle,pad=0.3', fc="white", alpha=0.2),
            arrowprops=dict(arrowstyle='->', connectionstyle="angle,angleA=-90,angleB=180,rad=0", alpha=0.2))

#label annotation
for label, x, y in zip(labels, df.iloc[:, 0], df.iloc[:, 1]):
   plt.annotate(
       label,
       xy=(x, y), xytext=(1, 1.5),
       textcoords='offset points')
plt.show()

print("\n\n for the petals data")
#define label collumns
df2=iris['v_short']

# Add petal data columns to Dataframe
df4=pd.DataFrame(iris,columns=['petal.length','petal.width'])
#print(df4)

kmeans= KMeans(n_clusters=3).fit(df4)
# getting centroids
centroids=kmeans.cluster_centers_
print("kmeans cluster_centers:", centroids)

# specise unknown data sample
plant=np.array([[4.6, 3.0, 1.5, 0.2],[6.2, 3.0, 4.1, 1.2]])
X = plant[:, [2, 3]]

# Predicting the lables of new inserted data
predicted=kmeans.predict(X)
print("predicted lables for petals data:", predicted)
labels=kmeans.labels_
#print("labels:", labels)

#scatter plot of sepal length ,sepal width ,unknown saple data with cetroids
plt.scatter(df4['petal.length'],df4['petal.width'],c=kmeans.labels_.astype(float),s=50,alpha=0.5)
plt.scatter(X[:,0],X[:,1],c='green',s=50)
plt.scatter(centroids[:,0],centroids[:,1],c='red',s=50)
plt.title("Scatter plot for Petal measurement")
plt.xlabel('Petal length')
plt.ylabel('Petal width')


#annotation
for df2, x, y,z in zip(df2, df4.iloc[:, 0], df4.iloc[:, 1],df3.iloc[:,0]):
    if z=='s':
        plt.annotate(
            df2,
            xy=(x, y), xytext=(-50, 20),
            textcoords='offset points', ha='right', va='bottom',
            bbox=dict(boxstyle='circle,pad=0.3',fc="white", alpha=0.2),
            arrowprops=dict(arrowstyle = '->', connectionstyle="arc,rad=40",alpha=0.2))
    if z == 've':
        plt.annotate(
            df2,
            xy=(x, y), xytext=(50, -50),
            textcoords='offset points', ha='right', va='bottom',
            bbox=dict(boxstyle='circle,pad=0.3', fc="pink", alpha=0.2),
            arrowprops=dict(arrowstyle='->', connectionstyle="arc,rad=40", alpha=0.2))
    if z == 'vi':
        plt.annotate(
            df2,
            xy=(x, y), xytext=(50, 50),
            textcoords='offset points', ha='right', va='bottom',
            bbox=dict(boxstyle='circle,pad=0.3', fc="white", alpha=0.2),
            arrowprops=dict(arrowstyle='->', connectionstyle="angle,angleA=-90,angleB=180,rad=0", alpha=0.2))

#label annotation
for label, x, y in zip(labels, df4.iloc[:, 0], df4.iloc[:, 1]):
   plt.annotate(
       label,
       xy=(x, y), xytext=(1, 1.5),
       textcoords='offset points')
plt.show()

