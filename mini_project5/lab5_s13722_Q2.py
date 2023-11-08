"""
 Author : Ayesha Sanahari
 Date : 20/Jan/2021
 Perform Independent two-sample t-test using python
 Input: “HornedLizards.csv” data file
 Output: Results of Statistical test performed
"""

import scipy.stats as sci
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.graphics.gofplots import qqplot

dataFrame2 = pd.read_csv("./La5_datasets/HornedLizards.csv")

# Drop rows with missing values
dataFrame2.dropna(inplace=True)
#print(dataFrame2)

# Group By survive
df_by_survive = dataFrame2.groupby("Survive")
print("\n\n")
print(df_by_survive.describe().head())


# seperate into two variables by grouping variables
servived = dataFrame2.loc[dataFrame2['Survive'] == "survived"]
survived_data=servived['Squamosal horn length']
#print(survived_data)

dead = dataFrame2.loc[dataFrame2['Survive'] == "dead"]
dead_data=dead['Squamosal horn length']
#print(dead_data)

# QQ plot for servived
qqplot(survived_data, line='s')
plt.title('QQplot for survived lizards')
plt.savefig('Q2 QQplot for survived.jpg')

# QQ plot for dead
qqplot(dead_data, line='s')
plt.title('QQplot for dead lizards')
plt.savefig('Q2 QQplot for dead.jpg')


# Histogram of survived lizards data
fig,ax = plt.subplots(figsize = (12,5))
sns.histplot(survived_data, ax=ax, bins=30, kde=True)
plt.title("Histogram of survived lizards data")
plt.xlabel("Squamosal horn length")
plt.ylabel("No. of Survived lizards")
plt.savefig('Q2 histogram survived.jpg')

# Histogram of dead lizards data
fig,ax = plt.subplots(figsize = (12,5))
sns.histplot(dead_data, ax=ax, bins=30, kde=True)
plt.title("Histogram of dead lizards data")
plt.xlabel("Squamosal horn length")
plt.ylabel("No. of dead lizards")
plt.savefig('Q2 histogram dead.jpg')


# two histograms on one plot
plt.hist([survived_data, dead_data], label=['servived', 'dead'])
plt.title('Histogram of survived and dead lizards')
plt.legend(loc='upper left')
plt.xlabel("Squamosal horn length")
plt.ylabel("No. of lizards")
plt.savefig('Q2 two histo.jpg')

#Shapiro test for survived
stat, p = sci.shapiro(survived_data)
print('P value for survived:',p)

#Shapiro test for dead
stat, p = sci.shapiro(dead_data)
print('P value for dead:', p)

# draw boxplots for two samples
fig,ax =plt.subplots(1,2, figsize=(12,3))
ax[0].boxplot(survived_data)
ax[0].set_title("Boxplot of survived lizards horn length")
ax[1].boxplot(dead_data)
ax[1].set_title("Boxplot of dead lizards horn length")
plt.savefig("q2 boxplot.jpg")

# draw violin plots for two samples
fig,ax =plt.subplots(1,2, figsize=(12,3))
ax[0].violinplot(survived_data)
ax[0].set_title("Violin plot of survived lizards horn length")
ax[1].violinplot(dead_data)
ax[1].set_title("Violin plot of dead lizards horn length")
plt.savefig("q2 violin plot.jpg")

#independed sample t-test
stat,p = sci.kruskal(survived_data,dead_data)
print('P value for independent two sample t-test',p)
