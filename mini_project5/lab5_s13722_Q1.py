"""
 Author : Ayesha Sanahari
 Date : 20/Jan/2021
 Perform one sample t test using python
 Input: Temperature.csv data file
 Output: Results of Statistical test performed
"""

import scipy.stats as sci
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import statsmodels.api as sm

dataFrame1 = pd.read_csv("./La5_datasets/Temperature.csv")
print(dataFrame1)
print(dataFrame1.describe())

# Draw Histogram
fig,ax = plt.subplots(figsize = (12,5))
sns.histplot(dataFrame1, ax=ax, bins=30, kde=True)
plt.title("Histogram for Body Temperature")
plt.xlabel("body Temperature (F)")
#plt.show()
plt.savefig('Q1 histogram.jpg')

# qqplot
plot2 = sm.qqplot(dataFrame1, line="s")
plt.title("Quantile-Quantile plot for body temperature")
plt.savefig('Q1 qqplot.jpg')
#plt.show()

#Shapiro test for normality test
stat,p =sci.shapiro(dataFrame1)
print("p value:", p)

#t test
checkvalue= 98.6
t,p = sci.ttest_1samp(dataFrame1,checkvalue)
print("t value:",t)
print("p value:",p)

