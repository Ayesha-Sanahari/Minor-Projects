"""
 Author : Ayesha Sanahari
 Date : 23/Jan/2021
 Perform Independent two-sample t-test using python
 Input: “BlackbirdTestosterone.csv” data file
 Output: Results of Statistical test performed
"""

import scipy.stats as sci
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from statsmodels.graphics.gofplots import qqplot

#read the datafile
dataFrame3 = pd.read_csv('./La5_datasets/BlackbirdTestosterone.csv')
print(dataFrame3)

# seperate columns
log_before=dataFrame3['log before']
log_after=dataFrame3['log after']
log_diff=dataFrame3['dif in logs']

print('log before:', log_before.describe())
print('log after:', log_after.describe())
print('log diff:', log_diff.describe())

# draw qq plot and histogram for log difference
qqplot(log_diff, line='s')
plt.title('QQplot for log difference')
plt.show()

plt.hist(log_diff)
plt.title("Histogram for log difference")
plt.xlabel('log difference value')
plt.ylabel('Count')
plt.show()

#Shapiro test for Normality
stat, p = sci.shapiro(log_diff)
print('P value of log difference:',p)

stat, p = sci.shapiro(log_before)
print('P value of log before:',p)
stat, p = sci.shapiro(log_after)
print('P value of log after:',p)

# Draw boxplots
fig,ax =plt.subplots(1,2, figsize=(12,3), sharey=True)
ax[0].boxplot(log_before)
ax[0].set_title("Boxplot of before high level of testosterone")
ax[1].boxplot(log_after)
ax[1].set_title("Boxplot of after high level of testosterone")
plt.show()

# Draw violin plots
fig,ax =plt.subplots(1,2, figsize=(12,3), sharey=True)
ax[0].violinplot(log_before)
ax[0].set_title("Violin plot of before high level of testosterone")
ax[1].violinplot(log_after)
ax[1].set_title("Violin plot of after high level of testosterone")
plt.show()

# Paired sample t-test
stat,p = sci.ttest_rel(log_before,log_after)
print('P value for paired sample t-test',p)