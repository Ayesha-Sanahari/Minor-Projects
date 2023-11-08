"""
 Author : Ayesha Sanahari
 Date : 22/Jan/2021
 Perform Chi-square contingency test using python
"""
import scipy.stats as sci
import matplotlib.pyplot as plt
import pandas as pd
import statsmodels.graphics.mosaicplot as mosaicplot

# Creating a Dataframe from the data and printing it
data = {'Uninfected': [1, 49],'Lightly infected': [10, 35],'Highly infected':[37,9]}
dataFrame4 = pd.DataFrame(data, index=['eaten by birds','not eaten by birds'])
print(dataFrame4)

# plot a mosaicplot for the dataframe
mosaicplot.mosaic(dataFrame4.stack(),title='mosaic plot')
plt.show()

# do a chi-square contingency test
chi_value,p,dof,expected = sci.chi2_contingency(dataFrame4)
print('Chi values :',chi_value)
print('P value :',p)
print('degree of freedom :',dof)

expected = pd.DataFrame(expected)
print(expected)