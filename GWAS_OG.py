"""
Module for some GWAS attempts.
"""

import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency

import matplotlib.pyplot as plt

#For each SNP create the 2x2 contingency matrix and estimate the p-value
def calc_pValue(group):
	#Number with NSAC and resistant
	counts = group.value_counts()

	#Do the susceptable first
	withoutNSAC_S = counts[0][0] if 0 in counts[0] else 0
	withNSAC_S = counts[0][1] if 1 in counts[0] else 0
	withBonus_S = counts[0][2] if 2 in counts[0] else 0

	withoutNSAC_R = counts[1][0] if 0 in counts[1] else 0
	withNSAC_R = counts[1][1] if 1 in counts[1] else 0
	withBonus_R = counts[1][2] if 2 in counts[1] else 0

	#If we have no 2 counts then ch2_contingency bails
	if withBonus_R or withBonus_S:
		chi2, pValue, dof, expected = chi2_contingency([[withoutNSAC_S, withNSAC_S, withBonus_S], [withoutNSAC_R, withNSAC_R, withBonus_R]])
	else:
		chi2, pValue, dof, expected = chi2_contingency([[withoutNSAC_S, withNSAC_S], [withoutNSAC_R, withNSAC_R]])

	return pValue

#Load the big matrix from file
#Transpose to make it easy group by resistance
df = pd.read_csv('/home/cryan//Veronica/Regression/allporinslevodecision.txt', sep='\t', index_col=0).transpose()

#Group by resistance and sum the number that are different from PAO1 
grouped = df.groupby('resistant')

#The resistant will just be a pain iterating below
df.pop('resistant')

#Calculate the p Value for each location
names = []
pValues = []
for nsac in df.columns:
	names.append(nsac)
	pValues.append(calc_pValue(grouped[nsac]))


#Make the Manhattan plot
def onpick(event):
	print(names[event.ind[0]])

plt.figure()
plt.imshow(df.values, aspect='auto')

figH = plt.figure()
axesH = figH.add_subplot(111)
axesH.scatter(np.arange(len(names)), -np.log10(pValues), picker=True)

#Plot the 0.05 significance line
#TODO: apply multiple testing critera
axesH.axhline(-np.log10(0.05), linestyle='--', color='k')

figH.canvas.mpl_connect('pick_event', onpick)

plt.show()

#TODO: plot the quantile-quantile plot to look for population stratification




	

