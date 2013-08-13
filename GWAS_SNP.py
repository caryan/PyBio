"""
Module for some GWAS attempts.
"""

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact

import matplotlib.pyplot as plt
import subprocess
from cStringIO import StringIO
from glob import glob



#Load all the SNP's for each file
fileNames = ['ARC4353.txt', 'ARC4359.txt']
fileNames = glob('ARC*.txt')
drugName = 'imipenem'
resistCutoff = 4
resistDF = pd.read_csv('resistances.txt', sep='\t', index_col='Strain')
strainNames = set([])
dfs = {}
for fileName in fileNames:
	strainName = fileName.split('.')[0]
	strainNames.add(strainName)
	#Use sed to clean-up the insertions/deletions and load into a dataframe
	dfs[strainName] = pd.read_csv(
			StringIO(subprocess.check_output(["sed", "-e", r"s/\^[0-9]*//g", "-e", r"s/\.\.[0-9]*//g",  fileName])), 
			sep='\t', usecols=['Region', 'Type', 'Reference', 'Allele'])
	#Inject the SNP key, strain name and resistance
	dfs[strainName]['strainName'] = strainName
	dfs[strainName]['resistance'] = resistDF[drugName][strainName] > resistCutoff

#Concatenate all the data frames into one big one for grouping and  pivoting
bigDF = pd.concat(dfs.values(), ignore_index=True)

#Group by SNP's
grouped = bigDF.groupby(['Region', 'Type', 'Reference', 'Allele'])

#For each SNP create the 2x2 contingency matrix and estimate the p-value
def calc_pValue(group):
	#Number with SNP and resistance
	withSNP = len(group)
	withSNP_R = np.count_nonzero(group['resistance'].values)
	withSNP_S = withSNP - withSNP_R
	withoutSNP = strainNames - set(group['strainName'].values)
	withoutSNP_R = 0; withoutSNP_S = 0;
	for strainName in withoutSNP:
		if resistDF[drugName][strainName] > resistCutoff:
			withoutSNP_R += 1
		else:
			withoutSNP_S += 1
	oddsRatio, pValue = fisher_exact([[withSNP_R, withSNP_S], [withoutSNP_R, withoutSNP_S]]) 
	return pValue

#Calculate the p Value for each SNP
SNPNames = []
SNPPos = []
pValues = []
for k, v in grouped:
	#Maybe filter for those that are at least present in two strains?
	if len(v) > 1:
		SNPNames.append('-'.join([str(x) for x in k]))
		SNPPos.append(k[0])
		pValues.append(calc_pValue(v))



#Make the Manhattan plot
def onpick(event):
	print(SNPNames[event.ind[0]])

figH = plt.figure()
axesH = figH.add_subplot(111)
axesH.scatter(SNPPos, -np.log10(pValues), picker=True)

#Plot the 0.05 significance line
#TODO: apply multiple testing critera
axesH.axhline(-np.log10(0.05), linestyle='--', color='k')

figH.canvas.mpl_connect('pick_event', onpick)

plt.show()

#TODO: plot the quantile-quantile plot to look for population stratification




	

