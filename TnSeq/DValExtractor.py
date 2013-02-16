from __future__ import division
 
import numpy as np
import pandas as pd

fileDirectory = '/home/cryan/Desktop/Colm NCTC8325/NEWAGG/'

drugList = ['MH', 'Cip', 'Ox']
DEFile = '/home/cryan/Desktop/Colm NCTC8325/DE.xlsx'
geneInfo = {}
essentials = {}
conditionals = {}

#Make a general DataFrame with all the data
geneInfo = pd.read_table(fileDirectory+'MHA.rdp', index_col='Locus').fillna(0).drop('intergenic')
geneInfo = geneInfo[['Length', 'ProteinID', 'Product']]

for drugName in drugList:
	#Load the worksheets
	tmpWorkSheetA = pd.read_table(fileDirectory+drugName+'A.rdp', index_col='Locus').fillna(0).drop('intergenic')
	tmpWorkSheetB = pd.read_table(fileDirectory+drugName+'B.rdp', index_col='Locus').fillna(0).drop('intergenic')

	#Calculate the average DVal
	geneInfo[drugName+'_DVal'] = 0.5*(tmpWorkSheetA['DvalGenome']+tmpWorkSheetB['DvalGenome'])

	#Calculate the average TN density from repeats
	geneInfo[drugName+'_TNDensity'] = 0.5*(tmpWorkSheetA['Sites']/tmpWorkSheetA['Length'] + tmpWorkSheetB['Sites']/tmpWorkSheetB['Length'])	

	# #Now filter into ECL
	geneInfo[drugName+'_ECL'] = 'L'
	geneInfo[drugName+'_ECL'][np.logical_or(geneInfo[drugName+'_DVal']<0.01, geneInfo[drugName+'_TNDensity']<1/430)] = 'E'
	geneInfo[drugName+'_ECL'][np.logical_and(geneInfo[drugName+'_DVal']>0.01, geneInfo[drugName+'_DVal']<0.1)] = 'C'

	#Add the sensitized and supergrower flags	
	if drugName != 'MH':
		geneInfo[drugName+'_Sensitized'] = np.logical_and(geneInfo['MH_ECL'] == 'L', np.logical_or(geneInfo[drugName+'_ECL'] == 'E', geneInfo[drugName+'_ECL'] == 'C'))
		geneInfo[drugName+'_SuperGrower'] = np.logical_and(geneInfo[drugName+'_ECL'] == 'L', np.logical_or(geneInfo['MH_ECL'] == 'E', geneInfo['MH_ECL'] == 'C'))


excelWriter = pd.ExcelWriter('DVals.xlsx')
excelReader = pd.ExcelFile(DEFile)

for drugName in drugList:
	if drugName != 'MH':
		#Pull out the pvalues from the DE file
		tmpDE = excelReader.parse(drugName, index_col='locus')
		tmpDF = geneInfo[['Length', 'ProteinID', 'Product', 'MH_DVal', 'MH_TNDensity', 'MH_ECL', drugName+'_DVal', drugName+'_TNDensity', drugName+'_ECL']]
		tmpDF[drugName+'_logFC'] = tmpDE['logFC']
		tmpDF[drugName+'_pValue_adj'] = tmpDE['PValue_adj']
		tmpDF.to_excel(excelWriter, sheet_name='MH vs '+drugName)
		tmpDF[geneInfo[drugName+'_Sensitized']].drop(['MH_ECL', drugName+'_ECL'], axis=1).to_excel(excelWriter, sheet_name=drugName+' Sensitized')
		tmpDF[geneInfo[drugName+'_SuperGrower']].drop(['MH_ECL', drugName+'_ECL'], axis=1).to_excel(excelWriter, sheet_name=drugName+' SuperGrower')
		
excelWriter.save()
