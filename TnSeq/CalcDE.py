# Copyright 2013 Colm Ryan colm@colmryan.org
# License GPL v3 (http://www.gnu.org/licenses/gpl.txt) 

'''
Script to calculate real dVals
'''
from __future__ import division

import numpy as np
import pandas as pd
import rpy2.robjects as rBridge

#File names

#Control
controlName = 'MH'

#Drugs
drugNames = ['Vancomycin']

#First let's use the the all insertion site sheet to figure out the expected read count density
allGenesControlWB_A = pd.ExcelFile(controlName + '_Aggregate_A.xlsx')
allGenesControlWB_B = pd.ExcelFile(controlName + '_Aggregate_B.xlsx')

geneDataControl_A = allGenesControlWB_A.parse(allGenesControlWB_A.sheet_names[0], index_col=1, na_values=['']).drop('intergenic')
geneDataControl_B = allGenesControlWB_B.parse(allGenesControlWB_B.sheet_names[0], index_col=1, na_values=['']).drop('intergenic')


#Let's first look at read density (reads/length) to determine overal trend 

#Remove the top/bottom 10% spikes
qcuts = pd.qcut(geneDataControl_A['ReadCount'], [0, 0.1, 0.9, 1])
goodCts = geneDataControl_A[qcuts==qcuts.levels[1]]

goodCts['readDensity'] = goodCts['ReadCount']/goodCts['Length']

# #Push the data into R for the loess fit
rBridge.r.assign('positions', rBridge.IntVector(goodCts['Start'].values))
rBridge.r.assign('counts', rBridge.FloatVector(goodCts['readDensity'].values))
rBridge.r('loessFit <- loess(counts~positions, data.frame(positions=positions, counts=counts), control = loess.control(statistics = c("approximate"),trace.hat = c("approximate"), iterations=3))')
rBridge.r('smoothed <- predict(loessFit, positions)')

goodCts['smoothed'] = np.array(rBridge.r('smoothed'), dtype=np.float64)

# # # #Now bin the positions into 1kb pair regions weighted by the read counts to estimate the local read density
# xPts = np.arange(0,goodCts['Position'].max(),1000)
# densityEstimate = np.histogram(goodCts['Position'], bins=xPts, weights=goodCts['TotalCount'])[0]
# densityLoess = np.histogram(goodCts['Position'], bins=xPts, weights=goodCts['smoothed'])[0]


# #Fit this to a quadratic (should) probably update to some non-parametric approach like LOESS)
# fitDensity = np.poly1d(np.polyfit(xPts[:-1], densityEstimate,2))

# #Now we can load the gene data and calculate expected reads for each gene
# aggregateSheet = open_workbook(aggregateFile, on_demand=True).sheet_by_index(0)

# geneNameCol = 2
# lengthCol = 11
# startCol = 8
# ctCol = 7

# df = pd.DataFrame({	'start':[int(tmpValue) if tmpValue is not '' else 0 for tmpValue in aggregateSheet.col_values(startCol,1)],
# 					'length':[int(tmpValue) if tmpValue is not '' else 0 for tmpValue in aggregateSheet.col_values(lengthCol,1)],
# 					'readCounts':[int(tmpValue) if tmpValue is not '' else 0 for tmpValue in aggregateSheet.col_values(ctCol, 1)]},
# 					index=aggregateSheet.col_values(geneNameCol, 1) )


# geneMidPts = df['start']+df['length']/2

# df['expectedCts'] = df['length']/1000.0*fitDensity(geneMidPts)

# df['dVals'] = (df['readCounts']/df['expectedCts'])

# foldFactor = np.log2(df['dVals'])