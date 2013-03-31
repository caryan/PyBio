# Copyright 2013 Colm Ryan colm@colmryan.org
# License GPL v3 (http://www.gnu.org/licenses/gpl.txt) 

from __future__ import division
 
#Let's write to SVG style graphics
#import matplotlib
#matplotlib.use('svg')

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.collections as mcollections

import numpy as np

from xlrd import open_workbook

#Worksheets with data
hopCountsFile_SampA = '/home/cryan/Desktop/Colm NCTC8325/Vancomycin Ring A Galaxy1174-(Hopcount_on_data_99_and_data_1172__Hop_Table).tabular.xlsx'
hopCountsFile_SampB = '/home/cryan/Desktop/Colm NCTC8325/vancomycin ring B Galaxy1176-(Hopcount_on_data_99_and_data_1173__Hop_Table).tabular.xlsx'
hopCountsAggr = '/home/cryan/Desktop/Colm NCTC8325/Vancomycin ring AGalaxy1178-(Aggregate_Hop_Table_on_data_99_and_data_1174__Result_Table).tabular.xlsx' 

#Region of interest
ROI = [0, 200000]

#Minimum number of counts to plot
minPlotCts = 10


sheetA = open_workbook(hopCountsFile_SampA).sheet_by_index(0);
sheetB = open_workbook(hopCountsFile_SampB).sheet_by_index(0);
sheetAggr = open_workbook(hopCountsAggr).sheet_by_index(0);

#Find the full set of insertion sites
colHeaders = [tmpCell.value for tmpCell in sheetA.row(0)]
posCol = colHeaders.index(u'Position')
plusCtsCol = colHeaders.index(u'PlusCount')
minusCtsCol = colHeaders.index(u'MinusCount')

insertionSitesA = [int(tmpValue) for tmpValue in sheetA.col_values(posCol,1)]
plusCtsA = {tmpSite:int(tmpValue) for tmpSite,tmpValue in zip(insertionSitesA,sheetA.col_values(plusCtsCol,1))}
minusCtsA = {tmpSite:int(tmpValue) for tmpSite,tmpValue in zip(insertionSitesA,sheetA.col_values(minusCtsCol,1))}

insertionSitesB = [int(tmpValue) for tmpValue in sheetB.col_values(posCol,1)]
plusCtsB = {tmpSite:int(tmpValue) for tmpSite,tmpValue in zip(insertionSitesB,sheetB.col_values(plusCtsCol,1))}
minusCtsB = {tmpSite:int(tmpValue) for tmpSite,tmpValue in zip(insertionSitesB,sheetB.col_values(minusCtsCol,1))}



allInsertionSites = set(insertionSitesA).union(set(insertionSitesB))

#Find the mean value of the insertions 
meanPlusCts = {}
meanMinusCts = {}
for tmpSite in allInsertionSites:
	tmpPlusA = plusCtsA[tmpSite] if tmpSite in plusCtsA else 0
	tmpPlusB = plusCtsB[tmpSite] if tmpSite in minusCtsB else 0
	meanPlusCts[tmpSite] = 0.5*(tmpPlusA + tmpPlusB)
	tmpMinusA = minusCtsA[tmpSite] if tmpSite in plusCtsA else 0
	tmpMinusB = minusCtsB[tmpSite] if tmpSite in minusCtsB else 0
	meanMinusCts[tmpSite] = 0.5*(tmpMinusA + tmpMinusB)


#Load the gene locations from the aggregate sheet
colHeaders = [tmpCell.value for tmpCell in sheetAggr.row(0)]
locusCol = colHeaders.index(u'Locus')
startCol = colHeaders.index(u'Start')
endCol = colHeaders.index(u'End')
strandCol = colHeaders.index(u'Strand')

geneNames = sheetAggr.col_values(locusCol, 2)
geneInfo = {tmpGene: {'start':int(tmpStart), 'stop':int(tmpEnd), 'strand':tmpStrand} for \
					 tmpGene, tmpStart, tmpEnd, tmpStrand in zip(geneNames, \
					 sheetAggr.col_values(startCol,2), sheetAggr.col_values(endCol,2), \
					 sheetAggr.col_values(strandCol,2)) }


#Now go through and add arrows for each gene
figH = plt.figure()
axesH = figH.add_subplot(111)

totLength = np.max(sheetAggr.col_values(endCol,2))
arrowPatches = []
for tmpGeneName, tmpGeneInfo in geneInfo.items():
	if tmpGeneInfo['stop'] > ROI[0] and tmpGeneInfo['start'] < ROI[1] :
	    #Convert into figure coordinate
	    startFigCoord = tmpGeneInfo['start']
	    stopFigCoord = tmpGeneInfo['stop']
	    axesH.text(0.5*(stopFigCoord+startFigCoord), 0, tmpGeneName, horizontalalignment='center', verticalalignment='center')
	    if tmpGeneInfo['strand'] == '+':
	    	arrowPatches.append(mpatches.FancyArrow(startFigCoord, 0, stopFigCoord-startFigCoord, 0,
	    		width = 0.1, length_includes_head=True, head_length=0.1*(stopFigCoord-startFigCoord), edgecolor='none'))
	    else:
	    	arrowPatches.append(mpatches.FancyArrow(stopFigCoord, 0, startFigCoord-stopFigCoord, 0,
	    		width = 0.1, length_includes_head=True, head_length=0.1*(stopFigCoord-startFigCoord), edgecolor='none'))

    
arrowCollection = mcollections.PatchCollection(arrowPatches)
axesH.add_collection(arrowCollection)
plt.xlim((ROI[0],ROI[1]))

#Add the bars for insertions
lineCollectionPlus = mcollections.LineCollection([((tmpInsertion, 0),(tmpInsertion, np.log10(meanPlusCts[tmpInsertion]))) 
					for tmpInsertion in allInsertionSites if (meanPlusCts[tmpInsertion] > minPlotCts and tmpInsertion>ROI[0] and tmpInsertion<ROI[1]) ], 
					linewidths=2.0, color='g')
lineCollectionMinus = mcollections.LineCollection([((tmpInsertion, 0),(tmpInsertion, -np.log10(meanMinusCts[tmpInsertion]))) 
					for tmpInsertion in allInsertionSites if (meanMinusCts[tmpInsertion] > minPlotCts and tmpInsertion>ROI[0] and tmpInsertion<ROI[1]) ], 
					linewidths=2.0, color='r')

axesH.add_collection(lineCollectionPlus)
axesH.add_collection(lineCollectionMinus)
plt.ylim((-6,6))
plt.show()


