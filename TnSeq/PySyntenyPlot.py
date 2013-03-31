# Copyright 2013 Colm Ryan colm@colmryan.org
# License GPL v3 (http://www.gnu.org/licenses/gpl.txt) 

'''
Created on Jul 5, 2011

Do some arrow synteny plots for Veronica

@author: caryan
'''

from __future__ import division
 
#Let's write to SVG style graphics
#import matplotlib
#matplotlib.use('svg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.collections as mcollections

import numpy as np

import csv

#Define a helper function to conver to figure coordinates.
def calcFigCoords(x, leftLimit, rightLimit):
    if leftLimit < rightLimit:
        return (x-leftLimit)/(rightLimit-leftLimit)
    else:
        return 1-(x-rightLimit)/(leftLimit-rightLimit)
    
def calcTickLines(leftLimit,rightLimit,step,vertShift,tickHeight):
    if leftLimit < rightLimit:
        tickHPos = np.arange(leftLimit,rightLimit-(step-1),step,dtype=np.float64)
        tickHPos = np.hstack((tickHPos,rightLimit))
        #Scale appropriately to 0..1
        tickHPos -= leftLimit
        tickHPos /= (rightLimit-leftLimit)
    else:
        tickHPos = np.arange(rightLimit,leftLimit-(step-1),step,dtype=np.float64)
        tickHPos = np.hstack((tickHPos,leftLimit))
        #Scale appropriately to 0..1
        tickHPos -= rightLimit
        tickHPos /= (leftLimit-rightLimit)
        tickHPos = 1-tickHPos
       
    tickLineSegs = []
    for tmpHPos in tickHPos:
        tickLineSegs.append(np.array([[tmpHPos, vertShift],[tmpHPos, vertShift+tickHeight]])) 

    return tickLineSegs
    

def createTickPatches(leftLimit, rightLimit, vertShift):
    #The tall ticks occur every 2000
    tallTickSegs = calcTickLines(leftLimit, rightLimit, 2000, vertShift, 0.1)

    #The medium ticks occur every 1000
    medTickSegs = calcTickLines(leftLimit, rightLimit, 1000, vertShift, 0.05)

    #The small ticks occur every 100
    smallTickSegs = calcTickLines(leftLimit, rightLimit, 200, vertShift, 0.02)

    totalTicks = []
    totalTicks.extend(tallTickSegs)
    totalTicks.extend(medTickSegs)
    totalTicks.extend(smallTickSegs)
    
    tickCollection = mcollections.LineCollection(totalTicks)
    tickCollection.set_color('k')
    tickCollection.set_linewidth(2)
    return tickCollection
    
if __name__ == '__main__':
    
    #Load the data from the gatt files
    gattFileName = 'combined4.gatt'
    with open(gattFileName,'r') as tmpFID:
        tmpReader = csv.reader(tmpFID, delimiter='\t')
        gattFileInfo = []
        for row in tmpReader:
            tmpInfo = {}
            tmpInfo['start'] = int(row[0])
            tmpInfo['stop'] = int(row[1])
            tmpInfo['name'] = row[3]
            tmpInfo['fragmentName'] = row[-1]
            gattFileInfo.append(tmpInfo)

    #Load the info from the fragment file
    fragFileName = 'van_operon_fragment_file.txt'
    with open(fragFileName,'r') as tmpFID:
        tmpReader = csv.reader(tmpFID, delimiter='\t')
        fragFileInfo = {}
        for row in tmpReader:
            tmpName = row[1]
            fragFileInfo[tmpName] = {}
            fragFileInfo[tmpName]['leftEnd'] = int(row[2])
            fragFileInfo[tmpName]['rightEnd'] = int(row[3])
        
    #Create the tick marker patches
    tickCollection = createTickPatches(fragFileInfo['VRSA11a_contig00074']['leftEnd'], fragFileInfo['VRSA11a_contig00074']['rightEnd'],0.5)
    
    #Now go through and add arrows for each gene
    figH = plt.figure()
    axesH = figH.add_subplot(111)

    arrowPatches = []
    for tmpGene in gattFileInfo:
        if tmpGene['fragmentName'] == 'VRSA11a_contig00074':
            tmpFragInfo = fragFileInfo[tmpGene['fragmentName']]
            #Convert into figure coordinate
            startFigCoord = calcFigCoords(tmpGene['start'], tmpFragInfo['leftEnd'], tmpFragInfo['rightEnd'])
            stopFigCoord = calcFigCoords(tmpGene['stop'], tmpFragInfo['leftEnd'], tmpFragInfo['rightEnd'])
            arrowPatches.append(mpatches.Arrow(startFigCoord, 0.4, stopFigCoord-startFigCoord,0, width=0.3))
            
    arrowCollection = mcollections.PatchCollection(arrowPatches)
    axesH.add_collection(arrowCollection)
    axesH.add_collection(tickCollection)
#    axesH.autoscale_view(True, True, True)
    plt.show()
    
            
        
              
    
    
    