# Copyright 2013 Colm Ryan colm@colmryan.org
# License GPL v3 (http://www.gnu.org/licenses/gpl.txt) 

"""
Plotting synteny functions. 

Requires: matplotlib; pandas; Biopython

"""

import pandas as pd
import numpy as np
import re
import subprocess
from cStringIO import StringIO


import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.collections as mcollections

def plot_synteny(geneInfo, contigs, ROI):
	"""
	Helper function to plot out gene information.  

	Parameters
	----------
	geneInfo : dictionary with strainName : dataframe of gene start/stops (use parse_prodigal)
	contigs : dictionary with strainName : which contig to plot
	ROI : tuple of (start, stop) : start/stop index to plot 
	"""
	figH = plt.figure()
	axesH = figH.add_subplot(111)

	pickedPatch = []

	def get_patch_picker(pickedPatch):
		def patch_picker(event):
			if pickedPatch:
				pickedPatch[0].set_facecolor('b')
				pickedPatch.pop()
			pickedPatch.append(event.artist)
			pickedPatch[0].set_facecolor('r')
			event.canvas.draw()
			vertices = event.artist.get_path().vertices
			vertShift = vertices[0,1]
			left = int(min(vertices[:,0]))
			right = int(max(vertices[:,0]))

			print('Left: {0:d}; Right {1:d}'.format(left, right))
			return True
		return patch_picker

	vertShift = 0

	for strain, geneDF in geneInfo.items():
		vertShift += 1
		#Pull out the good contig and region of interest
		geneDF = geneDF[(geneDF['Contig'] == contigs[strain]) & (geneDF['StartIdx'] < ROI[1]) & (geneDF['StopIdx'] > ROI[0])]
		plt.text(geneDF['StopIdx'].values[0], vertShift-0.125, strain, fontsize=16)
		for start, stop, strand in zip(geneDF['StartIdx'], geneDF['StopIdx'], geneDF['Strand']):
			if strand == '+':
				axesH.add_patch(mpatches.FancyArrow(start, vertShift, stop-start, 0,
					width = 0.1, length_includes_head=True, head_length=0.1*(stop-start), edgecolor='none', picker=5))
			else:
				axesH.add_patch(mpatches.FancyArrow(stop, vertShift, start-stop, 0,
					width = 0.1, length_includes_head=True, head_length=0.1*(stop-start), edgecolor='none', picker=5))
		axesH.add_collection(create_ticks(ROI, vertShift-0.25))
	

	plt.xlim(ROI[0], ROI[1])
	plt.ylim(0.5, vertShift+0.5)
	figH.canvas.mpl_connect('pick_event', get_patch_picker(pickedPatch))
	plt.show()

def parse_prodigal(fileName):
	"""
	Helper function to parse prodigal output and extract gene names and start/stop locations
	"""
	#We'll parse the GFF files from Prodigal
	#Pandas cannot yet deal with line comments so use sed to strip those
	cleanFile = subprocess.check_output(['sed', '/^#/d', fileName])

	#We can extract the contig number from the first column
	contigNumRE = re.compile('contig_(\d+)')
	contig_extractor = lambda s : int(contigNumRE.search(s).group(1))

	return pd.read_csv(StringIO(cleanFile), sep='\t', usecols=[0,3,4, 6], names=['Contig', 'StartIdx', 'StopIdx', 'Strand'], 
						converters={0: contig_extractor}, )


def create_ticks(ROI, vertPos):
	"""
	Helper function to create a line collection of ticks 
	"""
	#The tall ticks occur every 2000
	tallTicks = [np.array([[hPos, vertPos], [hPos, vertPos+0.1]]) for hPos in np.arange(ROI[0], ROI[1], 2000)]

	#The medium ticks occur every 1000
	mediumTicks = [np.array([[hPos, vertPos], [hPos, vertPos+0.05]]) for hPos in np.arange(ROI[0], ROI[1], 1000)]

	#The small ticks occur every 500
	smallTicks = [np.array([[hPos, vertPos], [hPos, vertPos+0.02]]) for hPos in np.arange(ROI[0], ROI[1], 500)]

	allTicks = []
	allTicks.extend(tallTicks)
	allTicks.extend(mediumTicks)
	allTicks.extend(smallTicks)

	tickCollection = mcollections.LineCollection(allTicks)
	tickCollection.set_color('k')
	tickCollection.set_linewidth(2)
	return tickCollection
	