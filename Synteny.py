# Copyright 2013 Colm Ryan colm@colmryan.org
# License GPL v3 (http://www.gnu.org/licenses/gpl.txt) 

"""
Plotting synteny functions. 

Requires: matplotlib; pandas; Biopython

"""

import pandas as pd
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

	vertShift = 0
	arrowPatches = []
	for strain, geneDF in geneInfo.items():
		vertShift += 1
		#Pull out the good contig and region of interest
		geneDF = geneDF[(geneDF['Contig'] == contigs[strain]) & (geneDF['StartIdx'] < ROI[1]) & (geneDF['StopIdx'] > ROI[0])]
		plt.text(geneDF['StopIdx'].values[0], vertShift-0.25, strain, fontsize=16)
		for start, stop, strand in zip(geneDF['StartIdx'], geneDF['StopIdx'], geneDF['Strand']):
			if strand == '+':
				arrowPatches.append(mpatches.FancyArrow(start, vertShift, stop-start, 0,
					width = 0.1, length_includes_head=True, head_length=0.1*(stop-start), edgecolor='none'))
			else:
				arrowPatches.append(mpatches.FancyArrow(stop, vertShift, start-stop, 0,
					width = 0.1, length_includes_head=True, head_length=0.1*(stop-start), edgecolor='none'))


	arrowCollection = mcollections.PatchCollection(arrowPatches)
	axesH.add_collection(arrowCollection)
	plt.xlim(ROI[0], ROI[1])
	plt.ylim(0.5, vertShift+0.5)
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


