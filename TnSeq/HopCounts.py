# Copyright 2013 Colm Ryan colm@colmryan.org
# License GPL v3 (http://www.gnu.org/licenses/gpl.txt) 

'''
Script to go from fasta files to hopcounts and aggregate hop counts for TnSeq experiments.

Usage: python fasta2hopcount.py FastaFile IndexName AnnotationFile OutputFile

Author: Colm Ryan 

'''
# import sys

import subprocess
import argparse
import numpy as np
import pandas as pd

def fastq_filter(inFile):
	#Use the fastq tools to trim and filter the fasta file
	print('Clipping adapter string from reads ... ')
	# /opt/fastx/bin/fastx_clipper -i $forward -Q 33 -o $forward.clp -a CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC -l 28 -v >> $log 2>&1
	subprocess.call(['fastx_clipper', '-i', inFile, '-o', '{0}.clp'.format(inFile), '-a' , 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC', '-l28', '-Q33'] )
	print('Finished')

	print('Trimming sequences to predefined length ... ')
	# /opt/fastx/bin/fastx_trimmer -i $forward.clp -o $forward.clp.trm -l 25 -Q 33 >> $log 2>&1
	subprocess.call(['fastx_trimmer', '-i', '{0}.clp'.format(inFile), '-o', '{0}.clp.trm'.format(inFile), '-l25', '-Q33'])
	print('Finished')

	print('Filtering for quality... ')
	#/opt/fastx/bin/fastq_quality_filter -i <fastq-file> -o <fastq-file.qc> -q 7 -p 95 -Q 33 -v
	# /opt/fastx/bin/fastq_quality_filter -i $forward.clp.trm -o $forward.clp.trm.qc -q7 -p 95 -Q 33 -v >> $log 2>&1
	subprocess.call(['fastq_quality_filter', '-i', '{0}.clp.trm'.format(inFile), '-o', '{0}.clp.trm.qc'.format(inFile), '-q7' , '-p95', '-Q33' ])
	print('Finished')

def run_bowtie(inFile, indexFile, outFile):
	#Use the arguments from the galaxy server
	#bowtie -q -p 4 -S -n 2 -e 70 -l 28 --maxbts 800 -k 1 -m 1 --best --phred33-quals S.aureus Sample_MHcontrol2.R1.fastq.clp.trm.qc > Sample_MHcontrol2.R1.fastq.sam
	print('Mapping with bowtie....')
	with open(outFile,'w') as outFID:
		subprocess.call(['bowtie', '-q', '-p4', '-S', '-n2', '-e70', '-l28', '--maxbts', '800', '-k1', '-m1', '--best', '--phred33-quals', indexFile, inFile], stdout=outFID)
	print('Finshed')

def hopcount_calc(samFile, outFile, annotationFile, refGenomeFile):
	#Use cut to reduce the SAM to the essential columns
	smallSAM = subprocess.Popen(['cut', '-f2,4', samFile], stdout=subprocess.PIPE)

	#Load the SAM data into a data frame
	try:
	    from cStringIO import StringIO
	except:
	    from StringIO import StringIO
	SAMData = pd.read_table(StringIO(smallSAM.communicate()[0]), skiprows=3, names=['direction', 'location'])

	#Filter for plus and minus directions
	plusData = SAMData[SAMData['direction'] == 0]
	minusData = SAMData[SAMData['direction'] == 16]
	#For some reason, minus data locations are shifted downy by 24 pts
	minusData['location'] += 24
	# del SAMData

	hopCounts = pd.DataFrame(plusData['location'].value_counts(), columns=['PlusCount'])
	hopCounts = hopCounts.join(pd.DataFrame(minusData['location'].value_counts(), columns=['MinusCount']), how='outer')

	#Clean-up
	hopCounts = hopCounts.fillna(0)
	hopCounts['PlusCount'] = np.int64(hopCounts['PlusCount'])
	hopCounts['MinusCount'] = np.int64(hopCounts['MinusCount'])

	#Add the total counts as the sum
	hopCounts['TotalCount'] = hopCounts['PlusCount']+hopCounts['MinusCount']

	#Add in the DNA strings from the reference file
	with open(refGenomeFile, 'r') as FID:
		FID.readline()
		referenceGenome = FID.read().replace('\n', '')
	#Positions seem to be off by one??? 
	hopCounts['Sequence'] = [referenceGenome[tmpPos-1:tmpPos+24] for tmpPos in hopCounts.index]

	#Add the extra columns of information for the genes
	hopCounts['Locus'] = ''
	hopCounts['ProteinID'] = ''
	hopCounts['Notation'] = ''
	hopCounts['Product'] = ''

	#Load the annotation information
	annotData = pd.read_table(annotationFile, index_col=3, names=['chr', 'Start', 'End', '', 'Gene', 'Notation', 'Product', 'ProteinID', 'junk'])
	del annotData['chr']
	del annotData['junk']
	
	for tmpLocus, tmpAnnot in annotData.iterrows():
		goodLocations = (hopCounts.index>tmpAnnot['Start']) & (hopCounts.index<tmpAnnot['End'])
		if np.any(goodLocations):
			hopCounts.ix[goodLocations, ['Locus']] = tmpLocus
			hopCounts.ix[goodLocations, ['ProteinID']] = tmpAnnot['ProteinID']
			hopCounts.ix[goodLocations, ['Notation']] = tmpAnnot['Notation']
			hopCounts.ix[goodLocations, ['Product']] = tmpAnnot['Product']

	hopCounts.to_csv('{0}.hopCount'.format(outFile), sep='\t', index_label='Position')

	#Now aggregate on genes
	presentAggHopCounts = hopCounts.groupby('Locus').agg({'MinusCount':'sum', 'PlusCount':'sum', 'TotalCount':'sum',
													'ProteinID':'first', 'Notation':'first', 'Sequence':'first'})
	presentAggHopCounts['Sites'] = hopCounts['Locus'].value_counts()

	allAggHopCounts = annotData.copy()
	allAggHopCounts['MinusCount'] = 0
	allAggHopCounts['PlusCount'] = 0
	allAggHopCounts['TotalCount'] = 0
	allAggHopCounts['Sites'] = 0
	allAggHopCounts['Length'] = allAggHopCounts['End'] - allAggHopCounts['Start']


	allAggHopCounts.update(presentAggHopCounts)

	#Add the DVal's
	totalCount = hopCounts['TotalCount'].sum()
	totalLength = len(referenceGenome)
	allAggHopCounts['ExpectedCount'] = allAggHopCounts['Length']*totalCount/totalLength
	allAggHopCounts['DVal'] = allAggHopCounts['TotalCount']/allAggHopCounts['ExpectedCount']

	#Make it pretty
	allAggHopCounts['PlusCount'] = np.int64(allAggHopCounts['PlusCount'])
	allAggHopCounts['MinusCount'] = np.int64(allAggHopCounts['MinusCount'])
	allAggHopCounts['TotalCount'] = np.int64(allAggHopCounts['TotalCount'])
	allAggHopCounts['Start'] = np.int64(allAggHopCounts['Start'])
	allAggHopCounts['End'] = np.int64(allAggHopCounts['End'])
	allAggHopCounts['Sites'] = np.int64(allAggHopCounts['Sites'])

	allAggHopCounts = allAggHopCounts.reindex(columns=['Sites', 'PlusCount', 'MinusCount', 'TotalCount', 'Start', 'End', 'Length', 'ExpectedCount', 'DVal', 'Product', 'ProteinID', 'Notation'])

	allAggHopCounts.to_csv('{0}.aggHopCount'.format(outFile), sep='\t', index_label='Locus')

def run_pipeline(fastaFile, indexName, refGenomeFile, annotationFile, outputFile):
	"""
	Helper function that runs pipeline by calling all three functions. 
	"""

	#Call the fastq clip, trim, filter set
	fastq_filter(fastaFile)

	#Call bowtie to do the mapping to the reference genomre
	run_bowtie('{0}.clp.trm.qc'.format(fastaFile), indexName, '{0}.sam'.format(fastaFile))

	#Call the hop count helper 
	hopcount_calc('{0}.sam'.format(fastaFile), outputFile, annotationFile, refGenomeFile)


if __name__ == '__main__':
	#Pull out the arguments.  

	parser = argparse.ArgumentParser()
	parser.add_argument('fastaFile', help='location of input fasta file')
	parser.add_argument('indexName', help='name of the bowtie index')
	parser.add_argument('refGenomeFile', help='fasta file for the reference genome')
	parser.add_argument('annotationFile', help='location of annotation file describing gene locations')
	parser.add_argument('outputFile', help='what to call the output hopcount files')
	args = parser.parse_args()

	run_pipeline(args.fastaFile, args.indexName, args.refGenomeFile, args.annotationFile, args.outputFile)

