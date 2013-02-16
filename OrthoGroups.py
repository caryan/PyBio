"""
Module for OrthoMCL manipulations.
"""

import os,glob
from collections import Counter
import subprocess

from Bio import SeqIO
import pandas as pd

def faacleanup(filein,fileout):
	"""
	change the output of a prodigal faa file to simiplified locus identifiers
	"""
	#read in all records from original file
	with open(filein,'r') as FID:
		allRecords = list(SeqIO.parse(FID,'fasta'))
	print('{0} records found'.format(len(allRecords)))
	#replace id with ARCXXXX_OOO1 etc.
	for index, item in enumerate(allRecords):
        	newName = '{1}_{0:04d}'.format(index+1, item.id[0:7])
		item.id=newName
		item.name=''
		item.description=''
	#put out modified records to new file
	with open(fileout,'w') as FID:
		SeqIO.write(allRecords,FID,'fasta')

def find_orphans(fastaDir, mclOutput, orphanDir):
	"""
	Finds all genes that aren't in any orthogroups

	Parameters
	----------
	fastaDir= all the clean fasta files (directory)
	mcloutput= output from mcl (the file)
	orphanDir= where the output goes (directory)
	"""

	#Load the set of all orthogroups genes as a set
	with open(mclOutput,'r') as FID:
		orthoGenes = set(FID.read().replace('\n','\t').split('\t'))

	for fastaFile in os.listdir(fastaDir):
		with open(os.path.join(fastaDir, fastaFile), 'r') as FID:
			allGenes = SeqIO.to_dict(SeqIO.parse(FID, 'fasta'))

		geneIds = set(allGenes.keys())

		#Find the orphan genes as those not in the orthoGenes group
		orphanIds = geneIds.difference(orthoGenes)

		#Write them to file
		orphanGenes = [allGenes[id] for id in orphanIds]
		print('Percent orphans for {0} is {1:.2f}'.format(fastaFile, (100.0*len(orphanGenes))/len(geneIds)))
		with open(os.path.join(orphanDir, fastaFile), 'w') as FID:
			SeqIO.write(orphanGenes, FID, 'fasta')


def find_common_groups(groupFile, fastaDir, groupDir, coverageCutoff=1):
	"""
	Find all orthogroups that have coverage above a coverage cutoff
	groupFile= group output from orthoMCL
	fataDir= CleanFasta files directory
	groupDir= output directory -orthogroups in all genomes and pulls in sequence information asscoiated with prodigal id
	coverageCutoff= scaling factor for number of genomes present in
	"""

	#Load the orthogroups
	with open(groupFile, 'r') as FID:
		groupNames = []
		groupGenes = []
		groupTaxons = []
		allTaxons = set([])
		for line in FID:
			groupName, genes = line.split(':')
			genes = genes.split()
			taxons = [x.split('|')[0] for x in genes]
			taxonCounts = Counter(taxons)
			allTaxons = allTaxons.union(taxons)

			#If we have any paralogs i.e. two genes from same taxon then skip
			if max(taxonCounts.values()) > 1:
				continue
			else:
				groupNames.append(groupName)
				groupGenes.append(genes)
				groupTaxons.append(taxons)

	#Create the data frame
	orthoDF = pd.DataFrame({'genes':groupGenes, 'taxons':groupTaxons}, index=groupNames)

	#Calculate the coverage
	orthoDF['coverage'] = orthoDF['taxons'].apply(len)

	#Filter for those above coverage cutoff
	orthoDF = orthoDF[orthoDF['coverage'] >= coverageCutoff*len(allTaxons)]

	#Load all the genes from the fasta files 
	allGenes = {}
	for fastaFile in os.listdir(fastaDir):
		with open(os.path.join(fastaDir, fastaFile), 'r') as FID:
			allGenes.update(SeqIO.to_dict(SeqIO.parse(FID, 'fasta')))
		
	#Now loop through all the orthogroups and write files
	for groupName, groupGenes in orthoDF['genes'].iteritems():
		with open(os.path.join(groupDir, groupName+'.fasta'), 'w') as FID:
			geneRecords = [allGenes[gene] for gene in groupGenes]
			SeqIO.write(geneRecords, FID, 'fasta')

	return orthoDF

def align_group_genes(muscleBin, groupDir, alignedDir):
	"""
	Use Muscle to align genes from each group.
	muscleBin= mucsle binary full path to '/users/kzdv345/local/bin/muscle3.8.31_i86linux64'
	groupDir= directory made in the last step that contains all orthogroups shared by all genomes
	alignedDir= directory I make to store my output in
	"""
	for fastaFile in os.listdir(groupDir):
		groupName = fastaFile.split('.')[0]
		subprocess.call([muscleBin, '-in', os.path.join(groupDir, fastaFile), '-out', os.path.join(alignedDir, groupName+'.afa')])


def trimAl(trimBin, alignedDir, trimmedDir):
	"""
	Use trimAl to clean-up aligned sequences from Muscle.
	trimBin= trimAI binary full path to '/users/kzdv345/local/bin/trimal
	trimmedDir=directory where my trimmed output comes
	"""
	for fastaFile in os.listdir(alignedDir):
		groupName = fastaFile.split('.')[0]
		subprocess.call([trimBin, '-in', os.path.join(alignedDir, fastaFile), '-out', os.path.join(trimmedDir, groupName+'.afa'), '-htmlout', os.path.join(trimmedDir, groupName+'.html'), '-automated1'])

def concat_seqs(trimmedDir, outputFile):
	"""
	Pull out all the trimmed sequences associated with one genome and concatenate into one sequence.
	Then put all the sequences together in one file. 
	"""
	alignedFiles = glob.glob(trimmedDir+'/*.afa')

	#Use the first one to get the list of taxons
	strainNames = [x.split('|')[0] for x in SeqIO.to_dict(SeqIO.parse(open(alignedFiles[0], 'r'), 'fasta')).keys()]

	#Create a dictionary of empty strings to contain the concatenated sequences
	concatSeqs = {k:'' for k in strainNames}


	for fileName in alignedFiles:
		with open(fileName,'r') as FID:
			seqs = SeqIO.to_dict(SeqIO.parse(FID, 'fasta'))
			for strainName, record in seqs.items():
				concatSeqs[strainName.split('|')[0]] += record.seq.tostring()

	with open(outputFile, 'w') as FID:
		for strain in strainNames:
			FID.write('>{0}\n'.format(strain))
			FID.write(concatSeqs[strain])
			FID.write('\n')




