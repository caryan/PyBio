# Copyright 2013 Colm Ryan colm@colmryan.org
# License GPL v3 (http://www.gnu.org/licenses/gpl.txt) 

"""
Module for OrthoMCL manipulations.
"""

import os,glob
from collections import Counter
import subprocess
import csv
import re
import difflib

from Bio import SeqIO
import pandas as pd


def pug_fixer(filein,fileout):
	"""
	changes the nucleotide file (fnr) file from pubmed to have locus tags at the start of the name of each gene
	"""
	#read in records from original files
	with open(filein,'r') as FID:
		allRecords=list(SeqIO.parse(FID,'fasta'))
	pat = re.compile('locus_tag=(\w+)')
	for rec in allRecords:
		geneName = re.findall(pat, rec.description)[0]
		rec.id = geneName
		rec.name = ''
		rec.description = ''

	#put out modified records to new file
	with open(fileout,'w') as FID:
		SeqIO.write(allRecords,FID,'fasta')


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

def pull_seqs(fastaDir, genes):
	"""
	Pull gene info out of fasta files given a list of gene names.

	Parameters
	------------
	fastaDir : directory of fasta files
	genes : iterable of gene names

	Output
	--------------
	list of Biopython sequence records

	"""
	
	#Load all the genes from the all fasta files 
	#TODO: only load the necessary files
	allGenes = {}
	for fastaFile in os.listdir(fastaDir):
		with open(os.path.join(fastaDir, fastaFile), 'r') as FID:
			allGenes.update(SeqIO.to_dict(SeqIO.parse(FID, 'fasta')))
	
	#Return the relevant records
	return [allGenes[gene] for gene in genes]

def find_common_groups(groupFile, fastaDir, groupDir, coverageCutoff=1):
	"""
	Find all orthogroups that have coverage above a coverage cutoff

	Parameters
	----------
	groupFile= group output from orthoMCL
	fastaDir= CleanFasta files directory
	groupDir= output directory -orthogroups in all genomes and pulls in sequence information asscoiated with prodigal id
	coverageCutoff= scaling factor for number of genomes present in
	"""

	#Load the orthogroups
	with open(groupFile, 'r') as FID:
		groupNames = []
		groupGenes = []
		groupStrains = []
		allStrains = set([])
		for line in FID:
			groupName, genes = line.split(':')
			genes = genes.split()
			strains = [x.split('|')[0] for x in genes]
			strainCounts = Counter(strains)
			allStrains = allStrains.union(strains)

			#If we have any paralogs i.e. two genes from same strain then skip
			if max(strainCounts.values()) > 1:
				continue
			else:
				groupNames.append(groupName)
				groupGenes.append(genes)
				groupStrains.append(strains)

	#Create the data frame
	orthoDF = pd.DataFrame({'genes':groupGenes, 'strains':groupStrains}, index=groupNames)

	#Calculate the coverage
	orthoDF['coverage'] = orthoDF['strains'].apply(len)

	#Filter for those above coverage cutoff
	orthoDF = orthoDF[orthoDF['coverage'] >= coverageCutoff*len(allStrains)]

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

	Parameters
	----------
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

	Parameters
	-----------
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

	#Use the first one to get the list of strains
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

def write_group_table(groupFile, cleanFastaDir, outputFile=None):
	"""
	Write out the orthogroup information in tabular format amenable to filtering in Excel.

	Basically we unfold the orthogroups and for reference genes we also print out info

	"""

	# #Dictionary of reference genomes keyed on taxon in group file and values of genome files.
	# refGenomeFiles = {'PAO1':'OriginalFastaFiles/PAO1PG.faa',
	#               'LES':'OriginalFastaFiles/LESB58PG.faa'}


	# #Load the reference genomes as a dictionaries keyed off gene id
	# refGenomes = {}
	# for shortName, origFastaFile in refGenomeFiles.items():
	# 	with open(origFastaFile, 'r') as FID:
	# 		refGenomes[shortName] = SeqIO.to_dict(SeqIO.parse(FID, 'fasta'), key_function=lambda rec : rec.id.split('|')[0])


	# #Pattern to remove extraneous info from description field
	# pat = re.compile('\[.*\]')

	# #Go through the orthogroup info
	# outWriter = open(outputFile,'w')
	# with open(groupFile, 'r') as inFID, open(outputFile,'w') as outFID:
	# 	outWriter = csv.writer(outFID, delimiter='\t')
	# 	for line in inFID:
	# 		groupName, genes = line.split(':')
	# 		genes = genes.split()
	# 		for gene in genes:
	# 			shortName, geneName = gene.split('|')
	# 			#If it is the reference list then give extra descritpion info.
	# 			if shortName in refGenomes.keys():
	# 				#Clean up description info
	# 				geneInfo = pat.sub('', refGenomes[shortName][geneName].description.replace(geneName+'|', ''))
	# 			else:
	# 				geneInfo = ''
	# 			outWriter.writerow([groupName, geneName, geneInfo ])

	#Use the clean fasta file directory to list all genomes
	genomeNames = [x.split('.')[0] for x in os.listdir(cleanFastaDir)]
	genes = {name:[] for name in genomeNames}

	#Go through the orthogroup info
	groupNames = []
	with open(groupFile,'r') as FID:
		for line in FID:
			groupName, groupGenes = line.split(':')
			groupGenes = groupGenes.split()
			strains = [x.split('|')[0] for x in groupGenes]
			strainCounts = Counter(strains)
			#If we have any paralogs i.e. two genes from same strain then skip
			if max(strainCounts.values()) > 1:
				continue

			groupNames.append(groupName)
			groupGenes = {k:v for k,v in [g.split('|') for g in groupGenes]}
			for name, geneGroups in genes.items():
				if name in groupGenes:
					geneGroups.append(groupGenes[name])
				else:
					geneGroups.append(None)

	df = pd.DataFrame(genes, index=groupNames)
	
	if outputFile:
		df.to_csv(outputFile, sep='\t')

	return df	


def trim_seqs(alignedDir, trimmedDir):
	"""
	A simple trimming program.  The claim is that with modern sequencing trimming in the middle of the 
	sequence is not as important and we should simply trim the beginning and ends.

	Parameters
	----------
	alignedDir : directory of input aligned fasta files
	trimmedDir : directory of output trimmed fasta files
	"""
	#Some regular expressions for matching '-'s at the beginning and end of strings
	frontPat = re.compile('^-*')
	endPat = re.compile('-*$')

	for fastaFile in os.listdir(alignedDir):
		groupName = fastaFile.split('.')[0]
		#Load all the sequences
		with open(os.path.join(alignedDir, fastaFile),'r') as FID:
			seqs = list(SeqIO.parse(FID, 'fasta'))
			#Now look for the most "-"s at the beginning or end of the file
			frontTrim = 0
			endTrim = 0
			for seq in seqs:
				frontTrim = max(frontTrim, re.search(frontPat, seq.seq.tostring()).end())
				endTrim = max(endTrim, len(seq.seq) - re.search(endPat, seq.seq.tostring()).start())

			#Re-loop and trim
			for seq in seqs:
				if frontTrim:
					seq.seq = seq.seq[frontTrim:]
				if endTrim:
					seq.seq = seq.seq[:-endTrim]

		#Write out to file
		with open(os.path.join(trimmedDir, fastaFile), 'w') as FID:
			SeqIO.write(seqs, FID, 'fasta')

def create_orthoSets(df):
	"""
	Create a dictionary of sets of orthogroups keyed off strains.

	Parameters
	-----------
	df : DataFrame of orthogroup membership
	"""
	orthoSets = {}
	for strain in df.columns:
		orthoSets[strain] = set(df[strain][pd.notnull(df[strain])].index)

	return orthoSets

def filter_orthoSets(filterList, orthoSets):
	"""
	Filter the orthogroups to find the exclusive members of a list of strains. 

	Parameters
	----------
	filterList : iterable of strains to filter for exclusive membership in orthogroups
	orthoSets : dictionary of sets of orthogroups for each strain
	"""
	#Take the union of the sets in the filterList
	inGroup = set([])
	outGroup = set([])

	for strain,strainSet in orthoSets.items():
		if strain in filterList:
			if inGroup:
				inGroup &= strainSet
			else:
				inGroup = strainSet
		else:
			outGroup |= strainSet

	return inGroup - outGroup

def find_aa_changes(inputDir, refStrain, outputDir):
	"""
	Enumerate all the amino acid changes in an orthogroup vs a reference gene.


	Parameters
	-------------
	inputDir : directory containing aligned orthogroup files
	refStrain : the reference strain to find differences with
	outputDir : directory to put output files in  
	"""

	#Iterate over each aligned fasta file in the directory
	for fastaFile in os.listdir(inputDir):
		groupName = fastaFile.split('.')[0]
		#Load all the sequences
		with open(os.path.join(inputDir, fastaFile),'r') as FID:
			seqs = SeqIO.to_dict(SeqIO.parse(FID, 'fasta'))

		#Find the reference strain gene
		refStrainGene = filter(lambda name : name.startswith(refStrain), seqs.keys())
		#If reference is theres
		if refStrainGene:
			refGeneSeq = seqs[refStrainGene[0]].seq.tostring()
			#Open the output file
			with open(os.path.join(outputDir, groupName+'_Changes.tsv'), 'w') as FID:
				#Write a header
				FID.write('\t'.join(['Gene', 'Start Idx.', 'Stop Idx.', 'Ref String', 'Replacement String']))
				FID.write('\n')
				#Iterate over strains
				for strainName, seq in seqs.items():
					if strainName != refStrainGene[0]:
						startdIdx = 0
						curIdx = 0
						compSeq = seq.seq.tostring()

						while curIdx < len(refGeneSeq):
							if refGeneSeq[curIdx] != compSeq[curIdx]:
								startIdx = curIdx
								curIdx += 1
								while curIdx < len(refGeneSeq):
									if refGeneSeq[curIdx] == compSeq[curIdx]:
										break
									curIdx += 1 
								FID.write('\t'.join([ strainName, str(startIdx+1), str(curIdx), refGeneSeq[startIdx:curIdx], compSeq[startIdx:curIdx] ]))
								FID.write('\n')
							curIdx += 1

def consolidate_aa_change_info(inputDir, drugInfoFile, drug, outputDir, totNumStrains):
	"""
	Consoldiate the ammino acid change information by grouping each change and listing involved strains and MIC info.

	Parameters
	------------
	inputDir : directory containing tsv files from find_aa_changes
	drugInfoFile : csv file containing MIC info for each strain and drug
	drug : which drug to print MIC info for
	outputDir : where to put the files
	totNumStrains : total number of strains for relative %'s
	"""

	#Load the drug info
	drugInfoDF = pd.read_csv(drugInfoFile, index_col='strain')
	
	#Loop over input files
	for inputFile in os.listdir(inputDir):
		aaChangeDF = pd.read_csv(os.path.join(inputDir, inputFile), sep='\t')

		#Group by the changes
		grouped = aaChangeDF.groupby(['Start Idx.', 'Ref String', 'Replacement String'])

		#Write the results out to file
		geneName = inputFile.split('_')[0]
		with open(os.path.join(outputDir,geneName), 'w') as FID:
			for change, strainNums in grouped.groups.items():
				FID.write(' '.join([str(x) for x in change]))
				FID.write(': {0:.1%} of strains\n'.format(float(len(strainNums))/totNumStrains))
				for strainNum in strainNums:
					strain = aaChangeDF['Gene'][strainNum].split('|')[0]
					try:
						MIC = drugInfoDF[drug][strain]
						country = drugInfoDF['country'][strain]
					except KeyError:
						MIC = 'no drug info.'
						country = ''
					FID.write('\t{0} {1} ; {2}\n'.format(strain, MIC, country))













