from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet
import pandas as pd
from progressbar import Percentage, Bar, ProgressBar, RotatingMarker, ETA
import glob
import cStringIO
import subprocess

# Wrangle multiple OrthoMCL results
import os.path

#path to folder containing all the other folders
folderPath = '/home/cryan/Downloads/resistome_set'

#Load the set of genes from PAO1 that we are interested in
with open(os.path.join(folderPath, 'porin_test.txt'),'r') as FID:
	porins = FID.read().splitlines()

#Load the set of strains that we are interested in
with open(os.path.join(folderPath, 'strains_nov_2013.txt'),'r') as FID:
	allStrains = FID.read().splitlines()

#Load the orthogroups
fileNames = glob.glob(os.path.join(folderPath, 'group_output', 'groups*.txt'))
OGs = []
for fileName in fileNames:
	with open(fileName,'r') as FID:
		OGs.append({})
		for line in FID:
			groupName, genes = line.split(':')
			OGs[-1][groupName] = set(genes.split())


#Load all the fasta files we'll need 
#This will take a while so throw up a progressbar
widgets = ['Loading fasta files: ', Percentage(), ' ', Bar(marker=RotatingMarker()),' ', ETA(), ' ']
pbar = ProgressBar(widgets=widgets, maxval=len(allStrains))
pbar.start()
seqs = {}
for ct, strainName in enumerate(allStrains):
	try:
		with open(os.path.join(folderPath, 'fastafiles', strainName + '.fasta'), 'r') as FID:
			seqs[strainName] = SeqIO.to_dict(SeqIO.parse(FID, 'fasta'))
			#Strip of '*' stop codon to avoid warning from muscle
			for rec in seqs[strainName].values():
				rec.seq = rec.seq.rstrip('*')

	except IOError:
		print("Unable to find fasta file for {}".format(strainName))
	pbar.update(ct)

goodStrains = set(seqs.keys())

#Now loop through the porin list; put together the sequences that contain it and write a new fasta file
for porin in porins:
	geneSet = set([])
	#Loop through the OGs and find ones where this porin is present
	for ogs in OGs:
		for genes in ogs.values():
			if 'PAO1|' + porin in genes:
				geneSet.update(genes)
				break

	if not geneSet:
		print('Failed to find any orthogroups for {}, not writing a fasta file.'.format(porin))
		continue

	#Filter out genes that are from strains in the strain list
	geneSet = set(filter(lambda g : g.split('|')[0] in goodStrains, geneSet))

	#Now we have all the genes go back to the sequences and put together a fasta file
	strainNames = set([])
	seqRecords = []
	maxLength = 0
	for gene in geneSet:
		strainName = gene.split('|')[0]
		strainNames.update([strainName])
		seqRecords.append(seqs[strainName][gene])
		maxLength = max(maxLength, len(seqRecords[-1]))

	#Add in dummy dash sequences for strains not represented
	# fillIns = goodStrains - strainNames
	# for fillIn in fillIns:
	# 	seqRecords.append(SeqRecord(Seq('-'*maxLength, SingleLetterAlphabet), id=fillIn+'|DummySeq'))

	fastaStr = cStringIO.StringIO()
	SeqIO.write(seqRecords, fastaStr, 'fasta')

	#Finally call muscle to align the sequences
	procHandle = subprocess.Popen(['muscle', '-quiet', '-diags', '-out', os.path.join(folderPath, 'aligned', porin+".afa")], stdin=subprocess.PIPE)
	procHandle.communicate(fastaStr.getvalue())
	procHandle.wait()



