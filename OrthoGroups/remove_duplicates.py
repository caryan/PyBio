"""
Remove duplicate strains from an aligned fasta file.  Keep the one that best matches a reference genome at the N terminus. 
"""
from Bio import SeqIO
from collections import Counter
import numpy as np

def remove_duplicates(fastaFile):

	#Load the original aligned fasta file
	with open(fastaFile, 'r') as FID:
		curRecords = SeqIO.to_dict(SeqIO.parse(FID, 'fasta'))
		strains = [_.split('|')[0] for _ in curRecords.keys()]
		strainCounts = Counter(strains)
		for strain, ct in strainCounts.items():
			if ct > 1:
				#Find gene names for this strain
				badGenes = filter(lambda g : g.startswith(strain), curRecords.keys())
				print('Strain {} has more than one entry: {}'.format(strain, '; '.join(badGenes)))
				#Assume that the one that starts with the least dashes is the best
				dashCounts = [len(curRecords[gene].seq) - len(curRecords[gene].seq.tostring().lstrip('-')) for gene in badGenes]
				bestOfBad = badGenes.pop(np.argmin(dashCounts))
				print('Choosing {}'.format(bestOfBad))
				#Remove badGenes for records dictionary
				for gene in badGenes:
				  curRecords.pop(gene)
	#Write the record dictionary back to file
	with open(fastaFile,'w') as FID:
		SeqIO.write(curRecords.values(), FID, 'fasta')
