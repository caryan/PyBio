import glob
from Bio import SeqIO

import pandas as pd

#Use the MIC list as the master list of strains
# with open("strains_dec_7_2013.csv", 'r') as FID:
# 	strains = [s.rstrip() for s in FID.readlines()]
# MICdf = pd.read_csv('antibiogram_july.csv', index_col=0)
# strains = MICdf.index.values
df = pd.read_csv("../meta_clc2.csv", index_col=0)
strains = ["ARC"+str(x) for x in df.index.values]


fileNames = glob.glob('*.afa')

get_strain_name = lambda seq: seq.id.split('|')[0]

porinCheck = pd.DataFrame(index=strains)

for fileName in fileNames:
	seqs = SeqIO.to_dict(SeqIO.parse(fileName, 'fasta'), key_function=get_strain_name)

	#Find the reference record first
	refName = 'PAO1'
	refSeq = seqs[refName].seq.tostring()
	#Figure out the length of the refernece sequence by stripping off leading/trailing '---'
	refSeqLength = len(refSeq.strip('-'))
	refSeqStart = len(refSeq) - len(refSeq.lstrip('-'))
	refSeqStop = len(refSeq.rstrip('-'))

	#Compare each to the reference strain
	#Record P if it there; A if it is absent and T for truncation
	groupName = fileName.split('.')[0]
	strainCheck = []
	for strain in strains:
		#If the strain is not here then note and move on
		if strain not in seqs.keys():
			strainCheck.append('A')
		#Otherwise, check for truncation by comparing to length of refSeq
		else:
			seqString = seqs[strain].seq.tostring()
			seqStart = len(seqString) - len(seqString.lstrip('-'))
			seqStop = len(seqString.rstrip('-'))
			myString = ''
			if (seqStart == refSeqStart) and (seqStop == refSeqStop):
				myString += 'P; '
			else:
				if seqStart < refSeqStart:
					myString += 'E_N{}; '.format(refSeqStart-seqStart)
				elif seqStart > refSeqStart:
					myString += 'T_N{}; '.format(seqStart-refSeqStart)
				if seqStop > refSeqStop:
					myString += 'E_C{}; '.format(seqStop-refSeqStop)
				elif seqStop < refSeqStop:
					myString += 'T_C{}; '.format(refSeqStop-seqStop)

			#Now check for changes within the sequences too
			start = max(refSeqStart, seqStart)
			stop = min(refSeqStop , seqStop)
			for ct, (refChar, seqChar) in enumerate(zip(refSeq[start:stop], seqString[start:stop])):
				if refChar == '-' and seqChar != '-':
					myString += "D-@{}; ".format(ct)
				elif seqChar == '-' and refChar != '-':
					myString += "D+@{}; ".format(ct)

			strainCheck.append(myString)
  	porinCheck[groupName] = pd.Series(strainCheck, index=strains)

#Output to a file
porinCheck.to_csv('BigMatrix.tsv', sep='\t')



	

