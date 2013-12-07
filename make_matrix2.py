import glob
from Bio import SeqIO

import pandas as pd

#Use the MIC list as the master list of strains
MICdf = pd.read_csv('antibiogram_july.csv', index_col=0)
strains = MICdf.index.values
fileNames = glob.glob('*.afa')

get_strain_name = lambda seq: seq.id.split('|')[0]

porinCheck = pd.DataFrame(index=strains)

for fileName in fileNames:
	with open(fileName, 'r') as FID:
		seqs = SeqIO.to_dict(SeqIO.parse(FID, 'fasta'), key_function=get_strain_name)

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
			if (seqStart == refSeqStart) and (seqStop == refSeqStop):
				strainCheck.append('P')
			else:
				myString = ''
				if seqStart < refSeqStart:
					myString += 'E_N{}'.format(refSeqStart-seqStart)
				elif seqStart > refSeqStart:
					myString += 'T_N{}'.format(seqStart-refSeqStart)
				if seqStop > refSeqStop:
					myString += 'E_C{}'.format(seqStop-refSeqStop)
				elif seqStop < refSeqStop:
					myString += 'T_C{}'.format(refSeqStop-seqStop)
				strainCheck.append(myString)

  	porinCheck[groupName] = pd.Series(strainCheck, index=strains)

#Output to a file
porinCheck.to_csv('BigMatrix.csv', sep='\t')

# #Concatentate into one DF
# totMatrix = pd.concat(matrixDFs.values(), axis=0)

# #We now sum along the rows and discard values close to zero or close the total number
# #in order to only look at interesting ones
# counts = totMatrix.sum(axis=1)

# #Uncomment to plot the histogram of counts
# import matplotlib.pyplot as plt
# plt.hist(counts.values, 100)
# plt.show()

# cutOff = 10
# #Use line below for pandas v0.11 and above
# #interestingMatrix = totMatrix.loc[(counts>cutoff) & (counts < totMatrix.shape[1]-cutOff), :]
# interestingMatrix = totMatrix[(counts>cutOff) & (counts < totMatrix.shape[1]-cutOff)]

# #Replace the NaN's with 2
# interestingMatrix.fillna(2, inplace=True)

# #Push the result out to file
# interestingMatrix.to_csv('BigMatrix.tsv', sep='\t')

	

