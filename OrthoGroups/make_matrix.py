import glob
from Bio import SeqIO
import pandas as pd
from progressbar import Percentage, ProgressBar, RotatingMarker, ETA

MICdf = pd.read_csv('antibiogram_july.csv', index_col=0)
fileNames = glob.glob('*.afa')

matrixDFs = {}
#This can take a long time so put up a progress bar
widgets = ['Working: ', Percentage(), ' ', Bar(marker=RotatingMarker()),
               ' ', ETA(), ' ']
pbar = ProgressBar(widgets=widgets, maxval=len(fileNames))
pbar.start()
for ct, fileName in enumerate(fileNames):
	pbar.update(ct)
	with open(fileName, 'r') as FID:
		seqs = list(SeqIO.parse(FID, 'fasta'))

	#Find the reference record first
	refName = 'PAO1'
## 	seqNames = []
## 	for seq in seqs:
## 		seqNames.append(seq.id.split('|')[0])
## 		if seqNames[-1] == refName:
## 			refSeq = seq.seq.tostring()
	for seq in seqs:
		if seq.id.split('|')[0] == refName:
			refSeq = seq.seq.tostring()
			break

	#Now loop through again and output the matrix
## 	with open(fileName.split('.')[0] + '.matrix', 'w') as FID:
## 		for seq in MICdf.index:
## 			FID.write(seq + ', ')
## 			FID.write(str(MICdf['meropenem'][seq])+', ')
## 			if seq in seqNames:
## 				FID.write(', '.join(['0' if (c1 == c2) else '1' for c1,c2 in zip(seqs[seqNames.index(seq)].seq.tostring(), refSeq)]))
## 			else:
## 				FID.write(', '.join(['2']*len(refSeq)))
## 			FID.write('\n')

	#Create a DF comparing each to the reference strain
	groupName = fileName.split('.')[0]
  	matrixDFs[groupName] = pd.DataFrame(
	{
	seq.id.split('|')[0]:
	pd.Series([0 if (c1 == c2) else 1
		   for c1,c2 in zip(seq.seq.tostring(), refSeq)],
		  index = ['V' + str(ct+1) + '_' + groupName for ct in range(len(refSeq))]
		  )
	for seq in seqs })

#Concatentate into one DF
totMatrix = pd.concat(matrixDFs.values(), axis=0)

#We now sum along the rows and discard values close to zero or close the total number
#in order to only look at interesting ones
counts = totMatrix.sum(axis=1)

#Uncomment to plot the histogram of counts
import matplotlib.pyplot as plt
plt.hist(counts.values, 100)
plt.show()

cutOff = 10
interestingMatrix = totMatrix.loc[(counts>cutOff) & (counts < totMatrix.shape[1]-cutOff), :]

#Replace the NaN's with 2
interestingMatrix.fillna(2, inplace=True)

#Push the result out to file
interestingMatrix.to_csv('BigMatrix2.tsv', sep='\t')

	

