import pandas as pd
from collections import defaultdict
import glob


SNPsStrains = defaultdict(list)
strains = []

# fileNames = glob.glob("*.xls")
# for fileName in fileNames:
#
# 	df = pd.read_excel(fileName, 0)
#
# 	strain = fileName.split()[0]
# 	strains.append(strain)
# 	for i,s in df.iterrows():
# 		#We'll store the SNP as a tuple of (RefPos, Type, Length, )
# 		SNP = (s["Reference Position"], s["Type"], s["Length"], s["Reference"], s["Allele"])
# 		SNPsStrains[SNP].append(strain)

excelFile = pd.read_excel("AZPAE12416.xlsx", sheetname=None)

for sheetName, sheet in excelFile.iteritems():
	strain = sheetName[:10]
	strains.append(strain)
	for i,s in sheet.iterrows():
		#We'll store the SNP as a tuple of (RefPos, Type, Length, )
		SNP = (s["Reference Position"], s["Type"], s["Length"], s["Reference"], s["Allele"])
		SNPsStrains[SNP].append(strain)

#Create new data frame
SNPs = SNPsStrains.keys()

#Fill the SNP columns
dfBis = pd.DataFrame()
dfBis["Ref. Pos."] = pd.Series([x[0] for x in SNPs], index=SNPs)
dfBis["Type"] = pd.Series([x[1] for x in SNPs], index=SNPs)
dfBis["Length"] = pd.Series([x[2] for x in SNPs], index=SNPs)
dfBis["Reference"] = pd.Series([x[3] for x in SNPs], index=SNPs)
dfBis["Allele"] = pd.Series([x[4] for x in SNPs], index=SNPs)

#Fill in the strain columns
for strain in strains:
	dfBis[strain] = pd.Series(['P' if strain in v else 'A' for v in SNPsStrains.values()], index=SNPs)

#Add a column for the total
dfBis["Num. Present"] = pd.Series([len(x) for x in SNPsStrains.values()], index=SNPs)

#Sort by position for ease of viewing
dfBis.sort_index(by="Ref. Pos.", inplace=True)

dfBis.to_csv("SNPCollection_PA.tsv", sep='\t', index=False)

#Pull out those that have multiple SNPs at the same location
SNPPosCounts = dfBis["Ref. Pos."].value_counts()

dfMulti = dfBis[dfBis["Ref. Pos."].map(lambda x: SNPPosCounts[x]) > 1]
dfMulti.to_csv("SNPMultiples_PA.tsv", sep='\t', index=False)

"""Repeat for table of changes"""

dfBis = pd.DataFrame()
dfBis["Ref. Pos."] = pd.Series([x[0] for x in SNPs], index=SNPs)
dfBis["Type"] = pd.Series([x[1] for x in SNPs], index=SNPs)
dfBis["Length"] = pd.Series([x[2] for x in SNPs], index=SNPs)
dfBis["Reference"] = pd.Series([x[3] for x in SNPs], index=SNPs)
dfBis["Allele"] = pd.Series([x[4] for x in SNPs], index=SNPs)

#Fill in the strain columns
for strain in strains:
	dfBis[strain] = pd.Series([dfBis["Allele"][k] if strain in v else float('nan') for k,v in SNPsStrains.items()], index=SNPs)

#Add a column for the total
dfBis["Num. Present"] = pd.Series([len(x) for x in SNPsStrains.values()], index=SNPs)

#Sort by position for ease of viewing
dfBis.sort_index(by="Ref. Pos.", inplace=True)

dfBis.to_csv("SNPCollection_Changes.tsv", sep='\t', index=False)

#Pull out those that have multiple SNPs at the same location
SNPPosCounts = dfBis["Ref. Pos."].value_counts()

dfMulti = dfBis[dfBis["Ref. Pos."].map(lambda x: SNPPosCounts[x]) > 1]
dfMulti.to_csv("SNPMultiples_Changes.tsv", sep='\t', index=False)
