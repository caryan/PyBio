import os
import subprocess

def run_flexbar(fastqFile, adapters, maxUncalled=2, minQualityScore=30, minReadLength=35, numThreads=2):
	"""
	Run flexbar to clean-up reads from a fastq file and a dictionary of adaptor sequences.
	"""
	flexBarDir = "/Users/caryan/Downloads/flexbar_v2.4_macosx"

	curDir = os.getcwd()
	os.chdir(flexBarDir)

	#Write the adapters to a fasta file
	with open("adapters.fasta", 'w') as FID:
		for k,v in adapters.items():
			FID.write('>' + k + '\n')
			FID.write(v + '\n')

	#Call flexbar 
	#Myseq supposedly puts out Illumina version 1.8 quality scores (http://seqanswers.com/forums/showthread.php?t=25450)
	output = subprocess.check_output(["./flexbar", "-r", fastqFile, "-f", "i1.8", "--pre-trim-phred", str(minQualityScore), 
						"-a", "adapters.fasta", "--min-read-length", str(minReadLength), "--max-uncalled", str(maxUncalled),
						"--threads", str(numThreads)])

	print(output)

def run_bowtie(inFile, indexFile, outFile, numThreads=2):
	#Original bowtie arguments from galaxy server
	#bowtie -q -p 4 -S -n 2 -e 70 -l 28 --maxbts 800 -k 1 -m 1 --best --phred33-quals S.aureus Sample_MHcontrol2.R1.fastq.clp.trm.qc > Sample_MHcontrol2.R1.fastq.sam
	
	# -q query input files are FASTQ .fq/.fastq (default)
	# -k 1 report up to <int> alns per read, i.e. only report the first position mapped for each read

	print('Mapping with bowtie....')
	subprocess.call(['bowtie2', "--threads", str(numThreads), "--phred33", "-k", str(1), "-x", indexFile, "-q", ,"-U", inFile, "-S", "mapped.sam"])
	print('Finshed')


adapters = {"RPI":"CAAGCAGAAGACGGCATACGAGAT",
			"RPI2":"GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA",
			"TRUESEQ1":"GATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
			"TRUESEQ2":"TATCTCGTATGCCGTCTTCTGCTTG",
			"TRUESEQ3":"ATCTCGTATGCCGTCTTCTGCTTG",
			"TRUESEQ_UNIVERSAL":"AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"}