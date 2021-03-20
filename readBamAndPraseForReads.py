import pysam as py 
import sys 
import numpy as np
import re
def alphaNumericSorting(l):
	convert = lambda text:int(text) if text.isdigit() else text
	alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)',key)]
	return sorted(l,key=alphanum_key)

def getChromosomeInfromationFromBAM(bamFile):
	chrom = {}
	bamFile=py.AlignmentFile(bamFile)
	for lines in bamFile.header['SQ']:
		chrom[lines['SN']] = int(lines['LN'])
	return chrom

def chromosomeIntegerPairing(l):
	chromosomeArray = l.keys()
	flag = 0
	intergerChromosomePair = {}
	for items in chromosomeArray:
		intergerChromosomePair[items] = flag
		flag = flag+1
	return intergerChromosomePair	

def readBamAndPutIntoNumpyArrays(bamFile,chromosome):
	bamChrom = []
	bamPositions = []
	chromosomeInfo = getChromosomeInfromationFromBAM(bamFile)
	chromosomeIntegerPair = chromosomeIntegerPairing(chromosomeInfo)

	bamRefChromosome=chromosomeIntegerPair[chromosome]
	bamFile=py.AlignmentFile(bamFile)
	for reads in bamFile:
		if reads.reference_id == bamRefChromosome:
			bamPositions.append(int(reads.pos))
	return np.asarray(bamPositions)

def readBAM(bamFile):
	return py.AlignmentFile(bamFile)
	
