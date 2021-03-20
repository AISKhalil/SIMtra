import random
import numpy as np
import pysam as py
import readBamAndPraseForReads as bamReader

def createBam(chromosome,positions,bamFile,readSize,outFile):
	chromInfo = bamReader.getChromosomeInfromationFromBAM(bamFile)
	header={'HD':{'VN':'1.0'},'SQ':[{'LN':chromInfo[chromosome],'SN':chromosome}]}
	######SAM columns; this will be used to create SAM file######
	chromosome = chromosome
	CIGAR = str(readSize)+'M'
	RNEXT = '*' 
	PNEXT = '0'
	TLEN= '0'
	DNA = ['A','T','G','C']
	qualityScores = ['>','?','@','A','B','C','D','E','F','G','H','I']
	flags = [0,16]
	MAPQs = [3,8,23,24,40,42]
	tags = (("XN",random.choice(range(0,15))),("XM",random.choice(range(0,20))),("XO",random.choice(range(0,3))),("XG",random.choice(range(0,9))),("NM",random.choice(range(0,15))),("MD",readSize),("YT","UU"))
	seq = ''.join(np.random.choice(DNA,readSize).tolist())
	qual = ''.join(np.random.choice(qualityScores,readSize).tolist())
	flag = random.choice(flags)
	mapq = random.choice(MAPQs)
	counter = 0 
	with py.AlignmentFile(outFile, "wb", header=header) as outf:
		for pos in positions: 
			a = py.AlignedSegment()
			a.query_name = "randomlyGenerated_CK_"+str(counter)
			a.query_sequence = seq
			a.flag = random.choice(flags)
			a.reference_id = 0
			a.reference_start = pos
			a.mapping_quality = mapq
			a.cigarstring = CIGAR
			a.next_reference_id = 0
			a.next_reference_start=0
			a.template_length=0
			a.query_qualities = py.qualitystring_to_array(qual)
			a.tags = (("XN",random.choice(range(0,15))),("XM",random.choice(range(0,20))),("XO",random.choice(range(0,3))),("XG",random.choice(range(0,9))),("NM",random.choice(range(0,15))),("MD",readSize),("YT","UU"))  
			outf.write(a)
			counter = counter + 1
