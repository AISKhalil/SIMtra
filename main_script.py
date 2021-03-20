import random
import pybedtools as pb
import matplotlib.pyplot as plt
import numpy as np
import readBamAndPraseForReads as bamReader
import argparse
import sys
import pysam as py
import createRandomLocations as rl
import generateArtificialReads as ge
import bamFileCreator as bw
import re

parser =argparse.ArgumentParser()
parser.add_argument('--bam',dest='input_bam_file',required=True, help='BAM/SAM file to artificially spike chromosome(s) with CNV(s)')
parser.add_argument('--medianRatio',dest='median_ratio_cut_off',required=True, help='Median Ratio cut-off to be used between the focal CNV and LSV')
parser.add_argument('--cnThreshold',dest='copy_number_cut_off',required=False, help='Minimun difference between copy number of LSV and focal',default="1")
#parser.add_argument('--noOfFocalDeletions',dest='focal_del_no',required=True, help='Number of deletions to add')
#parser.add_argument('--noOfFocalAmplifications',dest='focal_amp_no',required=True, help='Number of amplifications to add')
parser.add_argument('--outSAM',dest='out_sam_file',required=True, help='Name of output samfile to create')
parser.add_argument('--outSimulatedCNVBed',dest='out_bed_file',required=True, help='Name of the bed12 format that will have informations of the simulated CNVs')
parser.add_argument('--largecnvNumber',dest='large_num_cnvs',required=True, help='Number of Large Copynumber Variation to introduce into the chromosome, if multiple chromosome are entered, then give the number of LSV for each corresponding chromosome in a comma seperated list. eg, 4,4,2,1 where 1 being whole chromosome')
parser.add_argument('--focalcnvNumber',dest='small_num_cnvs',required=True, help='Number of focal CNVs to intersperse within the large cnv')
parser.add_argument('--chr',dest='chromosome_to_spike',required=True, help='Chromosome to introduce CNV into; Please follow the same chromosome description as provided in the BAM file, if multiple chromsome are provided then give in comma separated list. eg: chr1,chr3,chr8,chr11')
parser.add_argument('--readSize',dest='size_of_read',required=True, help='size of the reads to be generated; ideally should be the same as your real data')
parser.add_argument('--largecnvSize',dest='sizes_of_large_csv',required=False, help='Size range for large CNV(s); please provide in bases; for eg: 5000000-10000000(default)',default="5000000-10000000")
parser.add_argument('--largefocalSize',dest='sizes_of_focal_csv',required=False, help='Size range for focal CNV(s); please provide in bases; for eg: 50000-100000(default)',default="50000-100000")
parser.add_argument('--largecopyNumbers',dest='range_of_copyNumber',required=True, help='range of copy number to simulate, if multiple chromosome are entered, then give the number of copy for each corresponding lsv in each corresponding chromosome in a comma seperated list. eg: 1-4,1-4,3,4')
parser.add_argument('--insertionRange',dest='range_of_distance',required=False, help='range of distance to be given to the randomly inserted reads from the original read position;Please provide the minimum,maximum, and step seperated by a hyphen ;eg : 100-2000-50(default)',default="10-1000-10")
if len(sys.argv)==1:
	parser.print_help(sys.stderr)
	sys.exit(1)


######Getting Information from user
information = parser.parse_args()
bamFileName = information.input_bam_file
medianRatioCutOff = float(information.median_ratio_cut_off)
cnCutOff = int(information.copy_number_cut_off)
outputSAMFile = information.out_sam_file
outBEDFile = information.out_bed_file
largeCNVno = list(map(int,information.large_num_cnvs.split(',')))
focalCNVno = int(information.small_num_cnvs)
sizesForLargeCNV = list(map(int,information.sizes_of_large_csv.split('-')))
sizesForFocalCNV = list(map(int,information.sizes_of_focal_csv.split('-')))

copyNumberList = []
tmpArray = list(map(str,information.range_of_copyNumber.split(',')))
for items in tmpArray:
	if re.search("-",items) == None:
		copyNumberList.append(int(items))
	else:
		copyNumberList.append(range(int(items.split('-')[0]),int(items.split('-')[1])+1))
insertionRanges = range(list(map(int,information.range_of_distance.split('-')))[0],list(map(int,information.range_of_distance.split('-')))[1]+list(map(int,information.range_of_distance.split('-')))[2],list(map(int,information.range_of_distance.split('-')))[2])
chromosome = information.chromosome_to_spike.split(',')
readSize = int(information.size_of_read)

bedFile = open(outBEDFile,"w+")
zipped = zip(chromosome,largeCNVno,copyNumberList)

for chrm,LSV_no,LSV_copy_no in zipped:
       # print ("Working for chromosome: "+str(chrm)+" with ")

	bamObj = py.AlignmentFile(bamFileName,"r")
	chromosomeInformation = bamReader.getChromosomeInfromationFromBAM(bamFileName)
	chromosomeIntegerNumber = bamReader.chromosomeIntegerPairing(bamReader.getChromosomeInfromationFromBAM(bamFileName))[chrm]

	print ("Parsing BAM file for reads for "+str(chrm))
	bamPositions = []
	for reads in bamObj:
		if reads.reference_id == chromosomeIntegerNumber:
			bamPositions.append(reads.pos)
	if len(bamPositions) == 0:
		print (str(chrm)+" doesn't exist")
		continue
	bamPositions=np.asarray(bamPositions)
	print ("Total number of reads in bam File: "+str(len(bamPositions)))
	print ("Done parsing BAM file for reads for "+str(chrm))


	print ("Generating large CNVs for chromosome "+str(chrm)+" ....")
	listOfLocation = rl.chromosomeDivider(LSV_no,chromosomeInformation[chrm],chrm,bamPositions,LSV_copy_no)
	print ("Generating focal CNVs for chromosome  ....")
	listOfSmallRegions = rl.focalCNVgenerator(listOfLocation,sizesForFocalCNV,focalCNVno,chrm,chromosomeInformation[chrm],bamPositions,medianRatioCutOff,cnCutOff)
	print ("Generating artificial reads for LSV")
	deletedReadsArray = ge.lsvArtificialReadsGenerator(listOfLocation,bamPositions,chrm,insertionRanges,bedFile)
	print ("Generating artificial reads for focal")
	deletedReadsArray = ge.focalArtificialReadsGenerator(listOfSmallRegions,deletedReadsArray,chrm,insertionRanges,bedFile)
	print ("Total number of reads in after read removal File: "+str(len(deletedReadsArray)))

	print ("Writing BAM files")
	outputBAMName = outputSAMFile+chrm+"_Artificial_reads.bam"
	bw.createBam(chrm,deletedReadsArray,bamFileName,readSize,outputBAMName)
	print (str(chrm)+" is all done\n\n")
bedFile.close()
