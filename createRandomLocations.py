import random
import numpy as np
import readBamAndPraseForReads as bamReader
import pysam as py
import pybedtools as pb

def medianReadCounter(start,end,positions,binSize=5000):
	array = range(start,end+binSize,binSize)
	positions = np.sort(positions)
	readCounts = []
	indexToStart = 0
	for i in range(len(array)-1):
		data = positions[indexToStart:][np.logical_and(positions[indexToStart:] >= array[i],positions[indexToStart:] <= array[i+1])]
		readCounts.append(len(data))
		if readCounts[-1] != 0:
			indexToStart = np.where(positions == data[-1])[0].tolist()[0]
	return np.median(readCounts)

def countReadsPerBin(start,end,positions,binSize=5000):
	array = range(start,end+binSize,binSize)
	positions = np.sort(positions)
	readCounts = []
	indexToStart = 0
	print ("Binning in progress")
	for i in range(len(array)-1):
		data = positions[indexToStart:][np.logical_and(positions[indexToStart:] >= array[i],positions[indexToStart:] <= array[i+1])]
		readCounts.append(len(data))
		if readCounts[-1] != 0:
			indexToStart = np.where(positions == data[-1])[0].tolist()[0]
	return readCounts
	
def chromosomeDivider(numberOfRegions,sizeOfChromosome,chromosome,positions,copyNumber,endRegion = 100000):
	if numberOfRegions == 1:
		start = 0
		endPosition = sizeOfChromosome
		copyNumberSelected = copyNumber
		print ("Copy Number for whole chromosome amplification will be: "+ str(copyNumberSelected))
		medianReadCounts = medianReadCounter(start,endPosition,positions)
		return [[chromosome,int(start),int(endPosition),str(endPosition-start),str(copyNumberSelected), str(medianReadCounts)]]
	else:
		copyNumber = np.asarray(copyNumber)
		np.random.shuffle(copyNumber)
		if len(copyNumber) < numberOfRegions:
			numberOfRandomRegionCheck = numberOfRegions - len(copyNumber)
			previousCN = copyNumber[-1]
			tmpCNArray = []
			for i in range(numberOfRandomRegionCheck):
				flag = 0
				while flag < 1:
					y = random.choice(copyNumber)
					if y != previousCN:
						tmpCNArray.append(y)
						previousCN = y
						flag = 1
			copyNumberInLoop = np.append(copyNumber,tmpCNArray)
		else:
			copyNumberInLoop = copyNumber
		print ("copy number combinations will be : "+str(copyNumberInLoop))
		locationsOfrandomDivisions = []
		start = 0 
		step = 0
		binSize = int(sizeOfChromosome/numberOfRegions)
		copyNumberPrevious = 0
		for i in range(numberOfRegions):
			step = step + binSize
			a  = step - endRegion
			endPosition = random.randint(a,step)
			medianReadCounts = medianReadCounter(start,endPosition,positions)
			locationsOfrandomDivisions.append([chromosome,int(start),int(endPosition),str(endPosition-start),str(copyNumberInLoop[i]),str(medianReadCounts)])
			start = endPosition
	return	locationsOfrandomDivisions


def focalCNVgenerator(listOfLargeCNVs,sizesForFocalCNV,numberOfCNVs,chromosome,sizeOfChromosome,positions,medianRatioCutOff,copyNumberThreshold):
	copyNumber =np.array([0,1,2,3,4,5,6,7,8,9,10],dtype=float)
	listOfSmallRegions = []	
	unwantedRegions = pb.BedTool("unwantedRegions.gap.blacklisted.bed")
	selfOverlapChecker = [[chromosome,0,10]]
	
	for locations in listOfLargeCNVs:
		selfOverlapChecker.append([chromosome,locations[1]-5000,locations[1]+5000])
		selfOverlapChecker.append([chromosome,locations[2]-5000,locations[2]+5000])
	#print selfOverlapChecker
	for i in range(numberOfCNVs*100):
		start = random.randint(0,sizeOfChromosome)
		end = start+random.randint(sizesForFocalCNV[0],sizesForFocalCNV[1])
		listOfSmallRegions.append([chromosome,int(start),int(end)])
		
	focalRegions = pb.BedTool(listOfSmallRegions)
	regionsClearedOfGapAndBlacklisted = focalRegions.intersect(unwantedRegions,v=True)
	selfOverlap = [[chromosome,0,10]]
	nonSelfOverlapping = []

	for items in regionsClearedOfGapAndBlacklisted:
		things = pb.BedTool([[items[0],items[1],items[2]]])
		selfCheckObj = pb.BedTool(selfOverlap)
		overlapCounter = things.intersect(selfCheckObj).count()
		if overlapCounter == 0:
			medianCounts=medianReadCounter(int(items[1]),int(items[2]),positions)
			nonSelfOverlapping.append([str(items[0]),int(items[1]),int(items[2]),str(medianCounts)])
			selfOverlap.append([items[0],items[1],items[2]])
	nonSelfOverlapping = pb.BedTool(nonSelfOverlapping)
	largeAmplificationsIntersection = nonSelfOverlapping.intersect(pb.BedTool(listOfLargeCNVs),wo=True)
	medianPassedLocations = []
	
	for locations in largeAmplificationsIntersection:
		if float(locations[9]) != 0.0:
			if (abs(float(locations[3])-float(locations[9]))/float(locations[9])) <= medianRatioCutOff:
				copyNumberRandom = random.choice(copyNumber)
				if copyNumberRandom != int(locations[8]) and abs(copyNumberRandom - int(locations[8])) >= copyNumberThreshold:
					medianPassedLocations.append([str(locations[0]),int(locations[1]),int(locations[2]),float((copyNumberRandom/float(locations[8]))-1),int(copyNumberRandom),float(locations[3]),int(locations[8])])

	medianPassedLocations = random.sample(medianPassedLocations,numberOfCNVs)
	
	return medianPassedLocations
	
