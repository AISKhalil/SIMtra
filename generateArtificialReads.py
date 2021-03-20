import random
import numpy as np
import readBamAndPraseForReads as bamReader
import sys
def lsvArtificialReadsGenerator(listOfLocation,positions,chromosome,insertionRanges,writer):
	copyNumberHalf = []
	flag = 1
	copynumberAddedLocations = []
	
	sign = [-1,1]
	for locations in listOfLocation:
		
		writer.write(chromosome+"\t"+str(locations[1])+"\t"+str(locations[2])+"\t"+str(locations[4])+"\tlarge_CNVs_"+str(flag)+"\n")
		print ("Generating reads for LSV no "+str(flag)+" which has copy Number of "+str(locations[4]))
		#print chromosome+"\t"+str(locations[0])+"\t"+str(locations[1])+"\t"+str(locations[2])+"\t"+str(locations[3])+"\tlarge_CNVs_"+str(flag)+"\t"+str(locations[4])+"\n"
		randomlyGeneratedLocations=positions[np.logical_and(positions >= int(locations[1]),positions<=int(locations[2]))]
	
		if int(locations[4])>2:
	
			if int(locations[4]) == 3:
				
				halfSampledLocations = random.sample(list(randomlyGeneratedLocations),int(len(randomlyGeneratedLocations)/2))
				for pos in halfSampledLocations:
					checker = 0
				
					while checker < 1:
						newReadPos = pos + (random.choice(sign)*random.choice(insertionRanges))
						if newReadPos >= int(locations[1]) and newReadPos<=int(locations[2]):
							copynumberAddedLocations.append(newReadPos)
							checker = checker + 1



			elif int(locations[4]) == 4:
						
				for pos in randomlyGeneratedLocations:
					checker = 0
					#sign = 1
					while checker < 1:
						newReadPos = pos + (random.choice(sign)*random.choice(insertionRanges))
						if newReadPos >= int(locations[1]) and newReadPos<=int(locations[2]):
							copynumberAddedLocations.append(newReadPos)
							checker = checker + 1

		if int(locations[4]) == 1: 
			copyNumberHalf.append(locations)
		flag = flag+1

	copynumberAddedLocations = np.concatenate((positions,copynumberAddedLocations))

	indexToDelete = []
	if len(copyNumberHalf) > 0:
		for locations in copyNumberHalf:
			tmp = []
			condition = np.where(np.logical_and(copynumberAddedLocations>=int(locations[1]),copynumberAddedLocations<=int(locations[2])))
			tmp = []
			for x in condition:
				tmp += x.tolist()
			if tmp != None:
				indexToDelete += random.sample(tmp,int(len(tmp)/2))
	
	deletedReadsArray = np.delete(copynumberAddedLocations,indexToDelete)
	return deletedReadsArray

def focalArtificialReadsGenerator(listOfSmallRegions,positions,chromosome,insertionRanges,writer):
	#copyNumberZeros = []
	#copyNumberHalf = []
	flag = 1
	#copynumberAddedLocations = []
	focalRegionsToAdd = []
	indexToDelete = []
	for locations in listOfSmallRegions:
		writer.write(chromosome+"\t"+str(locations[1])+"\t"+str(locations[2])+"\t"+str(locations[4])+"\tsmall_CNVs_"+str(flag)+"\n")
		print ("Generating reads for focal region no "+str(flag)+" which has copy Number of "+str(locations[4])+" positioned within LSV of copy Number "+str(locations[6]))
		randomlyGeneratedLocations=positions[np.logical_and(positions >= int(locations[1]),positions<=int(locations[2]))]
		if float(locations[3]) > 0.0:
			numberOfReadsToGenerate = int((float(locations[3])*len(randomlyGeneratedLocations)))
			reSampling = [np.random.choice(randomlyGeneratedLocations) for _ in range(numberOfReadsToGenerate)]
			for pos in reSampling:
				x = 0
				sign = 1
				while x < 1:
					newReadPos = pos + (sign*random.choice(insertionRanges))
					if newReadPos >= int(locations[1]) and newReadPos<=int(locations[2]):
						focalRegionsToAdd.append(newReadPos)
						x = x + 1
					else:
						sign = -1

		if float(locations[3]) < 0.0:
			numberOfReadsToDownSample = abs(int((float(locations[3])*len(randomlyGeneratedLocations))))
			condition = np.where(np.logical_and(positions >= int(locations[1]),positions <= int(locations[2])))
			tmp = []
			for x in condition:
				tmp += x.tolist()
			if tmp != None:
				indexToDelete += random.sample(tmp,numberOfReadsToDownSample)
		flag = flag+1

	deletedReadsArray = np.delete(positions,indexToDelete)
	deletedReadsArray = np.concatenate((deletedReadsArray,focalRegionsToAdd))
	return deletedReadsArray

