source /media/ahmed/Data2/CNAtra_2019/simulatedDataEstPython/bin/activate
cd /media/ahmed/Data2/CNAtra_2019/simulatedData/simDataGenerator
# IMR32 generation SAM with artifical CNVs
python main_script.py --bam /media/ahmed/Data2/CNAtra_2019/realData/oldData/IMR32.bam --outSAM /media/ahmed/Data2/CNAtra_2019/simulatedData/Data/artificialChrOnly/16March/IMR32_ --outSimulatedCNVBed /media/ahmed/Data2/CNAtra_2019/simulatedData/Data/artificialChrOnly/16March/IMR32.artificialChrs.bed --medianRatio 0.05 --chr chr4,chr5,chr8,chr9,chr14,chr18,chr19,chr20,chr21,chr22 --readSize 101 --largecnvNumber 1,5,5,4,4,3,3,2,2,1 --largecopyNumbers 3,1-4,1-4,1-4,1-4,2-4,2-4,3-4,3-4,4 --focalcnvNumber 40


# CHP212 generation SAM with artifical CNVs
python main_script.py --bam /media/ahmed/Data2/CNAtra_2019/realData/oldData/CHP212.bam --outSAM /media/ahmed/Data2/CNAtra_2019/simulatedData/Data/artificialChrOnly/16March/CHP212_ --outSimulatedCNVBed /media/ahmed/Data2/CNAtra_2019/simulatedData/Data/artificialChrOnly/16March/CHP212.artificialChrs.bed --medianRatio 0.05 --chr chr4,chr8,chr12,chr13,chr18,chr19,chr21,chr22 --readSize 101 --largecnvNumber 1,5,5,4,4,3,2,1 --largecopyNumbers 3,1-4,1-4,1-4,1-4,2-4,3-4,4 --focalcnvNumber 40


python main_script.py --bam /media/ahmed/Data2/CNAtra_2019/realData/oldData/HepG2.bam --outSAM /media/ahmed/Data2/CNAtra_2019/simulatedData/Data/artificialChrOnly/16March/HepG2_ --outSimulatedCNVBed /media/ahmed/Data2/CNAtra_2019/simulatedData/Data/artificialChrOnly/16March/HepG2.artificialChrs.bed --medianRatio 0.05 --chr chr11,chr12,chr15,chr19,chr21,chr22 --readSize 76 --largecnvNumber 1,5,4,3,2,1 --largecopyNumbers 3,1-4,1-4,2-4,3-4,4 --focalcnvNumber 40


