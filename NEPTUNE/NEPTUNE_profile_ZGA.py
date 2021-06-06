import sys
import itertools
import numpy

MalePeak = "MalePSperm_peak"
FemalePeak = "FemalePMII_peak"

MaleProfile = ["MalePSperm_siteprof1","MalePSperm_siteprof2","MalePSperm_siteprof3","MaleP0.5_siteprof1","MaleP0.5_siteprof2", \
               "MaleP1_siteprof1","MaleP1_siteprof2","MaleP1_siteprof3","MaleP1.5_siteprof1","MaleP1.5_siteprof2", "MaleP1.5_siteprof3", \
               "MaleP2_siteprof1","MaleP2_siteprof2","MaleP2_siteprof3","MaleP3_siteprof1","MaleP3_siteprof2",\
               "MaleP4_siteprof1","MaleP4_siteprof2","MaleP4_siteprof3","MaleP6_siteprof1","MaleP6_siteprof2","MaleP6_siteprof3","MaleP6_siteprof4","MaleP6_siteprof5",\
               "MaleP8_siteprof1","MaleP8_siteprof2","MaleP8_siteprof3","MaleP8_siteprof4","MaleP12_siteprof1","MaleP12_siteprof2","MaleP12_siteprof3"]
FemaleProfile = ["FemalePMII_siteprof1","FemalePMII_siteprof2","FemalePMII_siteprof3","FemalePMII_siteprof4","FemalePMII_siteprof5","FemaleP0.5_siteprof1","FemaleP0.5_siteprof2","FemaleP0.5_siteprof3",\
                 "FemaleP1_siteprof1","FemaleP1_siteprof2","FemaleP1_siteprof3","FemaleP1.5_siteprof1","FemaleP1.5_siteprof2","FemaleP1.5_siteprof3",\
                 "FemaleP2_siteprof1","FemaleP2_siteprof2","FemaleP2_siteprof3","FemaleP3_siteprof1","FemaleP3_siteprof2","FemaleP3_siteprof3","FemaleP3_siteprof4",\
                 "FemaleP4_siteprof1","FemaleP4_siteprof2","FemaleP4_siteprof3","FemaleP4_siteprof4","FemaleP6_siteprof1","FemaleP6_siteprof2","FemaleP6_siteprof3","FemaleP6_siteprof4","FemaleP6_siteprof5",\
                 "FemaleP8_siteprof1","FemaleP8_siteprof2","FemaleP8_siteprof3","FemaleP8_siteprof4","FemaleP12_siteprof1","FemaleP12_siteprof2","FemaleP12_siteprof3"]

def readProfile(peak, profile, subset = 'all'):
	
	id = []
	for line in open(peak, 'r').xreadlines():
		line = line.strip().split('\t')
		id.append(line[3]+'_'+line[4])
		
	id_value = {}
	line_count = 0
	for line in open(profile, 'r').xreadlines():
		line_value = [float(t) for t in line.strip().split(',')]
		id_value[id[line_count]] = line_value
		line_count += 1
		
	flatID = list(itertools.chain(*id_value.values()))
	upperID = sorted(flatID)[int(len(flatID)*0.99)+1]
	outline = []
	if subset == 'all':
		for k in id_value.keys():
			for i in range(0,len(id_value[k])):
				if id_value[k][i] > upperID:
					id_value[k][i] = upperID
		for pos in zip(*id_value.values()):
			outline.append(numpy.mean(pos))
	else:
		outmatrix = []
		for line in open(subset,'r').xreadlines():
			line = line.strip().split('\t')
			if id_value.has_key(line[3]+'_'+line[4]):
				for i in range(0,len(id_value[line[3]+'_'+line[4]])):
					if id_value[line[3]+'_'+line[4]][i] > upperID:
						id_value[line[3]+'_'+line[4]][i] = upperID
				outmatrix.append(id_value[line[3]+'_'+line[4]])
		for pos in zip(*outmatrix):
			outline.append(numpy.mean(pos))
	
	return ','.join([str(t) for t in outline])

def generateScripts(gene):

	outf = open('profile_%s.r'%gene,'w')
	for male in MaleProfile:
		print >>outf, "%s<-c("%male+readProfile(MalePeak, male, gene)+")"
	for female in FemaleProfile:
		print >>outf, "%s<-c("%female+readProfile(FemalePeak, female, gene)+")"
	outf.close()

def main():

	generateScripts("all")
	generateScripts("mm9.minor_ZGA.bed")
	generateScripts("mm9.major_ZGA.bed")
	
if __name__ == "__main__":
	main()

