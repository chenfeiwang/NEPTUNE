import sys
import itertools
import numpy

MalePeak = "MalePN_peak"
FemalePeak = "FemalePN_peak"

MaleProfile = ["MalePN_siteprof1","MalePN_siteprof2","MalePN_siteprof3","MalePN_siteprof4","MalePN_siteprof5",\
				"MalePN_siteprof6","MalePN_siteprof7","MalePN_siteprof8","MalePN_siteprof9","MalePN_siteprof10",]
FemaleProfile = ["FemalePN_siteprof1","FemalePN_siteprof2","FemalePN_siteprof3","FemalePN_siteprof4","FemalePN_siteprof5",\
				"FemalePN_siteprof6","FemalePN_siteprof7","FemalePN_siteprof8","FemalePN_siteprof9","FemalePN_siteprof10",]

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
	generateScripts("data/mm9.minor_ZGA.bed")
	generateScripts("data/mm9.major_ZGA.bed")
	
if __name__ == "__main__":
	main()

