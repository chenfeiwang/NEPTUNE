import os
import sys
import numpy

def calculate_NDR_score(infile):

	NDR = []
	for line in open(infile).xreadlines():	
		line = [float(t) for t in line.strip().split(',')]
		minusone = max(line[75:90])
		plusone = max(line[110:125])
		center = numpy.mean(line[99:102])
		if numpy.mean(line) == 0:
			NDR.append(0)
		else:
			NDR.append((max(plusone,minusone)-center)/(max(line)-min(line)))
	return NDR

def calculate_POS_score(infile):
	
	POS = []	
	for line in open(infile).xreadlines():
		line = [float(t) for t in line.strip().split(',')]
		COR = numpy.corrcoef(line[100:150],line[109:159])[0][1]
		if numpy.isnan(COR):
			POS.append(0)
		else:
			POS.append(abs(COR))
		
	return POS

def main():

	N = 33826
		
	NDR_all = {}
	POS_all = {}
	for i in range(1,int(sys.argv[2])+1):
		NDR = calculate_NDR_score('%s_siteprof%d'%(sys.argv[1],i))
		POS = calculate_POS_score('%s_siteprof%d'%(sys.argv[1],i))
		for j in range(0,N):
			if NDR_all.has_key(j):
				NDR_all[j].append(str(NDR[j]))
				POS_all[j].append(str(POS[j]))
			else:
				NDR_all[j] = [str(NDR[j])]
				POS_all[j] = [str(POS[j])]
	
	pos = 0
	outf1 = open('%s_NDR.promoter_score.txt'%sys.argv[1],'w')
	outf2 = open('%s_POS.promoter_score.txt'%sys.argv[1],'w')
	for line in open('%s_peak'%sys.argv[1],'r').xreadlines():
		line = line.strip().split('\t')
		print >>outf1, line[3]+'\t'+line[4]+'\t'+'\t'.join(NDR_all[pos])
		print >>outf2, line[3]+'\t'+line[4]+'\t'+'\t'.join(POS_all[pos])
		pos += 1
	outf1.close()
	outf2.close()

if __name__ == "__main__":
	main()
