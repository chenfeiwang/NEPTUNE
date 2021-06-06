import sys

def binned_region(input):
	"""merge the nucleosome tags of 200bp bins into fixed length bins. Need input the bin length, from 1000bp+."""
	step = [1000,40000]
	for bin_length in step:
		nucleo = {}
		outf = open(input[:-4]+'.%sbp.txt'%bin_length,'w')
		for line in open(input, 'r').xreadlines():
			line = line.strip().split('\t')
			if line[0] == 'chr':
				print >>outf, '\t'.join(line)
			else:
				pos = int(line[1])/bin_length
				if nucleo.has_key(line[0]+'_'+str(pos*bin_length)):
					nucleo[line[0]+'_'+str(pos*bin_length)].append(float(line[3]))
				else:
					nucleo[line[0]+'_'+str(pos*bin_length)] = [float(line[3])]
	
		for k,v in nucleo.iteritems():
			print >>outf, k.replace('_','\t')+'\t'+str(sum(v))
		outf.close()	

def main():

	binned_region(sys.argv[1])

if __name__ == "__main__":
	main()
