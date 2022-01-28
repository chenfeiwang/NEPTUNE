import sys
import pysam
import pyfasta

chr_mm9 = ['chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19', \
		'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chrM','chrX','chrY']
		
chr_hg19 = ['chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22',\
		'chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chrM','chrX','chrY']
		
def uniq_sequenced_reads(infile, extension = 37, species = 'mm9'):
	"""Unique the sequenced reads for downstream analysis."""
	
	unique_reads = {}
	outf = open(infile[:-4]+'.uniq.bed','w')
	
	for line in pysam.Samfile(infile,'rb'):
		if species == 'mm9':
			chr = chr_mm9[line.rname]			# chr for mouse	
		elif species == 'hg19':
			chr = chr_hg19[line.rname]			# chr for human
		start = str(line.pos)					# start
		end = line.pos
		for cg in line.cigar:
			if cg[0] == 0:
				end += cg[1]					# end by the CIGAR tag
		end = str(end)					
		if start == end:
			continue
		if line.is_reverse:						# strand
			strand = '-'
		else:
			strand = '+'

		if chr+'_'+start+'_'+end+'_'+strand in unique_reads:
			continue
		else:
			unique_reads[chr+'_'+start+'_'+end+'_'+strand] = 0
			print('\t'.join([chr, start, end, str(line.tlen), '0', strand]),file=outf)
	
	outf.close()	
	
def extend_sequenced_reads(infile, extension = 37):
	"""Retrieve the nucleosome length reads, centered then extended by 37/73 bp."""

	nucleo_reads = open(infile[:-4]+'.nucleo.%dbp.bed'%extension,'w')
	nucleo_count = 0
	for line in open(infile,'r'):
		line = line.strip().split('\t')
		if abs(int(line[3])) > 120 and abs(int(line[3])) <=180:
			center = int(int(line[1]) + abs(int(line[3]))/2)
			nucleo_count += 1
			print('\t'.join([line[0],str(max(0,center-extension)),str(center+extension),str(abs(int(line[3]))),'0','+']),file=nucleo_reads)
	nucleo_reads.close()
	print(infile+'\t'+str(nucleo_count))

def rotational_positioning(infile, fasta, sample_reads = 1000000):
	"""Calculate the AA/TT/AT frequency based on the uniq reads."""
	
	genome = pyfasta.Fasta(fasta)
	aattat = [0]*251
	total = 0
	for line in open(infile,'r'):
		line = line.strip().split('\t')
		if total >= sample_reads:
			break
		if abs(int(line[3])) > 120 and abs(int(line[3])) <= 180:
			if line[5] == '+':
				begin = int(line[1]) - 51
				end = int(line[1]) + 200
				if begin < 0 or end > len(genome[line[0]]):
					continue
				total += 1
				for i in range(begin,end):
					seq = genome[line[0]][i:i+2].upper()
					if seq == 'AT' or seq == 'AA' or seq == 'TA' or seq == 'TT':
						aattat[i-begin] += 1
			elif line[5] == '-':
				begin = int(line[2]) - 201
				end = int(line[2]) + 50
				if begin < 0 or end > len(genome[line[0]]):
					continue
				total += 1
				for i in range(begin,end):
					seq = genome[line[0]][i:i+2].upper()
					if seq == 'AT' or seq == 'AA' or seq == 'TA' or seq == 'TT':
						aattat[end-i-1] += 1
	for i in range(251):
		aattat[i] = str(aattat[i]/float(total))
	outf = open(infile[:-4]+'.RP.log','w')
	print('nuc_code<-c('+','.join(aattat)+')',file=outf)
	outf.close()

def statistical_positioning(infile, fasta, sample_reads = 100000):
	"""
	Calculate the relative nucleosome to nucleosome positions based on the uniq reads.
	"""
	genome = pyfasta.Fasta(fasta)
	tagplus, tagminus = {}, {}
	total = 0
	for line in open(infile,'r'):
		line = line.strip().split('\t')
		if total > sample_reads:
			break
		if abs(int(line[3])) > 120 and abs(int(line[3])) <= 180:
			total += 1
			if line[5] == '+':
				position = str(int(line[1]) - 1)
				if not line[0] in tagplus:
					tagplus[line[0]] = {}
				if not position in tagplus[line[0]]:
					tagplus[line[0]][position] = 1
				else:
					tagplus[line[0]][position] += 1
			elif line[5] == '-':
				position = str(int(line[2]) - 1)
				if not line[0] in tagminus:
					tagminus[line[0]] = {}
				if not position in tagminus[line[0]]:
					tagminus[line[0]][position] = 1
				else:
					tagminus[line[0]][position] += 1			
	
	same_strand_value = [0]*2001
	oppo_strand_value = [0]*2001
	for line in open(infile,'r'):
		line = line.strip().split('\t')
		if line[5] == '+':
			begin = max(int(line[1]) - 1001, 0)
			end = min(int(line[1]) + 1000, len(genome[line[0]]))
			for i in range(begin, end):
				if line[0] in tagplus and str(i) in tagplus[line[0]]:
					same_strand_value[i-begin] += tagplus[line[0]][str(i)]
				if line[0] in tagminus and str(i) in tagminus[line[0]]:
					oppo_strand_value[i-begin] += tagminus[line[0]][str(i)]
		elif line[5] == '-':
			begin = max(int(line[2]) - 1001, 0)
			end = min(int(line[2]) + 1000, len(genome[line[0]]))
			for i in range(begin, end):
				if line[0] in tagplus and str(i) in tagplus[line[0]]:
					oppo_strand_value[end-i-1] += tagplus[line[0]][str(i)]
				if line[0] in tagminus and str(i) in tagminus[line[0]]:
					same_strand_value[end-i-1] += tagminus[line[0]][str(i)]
	
	outf = open(infile[:-4]+'.SP.log','w')
	print('same_strand<-c('+','.join([str(t) for t in same_strand_value])+')',file=outf)
	print('oppo_strand<-c('+','.join([str(t) for t in oppo_strand_value])+')',file=outf)
	outf.close()
	
def main():

	if sys.argv[1] == 'uniq':
		uniq_sequenced_reads(sys.argv[2],sys.argv[3])
	elif sys.argv[1] == 'extend':
		extend_sequenced_reads(sys.argv[2],int(sys.argv[3]))
	elif sys.argv[1] == 'rotational':
		rotational_positioning(sys.argv[2],sys.argv[3])
	elif sys.argv[1] == 'statistical':
		statistical_positioning(sys.argv[2],sys.argv[3])
	
if __name__ == "__main__":
	main()
