# Goal: Convert tabular files to wig files that can be used with the TRANSIT essential gene identification platform. 
# To run: python make_wigs.py TAList.txt data.tabular

import sys

def make_TAdict(TAfile):
	TAs = {}
	for line in open(TAfile,'r').readlines():
		if not line.startswith(';'):
			info = line.split('\t')
			pos = int(info[4])
			TAs[pos]=0
	return TAs

def tab_to_wig(tabular,TAdict):
	out = open(tabular[:-7]+'wig', 'w')
	out.write("variableStep  chrom=chr1\n")
	file = open(tabular,'r')
	for line in file: 
		info = line.split("\t")
		pos = 0 
	# # Make sure you're looking at a line with real data
		if len(info)>= 6 and info[1].isdigit():
		# # Look for reads in the positive direction and change site to TA site instead of sequencing start site. 
			if int(info[4]) != 0:
				pos = int(info[1])
				if pos+15 in TAdict:
					TAdict[pos+15] += int(info[4])
				elif pos+14 in TAdict:
					TAdict[pos+14] += int(info[4])
		##Do the same with the reads that are in the reverse direction and change site to TA site instead of sequencing start site.  
			if int(info[5]) != 0:
				pos = int(info[1])
				if pos-15 in TAdict: 
					TAdict[pos-15] += int(info[5])
				elif pos-16 in TAdict:
					TAdict[pos-16] += int(info[5])
	file.close()
	
	keys = TAdict.keys()
	keys.sort()

	## Everything gets written to the output IGV file.  
	for k in keys:
		out.write('%d\t%d\n' %(k,TAdict[k]))
	
	out.close()
	return
	
#RUN
#-------------------------
if __name__ == "__main__":
	TAs = make_TAdict(sys.argv[1])
	tab_to_wig(sys.argv[2],TAs)