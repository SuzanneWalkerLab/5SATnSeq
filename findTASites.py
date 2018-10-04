# Goal: Make a list of the TA sites in a nucleotide fasta file

import glob 

def find_tas(fasta):
	file = open(fasta,'r')
	out=open(fasta[:-6]+'_TASites.txt','w')
	out.write(";PatID\tStrand\tPattern\tSeqID\tStart\tEnd\tmatching_seq\tScore\n")
	ct = 0
	looking = False
	for line in file.readlines():
		if line.startswith('>'):
			acc = line.split(' ')[0][1:]
		else: 
			for nt in line:
				if nt != '\n':
					ct+=1
					if nt == 'T':
						looking = True
					elif looking and nt != 'A':
						looking = False
					elif looking and nt == 'A':
						looking = False
						out.write('\t'.join(['TA','D','TA',acc,str(ct-1),str(ct),'TA','1'])+'\n')
				
# # RUN
if __name__=="__main__": 
	fastas = glob.glob("*.fasta")
	for f in fastas: 
		find_tas(f)