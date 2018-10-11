# Goal: Convert a GFF file from Prokka into something that looks like NCBI gene table. 

import glob

def gff_to_ncbi(gff,output=None):
	file = open(gff,"r")
	if output: out = open(output,"w")
	else: 
		ncbi = gff[:-4]+"_NCBI.txt"
		out = open(ncbi,"w")
	out.write("#Replicon\tReplicon Accession\tStart\tStop\tStrand\tGeneID\tLocus\tLocus tag\tProtein product\tLength\tCOG(s)\tProtein name\n")
	ct = 0
	for line in file.readlines():
		if line.startswith("##FASTA"):
			break
		elif not line.startswith("#"):
			# GFF file doesn't specify chr vs plasmid
			out.write("-\t")
			info = line.split("\t")
			# Write the Accession, Start, Stop, and Strand. No Gene ID in the GFF files. 
			out.write("%s\t%s\t%s\t%s\t-\t" %(info[0],info[3],info[4],info[6]))
			details = info[-1].split(";")
			detdict = {}
			for d in details: 
				detdict[d.split("=")[0]]=d.split("=")[1]
			if "gene" in detdict: 
				out.write(detdict["gene"]+"\t")
			else: out.write("-\t")
			if "locus_tag" in detdict: 
				out.write(detdict["locus_tag"]+'\t-\t-\t-\t')
			else: out.write("-\t-\t-\t-\t")
			if "product" in detdict: 
				out.write(detdict["product"].rstrip()+"\n")
			else: out.write("-\n")
	file.close()
	out.close()
	
#RUN
#-------------------------
if __name__ == "__main__":
	files = glob.glob('*.gff')
	for f in files: 
		gff_to_ncbi(f)