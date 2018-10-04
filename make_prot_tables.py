# # Goal: Take gff output from Prokka and make .prot_table required input for TRANSIT Gumbel. 

# # Packages
import glob

# # Function
def make_prot_table(gff):
	input = open(gff,'r')
	output = open(gff[:-3]+"prot_table",'w')
	for line in input.readlines(): 
		if line.startswith("##FASTA"):
			break
		elif not line.startswith("#"):
			info = line.split("\t")
			if info[2] != 'repeat_region':
				start = info[3]
				stop = info[4]
				strand = info[6]
				details = info[-1].split(';')
				detdict = {}
				for d in details: 
					detdict[d.split("=")[0]]=d.split("=")[1]
				try:
					tag = detdict["ID"]
				except: print line
				if "Name" in detdict: 
					gene = detdict["Name"]
				else: gene = '-'
				prod = detdict["product"].rstrip()
				prodlen = int((int(stop)-int(start))/3)
				for thing in [prod,start,stop,strand,prodlen,0,0,gene,tag]:
					output.write(str(thing)+'\t')
				output.write('-\n')
# # RUN
if __name__=="__main__": 
	gffs = glob.glob("*.gff")
	for g in gffs: 
		make_prot_table(g)