# Packages: 
import pandas as pd 
import glob 

# Functions: 
# Read in a TRANSIT-style prot_table file. 
def read_prot_table(filename):
	prottable = pd.read_csv(filename,sep='\t',header=None,names=['Description','Start','Stop','Strand','ProtLength','col1','col2','Gene','Tag','col3'])
	return prottable
	
# Read in a TRANSIT-style wig file. 
def read_wig(filename):
	wig = pd.read_csv(filename,sep='\t',header=0,names=['TA','Reads'])
	return wig 

# Executable:
# Open the list of corresponding tags
tags = pd.read_csv('MasterTagList.csv')

# Collect the names of the prot_table files 
PTs = glob.glob('*.prot_table')

# Loop through the prot_table files and change the gene name to the Roary tag
for prottable in PTs: 
	strain = prottable.split('.')[0]
	pt = read_prot_table(prottable)
	for row in pt.itertuples():
		if row.Tag != '-':
			roarytag = tags.Name[tags[strain+'ProkkaTag']==row.Tag].values[0]
			pt.loc[row.Index,'Gene']=roarytag
	pt.to_csv('Roary_'+prottable,sep='\t',header=False,index=False)

# Gather the names of the wig files
wigs = glob.glob('*.wig')
strains = list(set([x.split('_')[0] for x in wigs]))

# For each wig file, label the TA sites with the gene they are found in, using the Roary tags
# We do this strain-wise because the labeling takes a while, but it should be the same for every file for the strain, so we can save time by just reusing the results for the first file in a strain. 
for strain in strains:
	genes = ''
	for w in wigs: 
		if strain in w:
			wig = read_wig(w)
			if type(genes)==str:
				wig['Gene']=''
				strain = w.split('_')[0]
				prot = read_prot_table('Roary_'+strain+'.prot_table')
				for gene in prot.itertuples():
					wig.loc[wig.index[(wig.TA>gene.Start) & (wig.TA<gene.Stop)],'Gene']=gene.Gene
				
				genes = wig.Gene
			else: wig['Gene']=genes
			wig.to_csv('Roary_'+w[:-3]+'csv',index=False)




