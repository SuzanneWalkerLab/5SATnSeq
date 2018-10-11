""" 
Provides a set of functions that can be used for different kinds of TnSeq analysis. 
"""

# # Setup
# # ----------
import decimal as d 
from math import *
import numpy as np
import pandas as pd

# # Functions
# # ----------
def find_TAs(fasta):
	'''Takes in a fasta file. Returns a dataframe containing the start and stop position of each TA site. DataFrame keys: Start, Stop. 
	'''
	# Open file 
	file = open(fasta,'r')
	# Set up the variable that will hold results
	tas={'Start':[],'Stop':[]}
	# Discard the heading line 
	file.readline()
	# ct is position in the genome 
	ct = 0
	# Looking is a variable that says whether the last letter we encountered was a T. 
	looking = False
	# Loop through all the letters in the fasta, not including line breaks
	for line in file.readlines(): 
		for nt in line:
			if nt != '\n':
				ct+=1
				# If you find a T, change looking to true 
				if nt == 'T':
					looking = True
				# If the last letter was a T, and this letter is not an A, we have failed to find a TA site. 
				elif looking and nt != 'A':
					looking = False
				# If the last letter was a T, and this letter is an A, we record the site. 
				elif looking and nt == 'A':
					looking = False
					tas['Start'].append(ct-1)
					tas['Stop'].append(ct)
	return tas
	
def read_ncbi(file):
	"""
	Reads a file that lists the genes in an organism downloaded from the NCBI website. 

	REQUIRED INPUTS: 
	files: The NCBI file path. 

	OUTPUTS: 
	A dataframe with the file data. 
	"""
	# Read in the NCBI, ignoring the header row. 
	ncbi = pd.read_table(file,header=None,sep='\t',comment='#')
	# NCBI files sometimes have the COGs field and sometimes don't... 
	if ncbi.shape[1]==12: 
		ncbi.columns = ['Replicon','RefNum','Start','Stop','Strand','GeneID','Gene','Tag','Protein','Length','COG','Description']
	else: ncbi.columns = ['Replicon','RefNum','Start','Stop','Strand','GeneID','Gene','Tag','Protein','Length','Description']
	ncbi['Description']=ncbi['Description'].fillna('')
	# Get rid of the commas and tabs in the descriptions in case they will cause trouble down the line 
	ncbi['Description']=ncbi.Description.str.replace('\t','T')
	ncbi['Description']=ncbi.Description.str.replace(',','')
	ncbi['Description']=ncbi.Description.str.replace(';',':')
		
	return ncbi

def map_genes_to_TAs(ncbifile,fasta,name = False,beg = 0.03,term = 0.03):
	'''
	Reads in an NCBI file and maps the genes to the TA sites in the strain. Note that if a TA site is found in multiple overlapping genes, the site will have multiple entries. 
	
	REQUIRED INPUTS: 
	ncbifile: A file downloaded from NCBI containing the locations of genes.
	fasta: The nucleotide fasta for the organism. ONE contig only
	
	OPTIONAL INPUTS: 
	name: The name of the output file. If not provided, it will be <strain>_TAMap.txt
	beg: The proportion of the beginning of the gene that should be cut off, since insertions at the very beginning of a gene might not affect gene expression if the start of the gene is misannotated. 
	end: The proportion of the end of the gene that should be cut off, since insertions in the end of the gene might be tolerated. 
	
	OUTPUTS: 
	A dataframe is returned and also saved as a tab-delimited text file with a heading, listing the TA start, stop, the gene if it is in a gene, the strand that gene is in, a description of the gene if the TA site is in a gene, and Copy - a field which is 0 if this is the TA site's first entry and 1,2, etc. for every subsequent entry (due to gene overlap). 
	'''
	# Read in the NCBI file 
	if type(ncbifile)==str:
		ncbi = read_ncbi(ncbifile)
	else: ncbi = ncbifile
	
	# Edit the locations of the genes based on the beg and term variables. 
	ncbi['GLength'] = ncbi.Stop-ncbi.Start
	ncbi['Beg']=np.round(ncbi.Start + ncbi.GLength*beg)
	ncbi['Term']=np.round(ncbi.Stop - ncbi.GLength*term)
	ncbi.loc[ncbi.index[ncbi.Strand == '-'],'Beg'] = np.round(ncbi.Start[ncbi.Strand == '-']+ncbi.GLength[ncbi.Strand == '-']*term)
	ncbi.loc[ncbi.index[ncbi.Strand == '-'],'Term'] = np.round(ncbi.Stop[ncbi.Strand == '-']-ncbi.GLength[ncbi.Strand == '-']*beg)
	# Use the Fasta to find the locations of TA sites 
	tadict = find_TAs(fasta)
	tadf = pd.DataFrame(tadict)
	# Add gene information for the TA sites that fall within a gene 
	tadf['Gene']=''
	tadf['Strand']=''
	tadf['Description']=''
	for gene in ncbi.itertuples():
		taind = tadf.index[((tadf['Start']>=gene.Beg) & (tadf['Stop']<=gene.Term))]
		tadf.loc[taind,['Gene']]+=gene.Tag+';'
		tadf.loc[taind,['Strand']]+=gene.Strand+';'
		tadf.loc[taind,['Description']]+=gene.Description+';'
	# If a TA site is found in multiple genes (genes separated by ;), add additional entries for the TA site, each with just a single gene's information. 
	tadf['NumGenes']=tadf.Gene.str.count(';')
	tadf['Copy']=0
	numsplits = tadf['NumGenes'][tadf['NumGenes'].notnull()].unique()
	for n in numsplits:
		if n >1:
			for x in range(0,int(n)):
				splitgenes = tadf[tadf['NumGenes']==n].copy()
				splitgenes['Gene']=splitgenes['Gene'].str.split(';').str[x]+';'
				splitgenes['Strand']=splitgenes['Strand'].str.split(';').str[x]+';'
				splitgenes['Description']=splitgenes['Description'].str.split(';').str[x]+';'
				splitgenes['NumGenes']=1
				splitgenes['Copy']=x
				tadf = tadf.append(splitgenes,ignore_index=True)
	# Get rid of all the old lines that had multiple genes listed 
	tadf = tadf[tadf['NumGenes']<=1]
	tadf = tadf.drop(axis = 1, labels = 'NumGenes')
	# Get rid of semicolons. 
	tadf['Gene']=tadf.Gene.str.replace(';','')
	tadf['Strand']=tadf.Strand.str.replace(';','')
	tadf['Description']=tadf.Description.str.replace(';','')
	# Sort and save as a tab-delimited text file 
	tadf = tadf.sort_values(by = ['Start','Gene'])
	if not name: 
		name = '.'.join(fasta.split('.')[:-1])+'_TAMap.txt'
	tadf.to_csv(name,sep='\t',index=False)
	return tadf

def map_promoters_to_TAs(ncbifile,fasta,size):
	'''
	Reads in an NCBI file and maps promoters to the TA sites in the strain. Promoters are defined as the 'size' bp region upstream of a gene. Note that if a TA site is found in multiple overlapping promoter regions, the site will have multiple entries in the dataframe. 
	
	REQUIRED INPUTS: 
	ncbifile: The file name of a NCBI gene table containing the locations of genes OR a dataframe from a NCBI file that has already been read in.
	fasta: The file name for the nucleotide fasta for the organism. ONE contig only 
	
	OUTPUTS: 
	A dataframe is returned and also saved as a tab-delimited text file with a heading, listing the TA start, stop, the gene if it is in a promoter window, the strand that gene is in, a description of the gene, and Copy - a field which is 0 if this is the TA site's first entry and 1,2, etc. for every subsequent entry (due to promoter overlap). 
	'''
	# Read in the NCBI file 
	if type(ncbifile)==str:
		ncbi = read_ncbi(ncbifile)
	else: ncbi = ncbifile
	# Use the NCBI file and size argument to figure out where the promoters should be. 
	ncbi.loc[ncbi.index[ncbi.Strand == '+'],'Stop']=ncbi['Start'][ncbi.Strand == '+']
	ncbi.loc[ncbi.index[ncbi.Strand == '-'],'Start']=ncbi['Stop'][ncbi.Strand == '-']
	ncbi.loc[ncbi.index[ncbi.Strand == '+'],'Start']=ncbi['Start'][ncbi.Strand == '+']-size
	ncbi.loc[ncbi.index[ncbi.Strand == '-'],'Stop']=ncbi['Start'][ncbi.Strand == '-']+size
	# Name the output file
	fname = '.'.join(fasta.split('.')[:-1])+'_ProMap'+str(size)+'.txt'
	# Map the promoters to TA sites. map_genes_to_TAs() also saves a tab-delimited text file. 
	prodf = map_genes_to_TAs(ncbi,fasta,fname,beg=0,term=0)
	return prodf

def remap_igv(igv,newmap):
	'''
	Relabels the IGV TA sites based on a new TA map. 
	
	REQUIRED INPUTS: 
	igv: a dataframe resulting from reading in an IGV file. See read_igv() function for details on the dataframe. 
	newmap: A pandas dataframe detailing the locations of the features you are mapping to the TA sites. For more details, see the map_genes_to_TAs() function. 

	OUTPUTS: 
	Returns a re-labeled IGV dataframe.  
	'''
	# Make a copy of the TA map dataframe to begin making the new IGV dataframe. 
	newigv = newmap.copy()
	# Set read counts to 0 and get the accession number.  
	newigv['Acc']=igv.iloc[0,0]
	newigv['Total'] = 0
	newigv['Plus'] = 0
	newigv['Minus'] = 0
	# Only pay attention to one copy of each TA site from the IGV 
	if 'Copy' in igv.columns: igv = igv[igv.Copy == 0]
	# Index the igv by the TA start location
	igv = igv.set_index('Start')
	# TA sites are sometimes mapped to multiple overlapping genes, indicated in the Copy column. Loop through the copies. 
	for copynum,copygroup in newigv.groupby('Copy'):
		# Save the original index, but then set the index to the TA start site. 
		origindex = copygroup.index
		copygroup = copygroup.set_index('Start')
		# Get the read counts from the IGV data frame.  
		copygroup['Total']=igv['Total']
		copygroup['Plus']=igv['Plus']
		copygroup['Minus']=igv['Minus']
		# Reset the index to what it originally was 
		copygroup = copygroup.set_index(origindex)
		# Add the read information to the remapped IGV 
		newigv.loc[copygroup.index, 'Total'] = copygroup.Total
		newigv.loc[copygroup.index, 'Plus'] = copygroup.Plus
		newigv.loc[copygroup.index, 'Minus'] = copygroup.Minus
	return newigv
	
def make_distributions(sample,nums):
	'''
	sample: a list-like array of numbers 
	nums: a list of integers 
	Draws and sums nums values from sample 200000 times. Returns a dict with nums for keys and a list of the 200,000 sums for the values
	'''
	distdict = {}
	sample = list(sample)
	for n in nums: 
		rnsamp = np.random.choice(list(sample),(n,200000))
		prosamp = np.sum(rnsamp,axis=0)
		distdict[str(n)]=prosamp
	return distdict 
	
# def reduce_outliers(reads,percent):
	# numreads = len(reads)
	# maxrank = np.ceil(numreads*percent)
	# readrank = reads.rank(method='max',ascending=False)
	# reads[readrank <= maxrank]=reads[readrank <= maxrank].min()
	# return reads
	
def read_tabular(tabularfile):
	'''
	Uses pandas to read in a tabular file as a dataframe. DataFrame has following columns: 
		Reference 
		Position 
		Locus
		Gene
		PlusCount 
		MinusCount 
		TotalCount 
		Product
		ProteinID
		Note 
		Sequence
	This is mostly just here so I can remember what the names of the columns are. 
	'''
	tabular = pd.read_csv(tabularfile,sep='\t',header=0)
	return tabular
	
def tabular_to_igv(tabular,TAmap,direction ='new'): 
	""" 
	Makes an IGV file from a tabular file.  
	
	REQUIRED INPUTS: 
	tabular: a tabular HopCounts file downloaded from Galaxy with TnSeq data. 
	TAMap.txt file: For each TA site, has the start, stop, gene, gene description, and whether the row of data is duplicated due to the TA site being in multple genes.  
	
	OPTIONAL INPUTS: 
	direction: Can make an IGV file for just plus direction, just the minus direction, the plus and minus summed, or all three. See OUTPUTS for more details. Options: 'plus','minus','sum','all','new'. Note that the 'new' files are much larger - if you have a space constraint, use 'all' instead. 
	
	OUTPUTS:
	Saves an IGV file with the same name as the tabular file, but with the .igv extension. (e.g. 'USA300_untr_blunt.tabular' makes an igv with the name 'USA300_untr_blunt.igv')
	The IGV file has the following data columns: 
	1. Genome reference number
	2. TA start site
	3. TA end site
	4. Number of reads (for direction='all' this will be total reads)
	5. Locus tag
	6. If 'all' or 'new', plus direction reads (put here for backwards compatibility)
	7. If 'all' or 'new', minus direction reads (put here for backwards compatibility)
	8. If 'new', gene strand (+ or -). 
	9. If 'new', description of the gene. 
	10. If 'new', an indicator of whether the TA site row is a duplicate due to the TA site being in multiple genes.  
	"""
	# Read in the tabular file 
	dat = read_tabular(tabular)
	# Get the accession number for the strain 
	ref = dat.iloc[0,0]
	# Start putting together the IGV. 
	igv = TAmap.copy()
	igv['Acc']=ref
	igv['Plus']=0
	igv['Minus']=0
	# Tabular data is listed as where the read sequence ends instead of where the TA site is, so we have to make that adjustment. 
	dat['Pos1']=dat['Position']+15
	dat['Pos2']=dat['Position']+14
	dat['Pos3']=dat['Position']-15
	dat['Pos4']=dat['Position']-16
	# We gotta loop through the copies of the TA sites (since some TA sites show up multiple times by being in overlapping genes. We don't want to miss those)
	if 'Copy' not in igv.columns: igv['Copy']=0
	for n in igv.Copy.unique():
		igv2 = igv[igv['Copy']==n]
		# save the igv index
		origindex=igv2.index
		# For each of the different positions above, create a new index for the igv that is the start of each TA site and do the same with the tabular data. Find the intersection, and add the tabular data to the igv for the appropriate columns of data. Then move on to the next tabular position.  
		igv2 = igv2.set_index('Start')
		dat2 = dat.set_index('Pos1')
		datintersect = igv2.index.intersection(dat2.index)
		igv2.loc[datintersect,'Plus']+= dat2['PlusCount'][datintersect]
		dat2 = dat.set_index('Pos2')
		datintersect = igv2.index.intersection(dat2.index)
		igv2.loc[datintersect,'Plus']+= dat2['PlusCount'][datintersect]
		dat2 = dat.set_index('Pos3')
		datintersect = igv2.index.intersection(dat2.index)
		igv2.loc[datintersect,'Minus']+= dat2['MinusCount'][datintersect]
		dat2 = dat.set_index('Pos4')
		datintersect = igv2.index.intersection(dat2.index)
		igv2.loc[datintersect,'Minus']+= dat2['MinusCount'][datintersect]
		# Go back to your initial indexing and make sure you save the readcounts you just acquired. 
		igv2 = igv2.set_index(origindex)
		igv.loc[igv.index[igv['Copy']==n],'Plus']+=igv2['Plus']
		igv.loc[igv.index[igv['Copy']==n],'Minus']+=igv2['Minus']
	# The sum of the plus and minus strand reads is called Total. 
	igv['Total']=igv['Plus']+igv['Minus']
	
	# Save the IGV in different ways, depending on the input argument. So if you have old scripts that depend on the old IGV style, you can still get it output that way. 
	outname = tabular[:-7]+'igv'
	if direction=='sum':
		igv = igv[igv.Copy == 0]
		igv.to_csv('OLD_'+outname,sep='\t',columns=['Acc','Start','Stop','Total','Gene'],index=False,header=False)
	elif direction=='all':
		igv = igv[igv.Copy == 0]
		igv.to_csv(outname,sep='\t',columns=['Acc','Start','Stop','Total','Gene','Plus','Minus'],index=False,header=False)
	elif direction=='plus': 
		igv = igv[igv.Copy == 0]
		igv.to_csv('PLUS_'+outname,sep='\t',columns=['Acc','Start','Stop','Plus','Gene'],index=False,header=False)
	elif direction=='minus':
		igv = igv[igv.Copy == 0]
		igv.to_csv('MINUS_'+outname,sep='\t',columns=['Acc','Start','Stop','Minus','Gene'],index=False,header=False)
	elif direction=='new':
		igv.to_csv(outname,sep='\t',columns=['Acc','Start','Stop','Total','Gene','Plus','Minus','Strand','Description','Copy'],index=False,header=False)
	else:
		print "Invalid direction. Direction must be 'plus','minus','sum', 'all', or 'new'."

	return
	
def read_igv(igvfile):
	'''
	Reads in an IGV file. We have made several kinds of IGVs over the years. The oldest ones just had 5 columns: the accession number, the start of the TA, the end of the TA, the total number of reads, and the gene locus tag. Newer ones have the plus strand reads and minus strand reads appended as the last two columns. The Newest one is like the newer one, but it has an additional final column of 'Description' which provides the description of what the gene encodes. The script reads in the igv, detects how many columns the igv has, and returns the igv as a properly labelled dataframe. 
	Returns the IGV as a dataframe. 
	'''
	igv = pd.read_csv(igvfile,sep='\t',header=None)
	ncol = igv.shape[1]
	if ncol == 5: igv.columns = ['Acc','Start','Stop','Total','Gene']
	elif ncol == 7: igv.columns = ['Acc','Start','Stop','Total','Gene','Plus','Minus']
	else: igv.columns = ['Acc','Start','Stop','Total','Gene','Plus','Minus','Strand','Description','Copy']
	return igv
	
def combine_igvs(igvlist,outname=False):
	''' 
	Takes in a list of igv files and adds the reads columns together. Saves the resulting igv. 
	'''
	igvlist=list(igvlist)
	igv1 = read_igv(igvlist[0])
	for igv in igvlist[1:]:
		igv2 = read_igv(igv)
		igv1['Total']+=igv2['Total']
		if igv1.shape[1]>5: 
			igv1['Plus']+=igv2['Plus']
			igv1['Minus']+=igv2['Minus']
	if not outname: 
		outname = '_'.join(igvlist[0].split('_')[:-1])+'_combined.igv'
	igv1.to_csv(outname,sep='\t',header=False,index=False)
	return outname

def hitcounts(igvlist,filepath=''):
	'''
	Takes in a list of igv files and makes a csv file that lists the file, the reads in the file, the number of TA sites hit, and the percent of TA sites hit in the file. The output is simply named 'HitCounts.csv'
	'''
	countdict = {'File':[],'Reads':[],'SitesHit':[],'PercentHit':[]}
	for igvfile in igvlist: 
		countdict['File'].append(igvfile)
		igv = read_igv(igvfile)
		if igv.shape[1]>7:
			igv = igv[igv.Copy == 0]
		countdict['Reads'].append(igv['Total'].sum())
		countdict['SitesHit'].append(sum(igv['Total']>0))
		countdict['PercentHit'].append(round(np.mean(igv['Total']>0)*100,2))
	countdf = pd.DataFrame(countdict)
	countdf.to_csv(filepath+'HitCounts.csv',index=False,columns=['File','Reads','SitesHit','PercentHit'])
	return 
	
def normalize(ctrl,exp):
	"""
	Normalizes the files in case one file is much larger than the other. The experimental file is always the one normalized. Normalized by the average number of reads per TA site hit (non-zero means)
	
	REQUIRED INPUTS: 
	ctrl: A dataframe of ctrl data 
	exp: A dataframe of experimental data

	OUTPUT: 
	The exp dataframe with an additional column normreads.  
	"""
	if ctrl.shape[1]>7:
		ctrltot = ctrl['Total'][ctrl.Copy==0]
	else: ctrltot = ctrl.Total
	if exp.shape[1]>7:
		exptot = exp['Total'][exp.Copy==0]
	else: 
		exptot = exp.Total
	normfact = float(ctrltot.sum())/sum(ctrltot>0)/(exptot.sum()/sum(exptot>0))
	exp['NormedReads']=exp['Total']*normfact
	return exp

def fdr(pvals):
	"""
	Returns Benjamini-Hochberg FDR q-values corresponding to p-values.

	This function implements an algorithm equivalent to L{bh_rejected} but
	yields a list of 'adjusted p-values', allowing for rejection decisions
	based on any given threshold.
	
	INPUTS: 
	pv: list of p-values 
	
	OUTPUT: 
	list of adjusted p-values (q-values)
	"""
	# Keeping in a dataframe will mean that we won't lose the order the pvalues came in with. 
	pqdf = pd.DataFrame()
	pqdf['pval']=pvals.sort_values()
	pqdf['rank']=None
	pqdf['qval']=None
	# Get rid of null values
	realpqdf = pqdf[pqdf['pval'].notnull()]
	# Rank the pvalues from smallest to largest. 
	realpqdf['rank']=realpqdf['pval'].rank()
	# m is the number of p-values you are comparing.  
	m = realpqdf.shape[0]
	# The q-val starts off being the pval*number of pvals being compared/rank
	realpqdf['qval']=realpqdf['pval']*m/realpqdf['rank']
	qvals = list(realpqdf['qval'])
	coeff = qvals[-1]
	# Loop through from the bottom (least significant)
	for i in range(m-1, -1, -1):
		# if any q-value is greater than the q-value for the neighboring higher p-value, make the q-value equal to the neighboring higher p-value. 
		if coeff < qvals[i]:
			qvals[i] = coeff
		else: coeff = qvals[i]
	realpqdf['qval']=qvals
	# Return the q-values and nulls 
	pqdf.loc[pqdf.index[pqdf['pval'].notnull()],'qval']=realpqdf['qval']
	return pqdf['qval']
	
