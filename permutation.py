import pandas as pd
import glob
import itertools
import random
import numpy as np
import time

# FUNCTIONS
# ----------

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

	
# ANALYSIS
# ----------
# Gather all of the annotated wig files, saved as csvs. Add together the ones for the same strain.  
wigs = glob.glob('Roary_*.csv')
strains = list(set([x.split('_')[1] for x in wigs]))
wigs2 = []
for s in strains:
	dat = ''
	for w in wigs: 
		if s in w: 
			if type(dat)==str: dat = pd.read_csv(w)
			else: 
				dat2 = pd.read_csv(w)
				dat['Reads']+=dat2['Reads']
	dat.to_csv(s+'_CombinedData.csv')
	wigs2.append(s+'_CombinedData.csv')

print wigs2
# Get all the combinations of files
filecombos = itertools.combinations(wigs2,2)
print filecombos

# Loop through the combinations 
for f in filecombos:
	print 'Performing permutation test on %s and %s' %(f[0],f[1])
	start = time.clock()
	# Read in the data 
	dat1 = pd.read_csv(f[0])
	dat2 = pd.read_csv(f[1])
	
	# Determine which one has more reads per TA hit
	non0mean1 = float(dat1.Reads.sum())/(sum(dat1.Reads>0))
	non0mean2 = float(dat2.Reads.sum())/(sum(dat2.Reads>0))
	
	# Normalize to the one with more hit TAs
	if sum(dat1.Reads>0) >= (sum(dat2.Reads>0)): 
		dat2['Reads']=dat2.Reads*non0mean1/non0mean2 
	else: dat1['Reads']=dat1.Reads*non0mean2/non0mean1
	
	# Filter to only include genes 
	dat1 = dat1[(dat1.Gene != '') & (dat1.Gene.notnull())]
	dat2 = dat2[(dat2.Gene != '') & (dat2.Gene.notnull())]
	
	# Find common genes  
	genes1 = set(dat1.Gene)
	genes2 = set(dat2.Gene)
	commongenes = list(genes1.intersection(genes2))
	
	# For each common gene, perform a permutation test to get a two-sided p-value.
	# Set up the columns for the output DF
	pvallist = []
	diflist = []
	numTAs1 = []
	numTAs2 = []
	sum1 = []
	sum2 = []
	# Loop through the genes 
	for g in commongenes: 
		# Get the reads from all of the TA sites in the gene the two files to then find the read sums for the gene 
		testdat1 = dat1.Reads[dat1.Gene==g]
		testdat2 = dat2.Reads[dat2.Gene==g]
		sum1.append(sum(testdat1))
		sum2.append(sum(testdat2))
		tot = sum(testdat1)+sum(testdat2)
		dif = abs(sum(testdat1)-sum(testdat2))
		# Store the difference in the gene's read sums for the two files 
		diflist.append(dif)
		# Pool the reads in the gene for both files. 
		pooleddat = list(testdat1)+list(testdat2)
		# Determine how many TA sites there are in the gene in the first file 
		num1 = len(testdat1)
		numTAs1.append(num1)
		numTAs2.append(len(testdat2))
		difdist = []
		# Perform the permutation to generate a distribution of gene read differences. 
		for i in range(20000):
			cursamp = random.sample(pooleddat,num1)
			samp1 = sum(cursamp)
			samp2 = tot-samp1
			difdist.append(abs(samp1-samp2))
		difdist = np.array(difdist)
		# Calculate the p-value as the proportion of distribution read differences greater (abs value) than the actual read difference
		pval = np.mean(difdist>=dif)
		pvallist.append(pval)
		
	# Arrange data in a DF
	outdf = pd.DataFrame({'Gene':commongenes,'TAs1':numTAs1,'TAs2':numTAs2,'Reads1':sum1,'Reads2':sum2,'Difference':diflist,'pVal':pvallist})
	# Correct for multiple hypothesis testing ala Benjamini-Hochberg
	outdf['qVal']=fdr(outdf.pVal)
	
	# Write results to file
	outdf.to_csv(f[0].split('_')[0]+'_vs_'+f[1].split('_')[0]+'.csv',index=False, columns = ['Gene','TAs1','TAs2','Reads1','Reads2','Difference','pVal','qVal'])
	
	# Status update
	timetaken = time.clock()-start 
	print 'Test took %d minutes and %d seconds.' %(round(timetaken/60),timetaken%60)
	
		
		