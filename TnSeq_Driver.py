"""
This performs Mann-Whitney U analysis on your data to find genes that are enriched in reads or depleted of reads. Then uses a bootstrapping approach to identify upregulated genes. 

WHAT YOU NEED: 
1. TABULAR OR IGV FILES 
Tabular files are tab-delimited files with the following columns (with header): 
	Reference - the accession number for the contig
	Position - the start location of the sequence
	Locus - the locus tag if the position is within a gene
	Gene - the gene name if one exists
	PlusCount - the plus strand sequencing read count
	MinusCount - the minus strand sequencing read count
	TotalCount - PlusCount+MinusCount
	Product - Description of the protein product
	ProteinID - the ID number for the protein 
	Note - a more detailed description of the of the protein product
	Sequence - The sequence - this field is usually empty. 
IGV files are tab-delimited files with no heading. The Walker lab has three versions of these floating around. They have the following columns: 
	Acc - The accession number for the genome (all versions)
	Start - The start site of the TA site (all versions)
	Stop - The stop site of the TA site (all versions)
	Total - The total number of reads at that site (all versions)
	Gene - The locus tag if the TA site is within a gene (all versions)
	Plus - The reads in the plus strand (versions 2+ only)
	Minus - The reads in the minus strand (versions 2+ only)
	Description - A description of the gene product (version 3 only)
	Strand - The strand the gene is in (version 3 only)
The names of the tabular and igv files must be styled as follows: organism-GroupID_treatment_tn
organism - This is simply the name of the organism or strain you are working with (HG003,USA300). This is mostly important for tabular files because it is how you match your tabular files to the NCBI files.  
GroupID - This is an identifier that helps you group things. For example, if you have two different controls from the same strain, one at 37 and one at 30 degrees, you can label your files USA300-30 and USA300-37. All files with the same organism-GroupID combination will be compared to the matching control file. 
treatment - for your control files, this HAS to be 'untr'. For your treated files, it can be just about anything (just no underscores).
tn - the transposon construct
NOTE: IF NO ORGANISM IS PROVIDED, IT IS ASSUMED TO BE HG003. HOWEVER, IF YOU HAVE A GROUP ID, YOU MUST INCLUDE AN ORGANISM.
ALSO NOTE: If the tabular files/IGV files are in a directory separate from where your script is located, you can provide a path using the argument -p. 
EG: >> python TnSeq_Driver.py -p C:\Katie\Documents\TnSeq\ 
	
2. IF USING TABULARS: A TAMap.txt FILE FOR EACH ORGANISM/STRAIN OR A FASTA AND NCBI FILE FOR EACH ORGANISM/STRAIN
The TAMap.txt file should be named <Organism>_TAMap.txt (e.g. HG003_TAMap.txt). Columns are: 
	Start - Start of the TA site 
	Stop - End of the TA site 
	Gene - Locus tag if the TA site falls in a gene
	Strand - The strand of the gene if the TA site is in a gene
	Description - Description of the protein product of the gene
	Copy - If a TA site is found in multiple overlapping genes, it is given one row for each gene it is found in. The copy indicates whether the row you are looking at is the first (0) or duplicate (1,2,etc.) row containing that TA site. 
If you do not have a TAMap file for your organism/strain, the script will make one if you have a nucleotide fasta (named <Organism>.fasta) and a gene table downloaded from NCBI (named <Organism>_NCBI.txt). 
NOTE: The TAMap/NCBI/Fasta files can be in either the directory where the script is located OR the folder where your IGVs/Tabulars are (specified with -p when you run the command) 

2. IF DOING BOOTSTRAP ANALYSIS: A PROMAP FOR EACH ORGANISM/STRAIN OR A FASTA AND NCBI FILE FOR EACH ORGANISM/STRAIN
The ProMap file should be named <Organism>_ProMap<size>.txt (e.g. HG003_ProMap500.txt). Columns are: 
	Start - Start of the TA site 
	Stop - End of the TA site 
	Gene - Locus tag if the TA site falls in the promoter region a gene (size bp upstream of the gene) 
	Strand - The strand of the gene if the TA site is in the promoter region of a gene
	Description - Description of the protein product of the gene
	Copy - If a TA site is found in multiple overlapping promoters, it is given one row for each gene it is found in. The copy indicates whether the row you are looking at is the first (0) or duplicate (1,2,etc.) row containing that TA site. 
If you do not have a ProMap file for your organism/strain, the script will make one if you have a nucleotide fasta (named <Organism>.fasta) and a gene table downloaded from NCBI (named <Organism>_NCBI.txt). Default size is 500 bp.  
NOTE: The TAMap/NCBI/Fasta files can be in either the directory where the script is located OR the folder where your IGVs/Tabulars are (specified with -p when you run the command)

COMMAND LINE ARGUMENTS:
-h	Displays a help message for running the script.
-i	Use this if you just want to make IGVs that you can use for some other purpose. 
-it	IGV type you want it to output. There are three types of IGVs floating around. Some scripts might not be compatible with newer IGV formats. 1: Has accession number, TA start, TA stop, locus tag, and total reads. 2. Has all type 1 fields plus plus strand reads and minus strand reads. 3. Has the same columns as 2, plus the strand of the gene, a description of the gene product, and a copy number for the TA site, as this type of IGV has multiple lines for TA sites that are found in multiple genes. 
-m	Use this if you only want to do the MWU analysis and don't want to do the upregulation analysis. 
-u	Use this if you only want to do the upregulation analysis and don't want to do the MWU analysis. 
-s <winsize>	Size of the promoter windows for the upregulation analysis. 
-p <path>	Use this if your script is stored in a separate file from your data files. <path> should be the path to your data folder. 

OUTPUTS : 
1. HitCounts text file that lists the the number of reads, and the percent of TA sites hit for each condition. 
2. IGV files for each tabular file, a combined igv file that sums the reads from all the transposons for each condition, and promoters igv file that sums the reads from the cap, pen, and tuf transposons for each condition. IGV is tab-delimited, no header. IGV fields are: 
	- NCBI genome reference number
	- Start of TA site
	- End of TA site
	- # Total Reads
	- Locus Tag (if in gene)
	- # Plus Direction Reads (only in newer igvs)
	- # Minus Direction Reads (only in newer igvs) 
	- Strand of the gene (only in newest igvs)
	- Description  (only in newest igvs) that describes the protein product
	- Copy  (Only in newest IGVs. see TAMap description above for details on this field)
	
3. A MWU csv file (starts with MWU_) for each experiment/control pair of blunt/combined files that has the following information: 
	- Gene: Locus Tag
	- Start: Start of gene
	- Stop: End of gene 
	- TAs: # of TA sites
	- CtrlReads: Control Reads
	- ExpReads: Experimental Reads
	- NormedExpReads: Normalized experimental reads - normalized by average number reads per TA site hit. 
	- Ratio: NormedExpReads+1 / CtrlReads+1 
	- PVal: MWU p-value
	- CorrPVal: FDR q-value
	- Description: protein product of the gene - This will only be present if you use new (type 3) IGVs or start from tabular files. 
4. An UPREG csv file (starts with UPREG) for each experiment/control pair of cap/pen/tuf/promoters files that has the following information: 
	- Gene: Locus tag 
	- Start: First TA site in promoter region 
	- Stop: Last TA site in promoter region 
	- Strand: The strand the gene is found on 
	- Hit: Whether the gene classifies as a hit. To be a plus strand hit, the plus strand q-value has to be less than 0.05, there has to be more than 1 TA site in the plus strand that has reads in both the control and treated file, and the read difference for the minus strand has to be below the 90th percentile. The converse is true for minus strand reads. 
	- Upreg: Whether the hit has a matching promoter direction and gene direction. 'YES' or ''
	- PlusTAs: Number of TA sites in the promoter region that have reads in the plus strand in both the control file and the treated file. 
	- CtrlPlus: Number of control plus strand reads in the qualifying plus strand TA sites. 
	- ExpPlus: Number of treated plus strand reads in the qualifying plus strand TA sites. 
	- PlusDiff: ExpPlus-CtrlPlus 
	- PlusPVal: Bootstrap-derived p-value for PlusDiff 
	- PlusQVal: Benjamini-Hochberg corrected p-value for plus strand
	- MinusTAs: Number of TA sites in the promoter region that have reads in the minus strand in both the control file and the treated file. 
	- CtrlMinus: Number of control minus strand reads in the qualifying minus strand TA sites. 
	- ExpMinus: Number of treated minus strand reads in the qualifying minus strand TA sites. 
	- MinusDiff: ExpMinus-CtrlMinus
	- MinusPVal: Bootstrap-derived p-value for MinusDiff 
	- MinusQVal: Benjamini-Hochberg corrected p-value for minus strand

EX: 
>> python TnSeq_Driver.py

"""

# # Setup
# # ----------
import argparse
import glob
import os.path
import numpy as np
import pandas as pd
import scipy.stats
import sys
import time
import TnSeqTools as tn
import warnings

# # Help
# # ----------
# # If you need help with the script, type the command  >> python TnSeq_Driver.py -h    to see the following help message: 
parser = argparse.ArgumentParser(description="Takes in tabular files or IGV files and conducts Mann-Whitney U (MWU) and false discovery rate (FDR) analysis to identify enriched/depleted reads within genes. Also performs bootstrapping analysis to identify upregulated genes, also followed by FDR. Requires that either the folder contains a set of IGV files or a set of tabular files, NCBI files, and TA site files or that you specify a path to a directory that contains these things. For descriptions of these files, see the associated guide. EX: >> python TnSeq_Driver.py")
parser.add_argument('-p','--path',help='The location of the IGV/Tabular files that you are analyzing if they are not in your current directory.',nargs=1,default=[''],type=str)
parser.add_argument('-i','--makeigvs',help='Type this if you only want the script to make IGVs for you and you do not care about the others functions of the script.',action="store_true",default=False)
parser.add_argument('-it','--igvtype',help='The type of IGV you want to make from the tabulars. Takes in a 1, 2, or 3. 1 is the Marina style of IGV, which has the accession number, the start of the TA site, the end of the TA site, the total number of reads, and the locus tag if the TA site is in a gene. If a TA site is in multiple genes, it is considered to be in the first gene. 2 is just like 1 except that it has the added columns of plus strand reads and minus strand reads. 3 has all the columns of 2, but if a TA site is found in multiple genes it appears as multiple lines in the IGV, one per gene. Also includes a description column and a Copy column that tells you whether this is the first entry for a TA site found in multiple genes (0) or a subsequent copy of the TA site (1,2, etc.). Default is 3.',nargs=1,default=[3],type=int)
parser.add_argument('-m','--MannWhitneyU',help='Type this if you only want to do Mann-Whitney U analysis and you do not care about the upregulation analysis.',action="store_true",default=False)
parser.add_argument('-u','--Upregulation',help='Type this if you only want to do Upregulation analysis and you do not care about the Mann-Whitney U analysis.',action="store_true",default=False)
parser.add_argument('-s','--PromoterSize',help='Use this to indicate how long (bp) you want your promoter windows to be. Default = 500.',nargs=1,default=[500],type=int)

# Parse arguments: 
args = parser.parse_args()

# Get path: 
filepath = args.path[0]
if (len(filepath)>0) & (not filepath.endswith('\\')):
	filepath+='\\'

# Decide what type of IGV you are outputting 
igvtype = args.igvtype[0]
if igvtype == 1: igvtype = 'sum'
elif igvtype == 2: igvtype = 'all'
else: igvtype = 'new'

# Figure out what kinds of analysis we are doing: 
if args.makeigvs and not args.MannWhitneyU and not args.Upregulation: 
	mwu = False 
	upreg = False 
elif args.MannWhitneyU and not args.Upregulation:
	mwu = True 
	upreg = False 
elif args.Upregulation and not args.MannWhitneyU:
	mwu = False 
	upreg = True
else: 
	mwu = True 
	upreg = True

# Determine size of promoter windows for upregulation analysis
prosize = args.PromoterSize[0]

# # Check filepath 
if filepath != '':
	if not os.path.exists(filepath): 
		print 'ERROR: The directory that you provided does not exist. Please check for typos.'
		sys.exit()

# # Get Ready!
# # ----------

# # Start the timer
now = time.strftime('%I:%M %p',time.localtime())
print "Started at %s." %(now)
starttime = time.clock()

# # Make sure the data is in the directory. 
igvs = glob.glob(filepath+'*.igv')
tabulars = glob.glob(filepath+'*.tabular')
if tabulars==[] and igvs==[]: 
	print "ERROR: No tabular or igv files found in directory."
	sys.exit()


# # Collect File Name Info
# # ----------

# # Find all the tabular/igv files in the folder. Determine what organisms are present (HG003, USA300, etc.), what groups are present, what treatments have been used, and what promoters are present. Expecting files to be named "Organism-GroupID_Treatment_Promoter.tabular". The GroupID is optional. An example is HG003-30degC_untr_blunt.tabular or USA300_oxacillin_pen.tabular

filedict = {'File':[],'Organism':[],'GroupID':[],'Treatment':[],'Tn':[],'Type':[]}
for file in tabulars+igvs: 
	filedict['File'].append(file)
	filedict['Type'].append(file.split('.')[-1])
	info = file.split('\\')[-1].split('/')[-1].split('_')
	if len(info)==2:
		filedict['Organism'].append('HG003')
		filedict['GroupID'].append('')
	elif len(info) <2 or len(info)>3: 
		print "ERROR: Tabular files are inappropriately named. Need to be either organism-GroupID_treatment_promoter.tabular or, if strain is HG003 and no GroupID, can be named treatment_promoter.tabular."
		sys.exit()
	else: 
		groupinfo = info[0].split('-')
		filedict['Organism'].append(groupinfo[0])
		if len(groupinfo)==1:
			filedict['GroupID'].append('')
		else: filedict['GroupID'].append('-'.join(groupinfo[1:]))
	filedict['Treatment'].append(info[-2])
	filedict['Tn'].append(info[-1].split('.')[0])

filedf = pd.DataFrame(filedict)

# Make sure you don't run into problems later finding matches for files: 
filedf['Tn']=filedf.Tn.str.lower()
filedf['Tn'][filedf.Tn == 'bcdept']='combined'
filedf['Tn'][filedf.Tn == 'cpt']='promoters'

# # MAKE IGVS
# # ----------
# For all tabulars, check if there is a corresponding IGV. If so, ignore the tabular. 
for tabular in filedf[filedf['Type']=='tabular'].itertuples(): 
	igvname = tabular.File[:-7]+'igv'
	if os.path.exists(igvname):
		filedf = filedf.drop(labels=tabular.Index)


# If an IGV doesn't exist for a tabular, make one:
# Loop by Organism
for org, orggroup in filedf.groupby('Organism'):
	# Check if there are tabulars that need to be converted for the organism
	if sum(orggroup['Type']=='tabular')>0:
		# If there are tabulars...
		# Look for a TAMap file in the script directory and data directory
		if os.path.exists(org+'_TAMap.txt'): 
			tas = pd.read_csv(org+'_TAMap.txt',sep='\t')
		elif os.path.exists(filepath+org+'_TAMap.txt'):
			tas = pd.read_csv(filepath+org+'_TAMap.txt',sep='\t')
		# If there is no TAMap for the organism in either location, make one using provided fasta and NCBI file. 
		elif (os.path.exists(org+'.fasta')) & (os.path.exists(org+'_NCBI.txt')): 
			tamapstart = time.clock()
			print 'Making TA map for %s. This will probably take a couple minutes...' %(org)
			tas = tn.map_genes_to_TAs(org+'_NCBI.txt',org+'.fasta')
			tamapstop = time.clock()
			timedif = tamapstop-tamapstart
			print 'This step took %d minutes and %d seconds.' %(round(timedif/60),timedif%60)
		elif (os.path.exists(filepath+org+'.fasta')) & (os.path.exists(filepath+org+'_NCBI.txt')): 
			tamapstart = time.clock()
			print 'Making TA map for %s. This will probably take a couple minutes...' %(org)
			tas = tn.map_genes_to_TAs(filepath+org+'_NCBI.txt',filepath+org+'.fasta')
			tamapstop = time.clock()
			timedif = tamapstop-tamapstart
			print 'This step took %d minutes and %d seconds.' %(round(timedif/60),timedif%60)
		# If you can't find the appropriate files, print an error statement. 
		else:
			print "ERROR: For each organism, you must have a TAMap.txt file or an NCBI and fasta file with which to make the TAMap. The file(s) must be found in the directory the script is located in."
			sys.exit()
		# Once the TAMap is found/created, convert the tabulars to IGVs.
		print 'Making IGV files from the tabular files for %s.' %(org)
		for tabfile in orggroup[orggroup.Type == 'tabular'].itertuples():
			tn.tabular_to_igv(tabfile.File,tas,igvtype)
			filedf.loc[tabfile.Index,'Type']='igv'
			filedf.loc[tabfile.Index,'File']=filedf.loc[tabfile.Index,'File'][:-7]+'igv'

# Combine appropriate IGVs 
for groupname, igvgroup in filedf.groupby(['Organism','GroupID','Treatment']):
	# If we are going to be doing MWU analysis, add matching blunt, cap, erm, dual, pen, and tuf files to get a combined file.
	if (sum(igvgroup['Tn']=='combined')==0) and mwu:
		combinefiles = igvgroup['File'][igvgroup['Tn'].isin(['blunt','cap','dual','erm','pen','tuf'])]
		if len(combinefiles)<6:
			promoters = list(igvgroup['Tn'][igvgroup['Tn'].isin(['blunt','cap','dual','erm','pen','tuf'])])
			promstring = ', '.join(promoters)
			print 'NOTE: For %s %s %s, not all data were found for the combined file. Only able to combine %s' %(groupname[0],groupname[1],groupname[2],promstring)
		if len(combinefiles)>1:
			fname = tn.combine_igvs(combinefiles)
			filedf=filedf.append({'File':fname,'Organism':groupname[0],'GroupID':groupname[1],'Treatment':groupname[2],'Tn':'combined','Type':'igv'},ignore_index=True)
	# If we are going to be doing upreg analysis, add matching cap, pen, and tuf files to get a promoters-only file.
	if (sum(igvgroup['Tn']=='promoters')==0) and upreg:
		combinefiles = igvgroup['File'][igvgroup['Tn'].isin(['cap','pen','tuf'])]
		if len(combinefiles)<3:
			promoters = list(igvgroup['Tn'][igvgroup['Tn'].isin(['cap','pen','tuf'])])
			promstring = ', '.join(promoters)
			print 'NOTE: For %s %s %s, not all data were found for the promoter combined file. Only able to combine %s' %(groupname[0],groupname[1],groupname[2],promstring)
		if len(combinefiles)>1:
			proigvname = filepath + groupname[0]
			if groupname[1] != '':
				proigvname += '-'+groupname[1]
			proigvname+='_'+groupname[2]+'_promoters.igv'
			fname = tn.combine_igvs(combinefiles,proigvname)
			filedf=filedf.append({'File':fname,'Organism':groupname[0],'GroupID':groupname[1],'Treatment':groupname[2],'Tn':'promoters','Type':'igv'},ignore_index=True)
			
print "IGVs have been created." 

# # PERFORM ANALYSES
# # ----------

# Count hit sites and reads in all of the files. 
print 'Making HitCounts file.' 
tn.hitcounts(filedf.File,filepath)

# Perform Mann-Whitney U analysis, if requested. 
if mwu: 
	print 'Performing Mann-Whitney U analysis and false discovery rate correction.'
	mwustart = time.clock()
	# Perform MWU and FDR on each set of matching files (same organism, same group id, same transposon). Only done for blunt and combined. 
	for groupinfo,groupfiles in filedf[filedf.Tn.isin(['blunt','combined'])].groupby(['Organism','GroupID','Tn']):
		# Read in the control file for the group 
		try: ctrl = tn.read_igv(groupfiles['File'][groupfiles['Treatment'].str.count('untr')==1].iloc[0])
		except: 
			print 'No control found for %s %s %s. Those files were skipped for MWU analysis.' %(groupinfo[0], groupinfo[1],groupinfo[2])
			continue
		# Loop through the non-control files in the group
		for igvfile in groupfiles.itertuples():
			if 'untr' not in igvfile.Treatment:
				# Read in the treated file, normalize it, and add the untr data to the treated data frame. 
				trt = tn.read_igv(igvfile.File)
				trt = tn.normalize(ctrl,trt)
				trt['CtrlReads']=ctrl['Total']
				# Make a TAs column so you can count the number of TAs in a gene 
				trt['TAs']=1
				# Get rid of TA sites not found in a gene 
				trt=trt[(trt.Gene != '') & pd.notnull(trt.Gene)]
				# Group TA sites by the gene they are in 
				genegroups=trt.groupby('Gene')
				# Aggregate the dataframe so that each gene has one line. Note that we have to do this a different way depending on whether this is an old or new format IGV. 
				if 'Description' in trt.columns: 
					mwuout = genegroups.agg({'Start':min,'Stop':max,'Total':sum,'CtrlReads':sum,'NormedReads':sum,'TAs':sum,'Description':max})
				else: mwuout = genegroups.agg({'Start':min,'Stop':max,'Total':sum,'CtrlReads':sum,'NormedReads':sum,'TAs':sum})
				# Do Mann-Whitney U analysis. 
				pvals = []
				for genename,genegroup in genegroups: 
					try: u,pv = scipy.stats.mannwhitneyu(genegroup['CtrlReads'],genegroup['NormedReads'],alternative='two-sided')
					except: pv = 1
					pvals.append(pv)
				mwuout['PVal']=pvals
				# Do FDR analysis on the MWU p-values 
				mwuout['CorrPVal']=tn.fdr(mwuout['PVal'])
				# Calculate the treated:untreated read ratio 
				mwuout['Ratio']=np.round((mwuout['NormedReads']+1)/(mwuout['CtrlReads']+1),4)
				# Save the dataframe as a csv 
				mwuout = mwuout.rename(index=str,columns={'NormedReads':'NormedExpReads','Total':'ExpReads'})
				outname = filepath+'MWU_'+igvfile.Organism 
				if igvfile.GroupID != '':
					outname+='-'+igvfile.GroupID
				outname+='_'+igvfile.Treatment+'_'+igvfile.Tn+'.csv'
				# Again, if you are using older IGV files, it won't have a Description column, but newer IGVs will. 
				if 'Description' in trt.columns:
					mwuout.to_csv(outname,columns=['Start','Stop','TAs','CtrlReads','ExpReads','NormedExpReads','Ratio','PVal','CorrPVal','Description'])
				else: mwuout.to_csv(outname,columns=['Start','Stop','TAs','CtrlReads','ExpReads','NormedExpReads','Ratio','PVal','CorrPVal'])
	mwutime = time.clock()-mwustart
	print 'Mann-Whitney U analysis completed. This step took %d minutes and %d seconds.' %(round(mwutime/60),mwutime%60)

# Perform upregulation analysis, if requested. 	
if upreg: 
	print 'Performing upregulation bootstrap analysis and false discovery rate correction.'
	upregstart = time.clock()
	# Perform bootstrap analysis on each set of matching files (same organism, same group id, same transposon). Done for all files except blunt and combined. 
	for groupinfo,groupfiles in filedf[~filedf.Tn.isin(['blunt','combined','erm'])].groupby(['Organism','GroupID','Tn']):
		org = groupinfo[0]
		# Look for a ProMap file in the script directory and data directory
		if os.path.exists('%s_ProMap%d.txt' %(org,prosize)): 
			pros = pd.read_csv('%s_ProMap%d.txt' %(org,prosize),sep='\t')
		elif os.path.exists('%s%s_ProMap%d.txt' %(filepath,org,prosize)):
			pros = pd.read_csv('%s%s_ProMap%d.txt' %(filepath,org,prosize),sep='\t')
		# If there is no ProMap for the organism in either location, make one using provided fasta and NCBI file. 
		elif (os.path.exists(org+'.fasta')) & (os.path.exists(org+'_NCBI.txt')): 
			promapstart = time.clock()
			print 'Making promoter map for %s. This will probably take a couple minutes...' %(org)
			pros = tn.map_promoters_to_TAs(org+'_NCBI.txt',org+'.fasta',prosize)
			promapstop = time.clock()
			timedif = promapstop-promapstart
			print 'This step took %d minutes and %d seconds.' %(round(timedif/60),timedif%60)
		elif (os.path.exists(filepath+org+'.fasta')) & (os.path.exists(filepath+org+'_NCBI.txt')): 
			promapstart = time.clock()
			print 'Making promoter map for %s. This will probably take a couple minutes...' %(org)
			pros = tn.map_promoters_to_TAs(filepath+org+'_NCBI.txt',filepath+org+'.fasta',prosize)
			promapstop = time.clock()
			timedif = promapstop-promapstart
			print 'This step took %d minutes and %d seconds.' %(round(timedif/60),timedif%60)
		# If you can't find the appropriate files, print an error statement. 
		else:
			print "ERROR: For each organism, you must have a promoter map txt file or an NCBI and fasta file with which to make the promoter map. The file(s) must be found in the directory the script is located in or the directory the data is located in."
			sys.exit()
		# Read in the control file for the group 
		try: ctrl = tn.read_igv(groupfiles['File'][groupfiles['Treatment'].str.count('untr')!=0].iloc[0])
		except: 
			print 'No control found for %s %s %s. Those files were skipped for upregulation analysis.' %(groupinfo[0], groupinfo[1],groupinfo[2])
			continue
		ctrl = tn.remap_igv(ctrl,pros)
		# Set cutoff for outliers: 
		newmax = ctrl.Total.sum()/10000
		
		# Loop through the non-control files in the group
		for igvfile in groupfiles.itertuples():
			if 'untr' not in igvfile.Treatment:
				# Read in the treated file, remap it to the promoters rather than genes, and add the untr data to the treated data frame. 
				trt = tn.read_igv(igvfile.File)
				trt = tn.remap_igv(trt,pros)
				trt['CtrlPlus']=ctrl['Plus']
				trt['CtrlMinus']=ctrl['Minus']
				
				# Only include TA sites hit in either file 
				# Make a TAs column so you can count the number of TAs in a promoter 
				trt['PlusTAs']=0
				trt.loc[trt.index[(trt.Plus>0)|(trt.CtrlPlus>0)],'PlusTAs']=1
				trt['MinusTAs']=0
				trt.loc[trt.index[(trt.Minus>0)|(trt.CtrlMinus>0)],'MinusTAs']=1
				
				# Reduce outliers:
				for col in ['Plus','CtrlPlus','Minus','CtrlMinus']:
					trt.loc[trt.index[trt[col]>newmax],col]=newmax
				
												
				# Get rid of TA sites not found in a promoter 
				trt=trt[(trt.Gene != '') & pd.notnull(trt.Gene)]
				
				# Get rid of reads for TA sites not in both files.
				trt.loc[trt.index[trt.PlusTAs==0],['Plus','CtrlPlus']]=0
				trt.loc[trt.index[trt.MinusTAs==0],['Minus','CtrlMinus']]=0
				
				# Find read differences
				trt['PlusDiff']=trt['Plus']-trt['CtrlPlus']
				trt['MinusDiff']=trt['Minus']-trt['CtrlMinus']
				
				# Make gene info prettier: (get rid of ;)
				trt['Gene']=trt.Gene.str.replace(';','')
				trt['Description']=trt.Description.str.replace(';','')
				trt['Strand']=trt.Strand.str.replace(';','')
				
				# Aggregate the dataframe so that each gene has one line. 
				bootout = trt.groupby('Gene').agg({'Start':min,'Stop':max,'Plus':sum,'Minus':sum,'CtrlPlus':sum,'CtrlMinus':sum,'PlusTAs':sum,'MinusTAs':sum,'PlusDiff':sum,'MinusDiff':sum,'Description':max,'Strand':max})
				
				# Do Bootstrap analysis. 
				bootout['PlusPVal'] = 1
				bootout['MinusPVal'] = 1
				plusdists = tn.make_distributions(trt['PlusDiff'][trt.PlusTAs==1],bootout['PlusTAs'].unique())
				minusdists = tn.make_distributions(trt['MinusDiff'][trt.MinusTAs==1],bootout['MinusTAs'].unique())
				
				for row in bootout.itertuples():
					if row.PlusTAs>0:
						curplus = plusdists[str(row.PlusTAs)]
						bootout.loc[row.Index,'PlusPVal']=np.mean(curplus>row.PlusDiff)
					if row.MinusTAs>0:
						curminus = minusdists[str(row.MinusTAs)]
						bootout.loc[row.Index,'MinusPVal']=np.mean(curminus>row.MinusDiff)
				
				# Do FDR analysis on the bootstrap p-values 
				bootout['PlusQVal']=tn.fdr(bootout.PlusPVal)
				bootout['MinusQVal']=tn.fdr(bootout.MinusPVal)
				bootout.loc[bootout.index[bootout.PlusQVal>1],'PlusQVal']=1
				bootout.loc[bootout.index[bootout.MinusQVal>1],'MinusQVal']=1
				
				# Define hits as those with a Qval in one strand less than 0.05 and a read difference in the other strand below the 90th percentile. 
				# Must also have more than one hit TA 
				bootout['Hit']=''
				if groupinfo[-1] != 'dual':
					plus90 = bootout.PlusDiff.quantile(0.90)
					minus90 = bootout.MinusDiff.quantile(0.90)
					bootout.loc[bootout.index[(bootout.PlusQVal<=0.05)&((bootout.MinusDiff==0)|(bootout.MinusDiff<=minus90))&(bootout.PlusTAs>1)],'Hit'] = 'PLUS'
					bootout.loc[bootout.index[(bootout.MinusQVal<=0.05)&((bootout.PlusDiff==0)|(bootout.PlusDiff<=plus90))&(bootout.MinusTAs>1)],'Hit'] = 'MINUS'
					bootout['Upreg']=''
					bootout.loc[bootout.index[(bootout.Hit == 'PLUS')&(bootout.Strand == '+')],'Upreg']='YES'
					bootout.loc[bootout.index[(bootout.Hit == 'MINUS')&(bootout.Strand == '-')],'Upreg']='YES'
				else: 
					
					bootout.loc[bootout.index[(bootout.PlusQVal<=0.05)&(bootout.PlusTAs>1)],'Hit'] = 'HIT'
					bootout.loc[bootout.index[(bootout.MinusQVal<=0.05)&(bootout.MinusTAs>1)],'Hit'] = 'HIT'
					bootout.loc[bootout.index[bootout.Hit == 'HIT'],'Upreg'] = 'YES'
				
				# Save the dataframe as a csv 
				bootout = bootout.rename(index=str,columns={'Plus':'ExpPlus','Minus':'ExpMinus'})
				outname = filepath+'UPREG%d_' %(prosize)+igvfile.Organism 
				if igvfile.GroupID != '':
					outname+='-'+igvfile.GroupID
				outname+='_'+igvfile.Treatment+'_'+igvfile.Tn+'.csv'
				bootout.to_csv(outname,columns=['Start','Stop','Strand','Hit','Upreg','PlusTAs','CtrlPlus','ExpPlus','PlusDiff','PlusPVal','PlusQVal','MinusTAs','CtrlMinus','ExpMinus','MinusDiff','MinusPVal','MinusQVal','Description'])
	upregtime = time.clock()-upregstart
	print 'Done with upregulation analysis. This step took %d minutes and %d seconds.' %(round(upregtime/60.),upregtime%60)
	
# # Final progress update  					
end = time.strftime('%I:%M %p',time.localtime())
endtime = time.clock()				
print "Completed at %s. Run time was approximately %d minutes and %d seconds." %(end,(endtime-starttime)/60, (endtime-starttime)%60)