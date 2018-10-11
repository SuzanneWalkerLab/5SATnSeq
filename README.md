# 5SATnSeq
Scripts and files needed to reproduce the analyses described in !!!paper citation!!!. 

Note: All python scripts were written for python 2.7 and require some or all of the following libraries: glob, itertools, numpy, pandas, random, time, sys

## Sequence Data Processing

FASTQ files obtained from sequencing facility are processed using Galaxy and in-house python scripts to ultimately obtain tabular files that map Tn sequencing reads to individual sequence sites. 


### Files Needed
- __*.FASTQ__ - The raw sequencing reads. Our FASTQs can be downloaded from the Sequence Reads Archive, accession ID !!!. 
- __Galaxy-Workflow-TnSeq.ga__ - The workflow used in Galaxy to process the FASTQ files into SAM files. 
- __Galaxy-SampleBarcodes.txt__ - An input for the Galaxy workflow that contains the six sample barcodes used to split the raw FASTQ from the pooled sample into the FASTQs for individual samples. 
- __Galaxy-TransposonBarcodes.txt__ - An input for the Galaxy workflow that contains the six transposon barcodes used to split each individual sample FASTQ by transposon construct. 
- __Filter1.txt__ - An input for the Galaxy workflow used to get rid of unmatched reads. 
- __Filter2.txt__ - An input for the Galaxy workflow used to get rid of unmatched reads.
- __Filter3.txt__ - An input for the Galaxy workflow used to get rid of unmatched reads.
- __*.fasta__ - A chromosome nucleotide fasta for the genome you are mapping to. See paper methods for the NCBI accession numbers for the strains used.  
- __sam_to_tabular.py__ - A python script that converts all SAM files in the directory to tab-delimited hop count files. 

### Process
1. Sign into usegalaxy.org and upload FASTQ, workflow, barcodes, filters, and fasta(s). 
2. Run workflow, following the prompts to select the appropriate inputs. If the FASTQ file has samples from multiple genomes, run the workflow once with each fasta file and only save the appropriate SAM files.
3. Download the resulting SAM files into a directory containing the sam_to_tabular.py script. 
4. Run the sam_to_tabular.py script `python sam_to_tabular.py`
5. Rename the tabular files so that they are compatible with downstream scripts: strain, underscore, treatment (untr if control), underscore, transposon construct (blunt, cap, dual, erm, pen, tuf). If you have multiple samples from the same strain but different control conditions, you can group the files to a particular control by adding an identifier after the strain, demarcated by a hyphen. Examples: USA300-RPMI_untr_blunt.tabular HG003_daptomycin-15_cap.tabular

## Essential Gene Analysis

We identified the essential genes using the TRANSIT software Gumbel method (1). The scripts below are those used to convert files so that they could be used with the TRANSIT software and to perform a permutation test to determine which genes have significant fitness differences across strains of *S. aureus*. 

### Files Needed
- __*.tabular__ - The tabular files acquired from processing the FASTQ files (see above).
- __*.fasta__ - A chromosome nucleotide fasta for the genome you are mapping to. See paper methods for the NCBI accession numbers for the strains used.
- __*.gff__ - A GFF file output from running Prokka (2) on the fasta files.
- __make_prot_tables.py__ - A python script that converts the GFF files from Prokka into prot_table files for TRANSIT. 
- __findTASites.py__ - For each fasta file in the directory, makes a list of all of the TA sites available.
- __make_wig.py__ - Makes a wig file recognized by the TRANSIT software from a tabular file and a list of TA sites created using findTASites.py.
- __labelWIGs.py__ - Using the master tag list, changes the gene names in the prot_table files into Roary tags that can be compared across all strains. Then it uses these new prot_tables to annotate the genes in the wig files. For this script to work, the prot_tables have to be named just with the name of the strain and the wig filename has to also include the strain. This script works on all prot_tables and wigs in the directory. The annotated WIGs are saved as csv files.  
- __MasterTagList.csv__ - A list of gene names and locations from different strains matched up using Roary (3). NCBI tags are also included, matched to the Prokka genes based on location. 
- __permutation.py__ - Adds together the reads from all annotated wig files for a strain (file name must start with Roary_*Strain* and must be a csv) and then for each pair of strains, performs a permutation test.  

### Process
1. Create Prot_Tables (runs on all GFFs in directory) `python make_prot_tables.py`
2. Create TA site files (runs on all fastas in directory) `python findTASites.py`
3. Create wig files (runs on a single fasta-tabular file pair) `python make_wig.py <XXX_TASites.txt> <XXX.tabular>`
4. Run TRANSIT (see <https://pythonhosted.org/tnseq-transit/index.html> for more information).
5. Make sure names are compatible with labelWIGs.py and permutation.py (see file descriptions above). 
6. Rename the genes in the Prot_Tables based on the Roary tags. `python labelWIGs.py` 
7. Perform the permutation test. By comparing the outputs of this script to the TRANSIT outputs, you can determine which differences in gene essentiality identified by TRANSIT are significant.  `python permutation.py`

## Identifying Depleted/Enriched/Upregulated Genes in Treated Files

We compared gene fitness in daptomycin-treated Tn-Seq library samples to untreated samples using Mann-Whitney U tests. We also looked for upregulation signatures using a bootstrapping approach. 

### Files Needed
- __TnSeqTools.py__ - A library of functions needed to run TnSeq_Driver.py. 
- __TnSeq_Driver.py__ - A script that will pair treated and untreated files in a directory (provided they are correctly named) and perform Mann-Whitney U tests to identify cases where gene fitness is significantly different in the treated sample and perform a bootstrapping test to find upregulation signatures. 
- __*.tabular__ - Tabular files created by processing SAM files. See Sequence Data Processing above. 
- __*.fasta__ - A nucleotide fasta for the chromosome of each organism. 
- __*.gff__ - A GFF file output from running Prokka (2) on the fasta files.
- __gff_to_ncbi.py__ - A script that converts from the Prokka gff output to an NCBI-style gene table. Processes all GFFs in a directory at once.  

### Process
1. Make sure all files are correctly named. Fasta files need to be named *Strain*.fasta. Tabular files must be named *Strain*-*GroupID*_*Treatment*_*Transposon*.tabular. Strain must match the fasta. Group ID must match between the appropriate control and the treated files you want to compare it with. Treatment must be "untr" if it is a control file. The recognized transposons are blunt, cap, dual, erm, pen, tuf, combined, and promoters. The gff file must be named *Strain.gff*. 
2. Convert the gff files to NCBI files. `python gff_to_ncbi.py`
3. Run the Mann-Whitney U and upregulation analyses `python TnSeq_Driver.py`

### Outputs
1. __*.igv__ - IGV files that list the reads at each TA dinucleotide site. Tab-delimited. No heading. Columns are 1. Accession number, 2. TA start, 3. TA stop, 4. Total reads, 5. Gene, 6. Plus strand reads, 7. Minus strand reads, 8. Gene strand, 9. Description, 10. Copy (If a TA site appears in multiple genes, it will have multiple entries in the table. 
2. __HitCounts.csv__ - A file that lists the reads in each IGV file, the number of TA sites with reads, and the percent of TA sites with reads. 
3. __MWU*.csv__ - The output of the Mann-Whitney U test. Hits are found by filtering the genes based on the q-value, read ratio, and number of reads (see manuscript for more details). 
4. __UPREG*.csv__ - The output of the Upregulation test. Hits are those with a YES in the Upreg column (see manuscript for more details). 
All other files created are related to mapping and do not contain data. 

## REFERENCES: 
1. DeJesus, Michael A., et al. "TRANSIT-a software tool for Himar1 TnSeq analysis." *PLoS computational biology* 11.10 (2015): e1004401.
2. Seemann, Torsten. "Prokka: rapid prokaryotic genome annotation." *Bioinformatics* 30.14 (2014): 2068-2069. 
3. Andrew J. Page, Carla A. Cummins, Martin Hunt, Vanessa K. Wong, Sandra Reuter, Matthew T. G. Holden, Maria Fookes, Daniel Falush, Jacqueline A. Keane, Julian Parkhill, "Roary: Rapid large-scale prokaryote pan genome analysis", *Bioinformatics*, 2015;31(22):3691-3693 