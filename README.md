# 5SATnSeq
Scripts and files needed to reproduce the analyses described in !!!paper citation!!!. 

Note: All python scripts were written for python 2.7 and require some or all of the following libraries: glob, numpy, pandas

## Sequence Data Processing

FASTQ files obtained from sequencing facility are processed using Galaxy and in-house python scripts to ultimately obtain tabular files that map Tn sequencing reads to individual sequence sites. 


### Files Needed
- __*.FASTQ__ - The raw sequencing reads. Our FASTQs can be downloaded from the Sequence Reads Archive, accession ID ~. 
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

We identified the essential genes using the TRANSIT software Gumbel method (1). The scripts below are those used to convert files so that they could be used with the TRANSIT software and to perform a permutation test to determine which genes have significant fitness differences across strains of S. aureus. 

### Files Needed
- __*.tabular__ - The tabular files acquired from processing the FASTQ files (see above)
- __*.fasta__ - A chromosome nucleotide fasta for the genome you are mapping to. See paper methods for the NCBI accession numbers for the strains used.
- __*.gff__ - A GFF file output from running Prokka (2) on the fasta files.
- __make_prot_tables.py__ - A python script that converts the GFF files from Prokka into prot_table files for TRANSIT. 
- __findTASites.py__ - For each fasta file in the directory, makes a list of all of the TA sites available.
- __make_wig.py__ - Makes a wig file recognized by the TRANSIT software from a tabular file and a list of TA sites created using findTASites.py.

### Process
1. Create Prot_Tables (runs on all GFFs in directory) `python make_prot_tables.py`
2. Create TA site files (runs on all fastas in directory) `python findTASites.py`
3. Create wig files (runs on a single fasta-tabular file pair) `python make_wig.py <XXX_TASites.txt> <XXX.tabular>`
4. Run TRANSIT (see <https://pythonhosted.org/tnseq-transit/index.html> for more information). 

REFERENCES: 
1. DeJesus, Michael A., et al. "TRANSIT-a software tool for Himar1 TnSeq analysis." PLoS computational biology 11.10 (2015): e1004401.
2. Seemann, Torsten. "Prokka: rapid prokaryotic genome annotation." Bioinformatics 30.14 (2014): 2068-2069. 