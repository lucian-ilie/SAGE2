SYSTEM REQUIREMENTS
=============================================================================
SAGE2 is a de novo genome assembler that uses C++/OpenMP.  It is recommended 
that the reads be corrected with RACER before assembling the genome with SAGE2.
SAGE2 is designed to run on a 64-bit Linux environment with gcc installed. 
Human genomes require approximately 250GB - 450GB of memory.

INSTALLATION
=============================================================================

Installing from Anaconda or Bioconda:
	
	conda install -c mmolnar sage2

Installing from GitHub:

	git clone https://github.com/lucian-ilie/SAGE2

Then change to the SAGE2 folder and type:

	make

RUNNING SAGE2
=============================================================================
The reads should be corrected first using RACER which can be installed from Anaconda or Bioconda:

	conda install -c mmolnar racer

Or from GitHub:

	https://github.com/lucian-ilie/RACER

Run RACER using the command:

	racer < inputReads > < correctedReads > < genomeLength >

< inputReads > is the input file containing the reads in FASTA or FASTQ format.  
< correctedReads > is the file name that will contain the corrected reads.  
< genomeLength > is the approximate length of the sampled genome in base pairs. If only parts of a genome were sequenced, then only the total length of those parts should be used (instead of the length of the total genome).  A precise value is not necessary, an approximation will work well.

SAGE2 can then be run with the corrected reads using the command:

	SAGE2 [options] -f <input file> -k <minimum overlap length>

Full list of options:  

	-f|--fileInput <string>
			Input file in FASTA/FASTQ format (interleaved). Mutually
			exclusive with option -l|--listInput. [Required]

	-l|--listInput <string>
			A list of input files in FASTA/FASTQ format. Mutually exclusive
			with option -f|--fileInput. [Required]

	-k|--minOverlap <int>
			The minimum length of an overlap to be considered for graph
			building. [Required]

	-o|--outputDir <string>
			Path to the output directory. [Default: null]

	-p|--prefix <string>
			Prefix of all output files in the output directory.
			[Default: "untitled"]

	-i|--inputPrefix <string>
			Prefix of all input files in the output directory.
			[Default: value of prefix]

	-m|--minStep <int>
			The first step which will be run by the program. [Default: 1]

	-M|--maxStep <int>
			The last step which will be run by the program. [Default: 7]

	-s|--saveAll
			Saves the results of all intermediate steps. Required if starting
			again from a step other than step 1.

	-d|--debug
			Outputs and saves more results for debugging.

	-h|--help
			Prints the usage and full list of options.


INPUT DATASET
=============================================================================
SAGE2 accepts FASTA and FASTQ input files. Paired-end reads should be placed 
one after another (interlaced) in the input file.  SAGE2 does not use single end reads. If using an input list of 
read files the forward and reverse files for paired-end reads must be followed
one after the other in the list with the forward reads first and then the 
reverse reads next.  The files must be listed as f1 and f2:

f1=forwardReads.fastq  
f2=reverseReads.fastq  

Multiple paired-end read files can be listed and interlaced paired-end read 
files can be used with the letter f:

f1=forwardReads1.fastq  
f2=reverseReads1.fastq  
f1=forwardReads2.fastq  
f2=reverseReads2.fastq  
f1=forwardReads3.fastq  
f2=reverseReads3.fastq  
f=interlacedReads1.fastq  
f=interlacedReads2.fastq  
f=interlacedReads3.fastq  
f=interlacedReads4.fastq

OUTPUT FILES
=============================================================================
SAGE2 produces several output files in the directory when using the -s option.
This option is required if you plan on rerunning the assembly from a minimum
step other than step 1. SAGE2 outputs a contig file and a scaffold file.  The 
contig file is named prefix_contig.fasta and the scaffold file is named 
prefix_scaffold.fasta.
