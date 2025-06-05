# RNASeq-DESeq-Pipeline
Pipeline for basic RNASeq and DESeq, created as part of Saha Lab

Introduction and File Calls 
RNA-Seq Data Processing Protocol
Matthew Luchs - 250225 - Good Luck!

This document outlines a complete RNA-Seq data processing pipeline consisting of sequential scripts for trimming raw reads, alignment to a reference genome, counting aligned reads, and converting the count data to CSV format. It also includes a script for initial DESeq and processes for Functional Analysis. 

Overview of Pipeline Steps

1. Trimming: Process raw FASTQ files to remove adapters and low-quality sequences
2. Alignment: Align trimmed reads to a reference genome
3. Feature Counting: Count reads mapping to genes
4. (Optional) CSV Conversion: Convert the count data to CSV format for downstream analysis

Usage Instructions AND Example File Calls

Step 1: Trimming
1. Run the script: `python 1_Trimming.py`
2. Enter the path to the folder containing your raw FASTQ files
3. Enter the output directory for trimmed files
4. File call Example:
-	python3 /Path/to/your/code/Trimming.py
-	Input: /Path/to/your/FASTQ_Folder
-	Output: /Path/to/where/you/want/trimmed/files/TRIMMED

Step 2: Alignment
1. Run the script: `python 2_Alignment.py`
2. Enter the path to the folder containing your trimmed FASTQ files
3. Enter the output directory for aligned BAM files
File call Example:
-	python3 /Path/to/your/code/Alignment.py
-	input: /Path/containing/trimmed/files/TRIMMED
-	output: /Path/to/where/you/want/aligned/files/ALIGNED


Step 3: Feature Counting
1. Run the script: `python 3_FeatureCounts.py`
2. Enter the path to the folder containing your BAM files
3. Enter the output directory for count files
4. File call Example
-	python3 /Path/to/your/code/FeatureCounts.py
-	input: /Path/containing/aligned/bam/files/ALIGNED
-	output: /Path/where/you/want/counts/COUNTS


Step 4: CSV Conversion
1. Run the script: `python 4_Counts_CSV.py [input_directory]`
   - Replace `[input_directory]` with the path to the folder containing your count files
2. File call Example
-	One Line: python3 /Path/to/your/code/Counts_CSV.py /Where/you/have/count/txt/files/TXT
-	Will simply output in the input folder


 
IMPORTANT NOTES 
*****Notes and Considerations*****

Any Issues you have will likely have to do with one of these issues. Hardcoded paths and file naming conventions are number one consideration, as it was written in the script. Easy fix to make but must be thorough. 

1. Dependencies:
   - Trim Galore (and its dependencies: Cutadapt and FastQC)
   - HISAT2
   - SAMtools
   - Subread (for featureCounts)
   - Python packages: pandas, pathlib

2. Hardcoded Paths:
   - Reference genome index: `/Users/miseq/Desktop/Luchs/Reference_Genome2.0/genome/genome_index`
   - GTF annotation file: 
`/Volumes/CaChannel/Luchs/Reference_gtf/genomic.gtf`
   - Consider modifying these paths if your files are located elsewhere
- They will be in a different place so please don’t forget to change these in the Alignment and Feature Counts Scripts

3. Resource Usage:
   - The alignment step uses 24 threads (`-p 24` for HISAT2 and `-@ 24` for SAMtools)
   - The feature counting step uses 2 threads (`-T 2`)
   - Adjust the above according to the computers capabilities
   - This all requires a lot of space, make sure your system has required storage or files will be improperly processed

4. File Naming Conventions:
   - Input FASTQ files: `*_R1_001_repair.fastq.gz` and `*_R2_001_repair.fastq.gz`
   - Trimmed FASTQ files: `*_R1_001_val_1.fq.gz` and `*_R2_001_val_2.fq.gz`

-	If your files are not named with this convention, the above step will lead to the script not running. Make sure these match your files. Guaranteed number one issue with these scripts, ran into myself one million times. 
-	
   - BAM files: `{sample_name}.bam`
   - Count files: `{sample_name}_counts.txt` and `{sample_name}_counts.csv`

5. Organization:
-	Name everything with clear call to contents and date
-	There is a lot of stuff being input/output in this process so make it easy for yourself  
Example Outputs 
Trimming 
 
Alignment 
Feature Counts 
CSV
  
1. Trimming 
Step 1: Trimming Raw Reads
The first script (`Trimming.py`) processes raw FASTQ files using Trim Galore to remove adapters and low-quality sequences.


Script Explanation:

1. Function `trim_and_fastqc(folder_path, output_dir)`:
   - Uses `glob` to find all paired-end FASTQ files in the specified folder (files with pattern `*_R1_001_repair.fastq.gz` and `*_R2_001_repair.fastq.gz`)
   - Sorts the files to ensure proper pairing of forward (R1) and reverse (R2) reads
   - Verifies that the number of R1 and R2 files match
   - For each pair of files:
     - Extracts a base name from the R1 file to create output file names
     - Runs Trim Galore with these parameters:
       - `--fastqc`: Generates FastQC quality reports
       - `--paired`: Processes paired-end data
       - `--length 20`: Discards reads shorter than 20 bases
       - `--clip_R2 15`: Trims 15 bases from the 5' end of R2 reads
       - `--three_prime_clip_R1 15`: Trims 15 bases from the 3' end of R1 reads

2. Function `check_trimming()`:
   - Prompts the user for input and output directories
   - Calls the `trim_and_fastqc` function with these directories
   
Expected Inputs and Outputs:

- Inputs: Paired-end FASTQ files (gzipped)
- Outputs: 
  - Trimmed FASTQ files 
  - FastQC quality reports for each file

 
2. Reference Genome 
From Molly Estes RNASeq Processing Protocol, No Script 

REFERENCE GENOME CREATION 
Before you can align your reads to a genome you need to build a reference genome. This will format the genome (fasta format) into a usable form for the Hisat2 align function. 
1. Navigate to your species of interests NCBI page(or anywhere you can find a genome). Here is an example for Zebra Finch 
https://www.ncbi.nlm.nih.gov/datasets/taxonomy/59729/ Then download the reference genome. You will want to download a gtf or gff file which will be used later as well as the fasta with extension .fasta or .fna. 
2. This command will make a reference genome 
“hisat2-build --large-index -p 16 
/Volumes/CaChannel/ZebraFinchBrain/GCF_003957565.2/ncbi_dataset/data/GCF_0039 57565.2/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna 
/Volumes/CaChannel/ZebraFinchBrain/reftwo/genome_index" 
hisat2-build: Invokes HISAT2 to build the index. 
--large-index: This option tells HISAT2 to use more memory for the index building process, potentially allowing for larger genomes to be indexed. 
-p 16: Specifies the number of threads (parallel processes) to use during index building. This time, it's set to 16, meaning the index building process will utilize 16 CPU cores. 
/Volumes/CaChannel/ZebraFinchBrain/GCF_003957565.2/ncbi_dataset/d ata/GCF_003957565.2/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna: Path to the input FASTA file containing the genomic sequence data. INPUT 
/Volumes/CaChannel/ZebraFinchBrain/reftwo/genome_index: Specifies the directory and prefix name for the index files that will be generated. HISAT2 will create
several files with this prefix in the specified directory containing the index information. OUTPUT 
3. You may get this error "Assertion failed: (0), function HGFM, file ./hgfm.h, line 2127 Abort trap: 6". This can be ignored as the genome is most likely successfully created. The output will be multiple files with the refix genome_index and the suffix .ht2l. Run the next command and if it works then the reference genome is built.  
3. Alignment 
Step 2: Alignment to Reference Genome
The second script (`Alignment.py`) aligns the trimmed reads to a reference genome using HISAT2 and converts the output to sorted BAM files using SAMtools.

Script Explanation:

1. Function `automated_alignment(input_folder, output_folder)`:
   - Recursively walks through the input folder and its subfolders
   - Identifies trimmed R1 files (ending with `_R1_001_val_1.fq.gz`)
   - For each R1 file:
     - Extracts the sample name
     - Checks if the corresponding R2 file exists
     - Creates and executes a pipeline command that:

       - Runs HISAT2 alignment with these parameters:
         - `-t`: Print timing information
         - `--rna-strandness RF`: Specifies the strand-specific library type (RF for Illumina TruSeq)
         - `--summary-file`: Saves alignment summary
         - `-p 24`: Uses 24 threads for processing
         - `-x`: Path to the reference genome index

       - Pipes the output to `samtools view` to convert SAM to uncompressed BAM
       - Pipes to `samtools sort` to sort the BAM file by read name (`-n` flag)
       - Outputs a sorted BAM file

2. Function `check_alignment()`:
   - Prompts the user for input and output directories
   - Calls the `automated_alignment` function

Expected Inputs and Outputs:

- Inputs: Trimmed FASTQ files from Step 1 
- Outputs:
  - Name-sorted BAM files (named as `{sample_name}.bam`)
  - Alignment summary text files (named as `{sample_name}.txt`)
 
4. Feature Counts 
Step 3: Feature Counting
The third script (`FeatureCounts.py`) uses the featureCounts tool to count reads mapping to genomic features (genes) from the aligned BAM files.


Script Explanation:

1. Helper Functions:
   - `check_file_exists(filepath, description)`: Verifies that a file or directory exists
   - `ensure_directory_exists(directory)`: Creates a directory if it doesn't exist

2. Function `run_feature_counts(gtf_file, input_folder, output_folder)`:
   - Validates inputs and creates the output directory if needed
   - Finds all BAM files in the input folder
   - For each BAM file:
     - Extracts the sample name from the filename
     - Verifies that the BAM file exists and is not empty

     - Builds a command for featureCounts with these parameters:
       - `-a`: Path to the GTF annotation file
       - `-T 2`: Uses 2 threads
       - `-s 2`: Specifies strand-specific counting mode (reverse)
       - `-Q 0`: No mapping quality filter
       - `-t gene`: Count at the gene level
       - `-g gene_id`: Use gene_id for feature identification
       - `-p`: Paired-end reads
       - `-C`: Exclude chimeric fragments

     - Executes the command and captures stdout/stderr
     - Tracks success/failure status and validates output
   - Prints a summary of successful and failed runs

3. Function `fcounts_pipeline()`:
   - Uses a hardcoded path to the GTF annotation file
   - Prompts the user for input and output directories
   - Calls `run_feature_counts` and handles errors

Expected Inputs and Outputs:

- Inputs:
  - BAM files from Step 2
  - GTF annotation file (hardcoded path): 
- Outputs:
  - Count files in tabular format (named as `{sample_name}_counts.txt`)
  - Summary of successful and failed runs printed to the console



 
5. TXT to CSV 
OPTIONAL Step 4: Converting Count Data to CSV
The final script (`Counts_CSV.py`) converts the featureCounts output files to CSV format for easier downstream analysis. Not necessary, just easier to work with than txt files.

Script Explanation:

1. Function `convert_to_csv(file_path)`:
   - Reads a featureCounts output file, skipping the first comment line
   - Renames the last column (which contains the actual counts) to "Count"
   - Creates a new filename with .csv extension
   - Saves the data as CSV
   - Prints confirmation of the conversion

2. Function `main()`:
   - Uses `argparse` to handle command-line arguments
   - Requires the input directory as an argument
   - Finds all files matching the pattern `*_counts.txt`
   - Calls `convert_to_csv` for each file

Expected Inputs and Outputs:

- Inputs: Count files from Step 3 (named as `{sample_name}_counts.txt`)
- Outputs: CSV-formatted count files (named as `{sample_name}_counts.csv`)
 
6. DESeq 
Overview
These scripts are designed to run differential gene expression analysis using DESeq2. The code consists of two main files:
1.	An R script (deseq2_single_comparison.R) that performs a single differential expression analysis between two groups
2.	A Python wrapper (DESeq2_Comparison_Runner.py) that can run multiple comparisons by calling the R script repeatedly
What the Scripts Do
deseq2_single_comparison.R
●	Takes RNA-seq count data from two groups (conditions)
●	Normalizes the data using DESeq2's methods
●	Identifies differentially expressed genes (DEGs) between the two groups
●	Generates several visualizations (MA plot, PCA plot, p-value histogram)
●	Creates comprehensive output files with the results
Key features:
●	Reads count data from CSV files
●	Performs filtering of low-expression genes
●	Calculates statistical significance with adjusted p-values
●	Creates data visualizations for quality control and results interpretation
The Python Wrapper Script
●	Manages multiple DESeq2 comparisons in a single run
●	Can run comparisons in parallel to save time
●	Tracks the progress and success of each comparison
●	Generates a summary HTML report of all the analyses
●	Provides detailed logging




How to Use
Basic Usage
1.	Make sure you have R and Python installed with the required packages:
○	R: DESeq2, edgeR, ggplot2, ggrepel
○	Python: standard libraries should be sufficient
2.	Organize your data:
○	Put count files (CSV format) in separate directories for each group/condition
○	Count files should contain gene IDs and counts
3.	Run the Python wrapper script (Here is what I used): 
python3 /Users/miseq/Desktop/Luchs/Code/DESeq/python-runner.py --base-dir /Users/miseq/Desktop/Luchs/ForDESeq --output-dir /Users/miseq/Desktop/Luchs/Comparisons --r-script /Users/miseq/Desktop/Luchs/Code/DESeq/DESeqSinglev2.R
Command Line Arguments
●	--base-dir: Root directory containing all your sample group directories
●	--output-dir: Where to save all results
●	--r-script: Path to the R script
Optional: 
●	--parallel: Number of comparisons to run simultaneously (default: 1)
●	--selected: Optional list of specific comparison names to run (instead of all)
●	--log-file: Where to save the detailed log
●	Output
The scripts will generate:
●	CSV files with differential expression results
●	PDF visualizations (MA plots, PCA plots, p-value histograms)
●	Summary text files with key statistics
●	An HTML summary report that links to all individual results
●	Detailed log files for troubleshooting


7. Functional Analysis 
FUNCTIONAL ANALYSIS (From Bjorn Shockey)

Congratulations! You are almost done with your RNA Sequencing! It will soon be time to start writing the paper! The first step is to see which pathways your differential expressed genes are related to. We will use DAVID’s Functional Analysis Webtool to do this. 
1. As part of the above code you will receive a csv that has the lists of DE genes, the counts associated with each RNA file, the Log2FoldChange(Upregulated or Downregulated), and a padj value(is your gene statistically significantly DE). 
 2. Highlight the gene column (not including the word or prefix gene. Use Excel command to delete the prefixes here) and go to https://delim.co/ to separate your gene names by commas.  
3. Next go to https://david.ncifcrf.gov/tools.jsp and click start analysis then upload.  4. Now enter the genes in the paste a list block. Select OFFICIAL_GENE_SYMBOL unless your list is in some other form(unlikely). Then select your species, in my case the Zebra Finch is called Taeniopygia Guttata. If you are using E. Coli you will most likely select “Escherichia coli str. K-12 substr. MG1655” unless you are using a different strain. Select
‘gene list’. Then Submit list. 
 5. Select Functional Analysis Chart. This will give you the list of terms, with the number of genes associated with each (count) and the benjamini score (essentiall an adjusted pvalue of significance). You will most likely be interested in analyzing the statistically significant terms Benjamini < 0.05. 
 
6. Note: at the bottom it will say “ x gene(s) from your list are not in the output.” This means that for these genes they were either not recognized or the pathways they are
associated with are not well represented enough to be in the output. The minimum count settings and pval can be edited in the options + symbol at the top of the page to get these genes in your output. Now you may celebrate completing RNA Seq and write the paper.


