#___________________________________________________________________________________________________________________

jmi.sh : SNP Annotation Pipeline
#___________________________________________________________________________________________________________________

./jmi.sh -a data/sample-1_S1_L001_R1_001.fastq -b data/sample-1_S1_L001_R2_001.fastq -c data/sample-2_S6_L001_R1_001.fastq -d data/sample-2_S6_L001_R2_001.fastq 
-r data/GCF_000005845.2_ASM584v2_genomic.fna -g data/GCF_000005845.2_ASM584v2_genomic.fna -o ann.vcf
All input files must be in data folder.

-a <reads_file1>: Parent Input reads file pair 1 (Required)
-b <reads_file2>: Parent Input reads file pair 2 (Required)
-c <reads_file1>: Daughter Input reads file pair 1 (Required)
-d <reads_file2>: Daughter Input reads file pair 2 (Required)
-r <reference_file>:    Input reads reference genome file (Required)
-g <reference_gfffile> : Input reads reference GFF file 
-o <final> :    Output vcf file (Required)
-h  Print usage information

#___________________________________________________________________________________________________________________

Structure of pipeline directories:

|------ jmi.sh
	|------ data/ [ALL input files must be in a data folder]
		|------ Parent R1 file
		|------ parent R2 file
    	|------ daughter R1 file
    	|------ daughter R2 file
		|------ Reference file
		|------ Reference GFF file
	|------ temp/ [where all temp file(s) are stored]
	|------ output/ [where output files are stored]
		

#_____________________________________________________________________________________________________________________

Output file of jmi.sh is a VCF file of annotated daughter SNPs which is aligned to Parent Consensus Sequence.
VCF file has information about SNP,INDEL and structural variations in the input file in correspondence with Parent Consensus genome.

Steps followed:
1. Aligning the Parent paired end reads with Reference genome.This is done using bwa tools
2. Sorting the bam file using Samtools
3. Assembling the Parent Consensus sequence using samtools and bcftools. 
4. Aligning the Daughter paired end reads with Parent consensus sequence. This is done by bwa tools
5. Calling varients from aligned daughter sequence with respect to Parent sequence
6. Variant calling using Samtools.
7. Annotation of the daughter SNP file using SnpEff and GFF file.

Importants steps to be completed before running the pipeline: 

1. Create a directory with all the mentioned directories 
2. Place all the input files in the data directory accordingly
3. Install SnpEff in the home directory (~)
4. Open the snpEff.config file and add the genome of interest as "new" and save the file
5. Go back to the directory where all the files are present 
6. Run jmi.sh with the command lines as given 

