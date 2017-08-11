#!/bin/bash
##pipeline

#---defining usage-----

usage()
{
echo -e "\nUsage: $0:" >&2
echo -e "$0 [-a <reads_file1>][-b <reads_file2>] [-c <dreads_file1>] [-d <dreads_file2>][-r <reference_file>] [-g <reference_gfffile>][-o <final>][-h]" >&2
echo -e "options:
\n\t-a :Necessary: specify input reads file 1 for parent
\n\t-b :Necessary:specify input reads file 2 for parent 
\n\t-c :Necessary: specify input reads file 1 for daughter
\n\t-d :Necessary:specify input reads file 2 for daughter
\n\t-r :Necessary:specify reference genome filename
\n\t-g :Necessary:specify reference GFF filename
\n\t-o :the output vcf file 
\n\t-h : Help " >&2
exit 1
}

#checking if no parameter is specified

#--------------------------------------

if [ $# == 0 ]; then

echo "Error! No parameter found" >&2

usage

fi

# Receiving command line argument and options
#------------------------------------------------------------------------------------
#fetching input options

verbose=false
while getopts ":a:b:c:d:r:g:o:h" option; 
do
## Getting arguments from command line
case "$option" in
a)
reads_file1=$OPTARG
echo "input read file 1 for parent: $reads_file1";;
b)
reads_file2=$OPTARG
echo "input read file 2 for parent: $OPTARG";;
c)
dreads_file1=$OPTARG
echo "input read file 1 for daughter: $OPTARG";;
d)
dreads_file2=$OPTARG
echo "input read file 2 for daughter: $OPTARG";;
r)
reference_file=$OPTARG;;
g)
reference_gfffile=$OPTARG;;
o)
final=$OPTARG
echo "output final: $OPTARG";;  
h)
usage ;;
\?)
echo -e "\nError! Unknown option:" >&2
usage ;;
:) 
echo "option -$OPTARG requires an argument.";;

esac
done

echo "checking for files starts now" 

# Check for existence of files
#-------------------------------------------------------------------------------------------

if [ ! -f $reads_file1 ]; then

echo -e "\nError! $reads_file1 not found or -a option specified improperly" >&2

exit 1

fi

if [ ! -f $reads_file2 ]; then

echo -e "\nError! $reads_file2 not found or -b option specified improperly" >&2

exit 1

fi
if [ ! -f $dreads_file1 ]; then

echo -e "\nError! $dreads_file1 not found or -c option specified improperly" >&2

exit 1

fi
if [ ! -f $dreads_file1 ]; then

echo -e "\nError! $dreads_file1 not found or -d option specified improperly" >&2

exit 1

fi
#check if genome file present in pwd

if [ ! -f $reference_file ]; then

echo -e "\nError! $reference_file not found or -r option specified improperly" >&2

exit 1

fi
#check if VCF file exists if yes overright or exit the program.

if [ -f $final ]; then

#if present ask to overwrite

echo -e "$final already present\n Do you wish to Overwrite(o) or Exit(e)? "
read opinion

#if overwrite permission not given exit 

if [ "$opinion" == "e" ]; then

exit 1

# remove existing file

else

rm -r $final

echo -e "\nOutput file will be overwritten: $final"

fi
fi

#pipeline begins

#------------------------------
#---step 1------#
#-------------Mapping using BWA----------------#
#----indexing------#

echo "indexing reference file";
bwa index $reference_file
echo "indexing finished!";

#----mapping reads to create sam file----#

echo "mapping reads to indexed reference file";
bwa mem -R '@RG\tID:foo\tSM:bar\tLB:library1' $reference_file $reads_file1 $reads_file2 > temp/parent.sam
echo "mapping completed";

#-----cleaning up read pairing information and flags------#
echo "cleaning up read pairing information and flag...";
samtools fixmate -O bam temp/parent.sam temp/parent.bam

#----sorting------#

echo "sorting...";
samtools sort -O bam -o temp/parent_sorted.bam -T /temp/parent_temp temp/parent.bam
echo "sorting done";

#----step 2 --------# 

#-----creating consensus parent sequence--------# 
echo "calling mpileup";
samtools mpileup -E -uf $reference_file temp/parent_sorted.bam -d 8000 | bcftools call -mv -Oz -o temp/parent.vcf.gz
echo "tabix indexing vcf file";
tabix temp/parent.vcf.gz
echo "creating parent consensus file";  
cat $reference_file | bcftools consensus temp/parent.vcf.gz > output/parentconsensus.fasta
echo "consensus file created!";

#----- step 3 ------# 
#-----align the daughter to parent consensus----#
echo "aligning the daughter to parent consensus"

# ---mapping ---#

echo "indexing parent consensus file";
bwa index output/parentconsensus.fasta
echo "indexing finished!";

#----mapping reads to create sam file----#

echo "mapping reads to indexed parent consensus file";
bwa mem -R '@RG\tID:foo\tSM:bar\tLB:library1' $reference_file $dreads_file1 $dreads_file2 > temp/daughter.sam
echo "mapping completed for daughter sample";

#-----cleaning up read pairing information and flags------#
echo "cleaning up read pairing information and flag...";
samtools fixmate -O bam temp/daughter.sam temp/daughter.bam

#----sorting------#

echo "sorting daughter bam file...";
samtools sort -O bam -o temp/daughter_sorted.bam -T temp/daughter_temp temp/daughter.bam
echo "sorting done";

#-----calling varients in daughter wrt parent----#
samtools mpileup -ugf output/parentconsensus.fasta temp/daughter_sorted.bam | bcftools call -vm -o output/daughterparent.vcf

echo "daughter parent vcf file created"
#----step 4-----# 
echo "callin SNPeff";

#-----create database with gff file------#

#----downloading human reference genome to create the data folder-----# 
echo "downloading data folder";

java -Xmx4g -jar ~/snpEff/snpEff.jar download -c ~/snpEff/snpEff.config -v GRCh37.75

echo "data folder created";

#---------creating the new genome folder in data folder ---------#

echo "creating new database";

mkdir ~/snpEff/data/new


#----add the new genome info in config file -----# 
#--------adding the required files for creating database from GFF file -------#
#----copy the GFF and Reference files to the ~/snpEff/data/ecoli -----#
echo "copying genome and gff files";

cp $reference_file temp/
cp $reference_gfffile temp/
cd temp
ref=$(basename "$reference_file" )
refgff=$(basename "$reference_gfffile")

echo "renaming and compressing"; 

mv $refgff genes.gff
mv $ref sequences.fa
#-----If GFF files donot have sequence information -------# 
echo "adding sequence info to gff file";
echo "###" >> genes.gff
echo "##FASTA" >> genes.gff
cat sequences.fa >> genes.gff

cp genes.gff ~/snpEff/data/new
cp sequences.fa ~/snpEff/data/new
cd .. 

#------------------------------------##
gzip -r ~/snpEff/data/new

echo "creating database"; 

java -jar ~/snpEff/snpEff.jar build -gff3 -v new 

echo "database created";  

#-----annotating daughter vcf file generated ------#

echo "annotating vcf file"; 

java -Xmx4g -jar ~/snpEff/snpEff.jar -v new output/daughterparent.vcf > $final

echo "the annotated files are here!!!!!";

#####----------------------------------------------------#####
#####----------------------------------------------------#####
