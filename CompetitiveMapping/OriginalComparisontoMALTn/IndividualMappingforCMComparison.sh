#!/bin/bash
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=12:00:00
 #PBS -l pmem=6gb
 #PBS -A open
 #PBS -o /gpfs/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/IndividualAlignments
 #PBS -e /gpfs/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/IndividualAlignments

### FEB 2023 ###
# Individual mapping and table for comparison to competitive mapping (species did not yet have data for/same mapping parameters).
# Calculates analysis ready reads, total mapped, total Q30 mapped, total unique, average length of mapped reads, breadth of coverage, and depth of coverage and puts results in a summary table.

##########################################

# Reference (indexed with bwa and samtools)
ref=/gpfs/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/IndividualAlignments/References/Streptococcus_gordonii/Streptococcus_gordonii_NCTC10231

# Path to directory with Samples (each in their own directory)
PATHS=/gpfs/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/IndividualAlignments/Samples

# Species mapping to to name summary file, name out directory
species=Streptococcus_gordonii
outdir=Streptococcus_gordonii

##############################################################################

# Load modules
module load gcc bwa samtools

# Prepare the table
echo -e "Sample \t Analysis-ready reads \t Total mapped reads \t Total Q30 mapped reads \t  Total after removing PCR duplicates \t Total unique reads \t Average read length \t Total SNPs \t Breadth of coverage \t Depth of coverage" > $PATHS/../"${species}"_MappingResults.txt

# Start of loop
for Sample in $PATHS/Sample-*;
 do
  cd $Sample

  ID=$(basename $Sample)
  NAME=$(echo $ID |cut -d "-" -f2)
  echo "Running mapping on "${NAME}"..."
  
  mkdir $outdir
  
###########################################

# Alignment with bwa aln
  bwa aln -l 1024 -n 0.01 -o 2 -t 8 $ref $PATHS/Sample-"${NAME}"/*fastq.gz > $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_firstalignment.sai

# Alignment in sam file format
  bwa samse $ref $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_firstalignment.sai $PATHS/Sample-"${NAME}"/"${NAME}"*fastq.gz > $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_firstalignment.sam
  
###########################################
# Filtering, sorting, removing duplicates and non-unique reads using SAMtools

# SAM to BAM, keep header
  samtools view -bSh $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_firstalignment.sam > $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_firstalignment.bam

# Filter out unmappped and low quality reads
  samtools view -bh -F4 $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_firstalignment.bam > $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_mapped.bam
  
# Sort mapped alignments by leftmost coordinates
  samtools sort -o $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_mapped_sorted.bam $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_mapped.bam

# Filter out by mapping quality greater than 30
  samtools view -bh -q 30  $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_mapped_sorted.bam > $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_mapped_q30_sorted.bam

# Remove PCR duplicates
  samtools rmdup -s $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_mapped_q30_sorted.bam $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_mapped_q30_sorted_rmdup.bam
  
# Remove reads with multiple mappings
  samtools view -h $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_mapped_q30_sorted_rmdup.bam |grep -v 'XT:A:R'|grep -v 'XA:Z' |grep -v 'XT:A:M' |awk '{if($0~/X1:i:0/||$0~/^@/)print $0}' | samtools view -bS - > $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_mapped_q30_sorted_rmdup_uniq.bam

###########################################
# Gather and print alignment statistics

# Create a report summary file
  echo " " > $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_MappingReport.txt
  echo "${NAME}" >> $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_MappingReport.txt
  echo "---------------------------------------" >> $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_MappingReport.txt

# Print general data
# Analysis-ready reads
 echo "Analysis-ready reads" >> $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_MappingReport.txt
  samtools view -c $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_firstalignment.bam >> $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_MappingReport.txt
# Total mapped reads
  echo "Total mapped reads" >> $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_MappingReport.txt
  samtools view -c $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_mapped.bam >> $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_MappingReport.txt
# Total Q30 mapped reads
  echo "Total Q30 mapped reads" >> $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_MappingReport.txt
  samtools view -c $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_mapped_q30_sorted.bam >> $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_MappingReport.txt
# Total after removing duplicates
  echo "Total after removing PCR duplicates" >> $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_MappingReport.txt
  samtools view -c $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_mapped_q30_sorted_rmdup.bam >> $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_MappingReport.txt
# Total after removing non-unique reads
echo "Total unique reads" >> $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_MappingReport.txt
  samtools view -c $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_mapped_q30_sorted_rmdup_uniq.bam >> $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_MappingReport.txt
# Average length of reads
  echo "Average length of reads" >> $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_MappingReport.txt
  samtools view $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_mapped_q30_sorted_rmdup_uniq.bam |awk '{SUM+=length($10);DIV++}END{print SUM/DIV}' >> $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_MappingReport.txt
  
# Bring to variables for later printing in summary table
AnalysisReady=`grep -A 1 "Analysis-ready reads" $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_MappingReport.txt | grep -v "Analysis-ready reads"`
TotalMapped=`grep -A 1 "Total mapped reads" $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_MappingReport.txt | grep -v "Total mapped reads"`
TotalQ30=`grep -A 1 "Total Q30 mapped reads" $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_MappingReport.txt | grep -v "Total Q30 mapped reads"`
TotalNoDuplicates=`grep -A 1 "Total after removing PCR duplicates" $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_MappingReport.txt | grep -v "Total after removing PCR duplicates"`
TotalUnique=`grep -A 1 "Total unique reads" $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_MappingReport.txt | grep -v "Total unique reads"`
AvgLength=`grep -A 1 "Average length of reads" $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_MappingReport.txt | grep -v "Average length of reads"`

# Calculate and print coverage statistics
# Use SAMtools Coverage
samtools coverage $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_mapped_q30_sorted_rmdup_uniq.bam -o $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_Coverage.txt
# For breadth of coverage
covbases=`grep -v "name" $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_Coverage.txt | cut -f5 | paste -sd+ | bc `
numberofbases=`grep -v "name" $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_Coverage.txt | cut -f3 | paste -sd+ | bc `
breadthofcoverage=`echo $covbases / $numberofbases | bc -l`
# For depth of coverage
totalmappinglength=`echo $AvgLength \* $TotalUnique | bc -l`
depthofcoverage=`echo $totalmappinglength / $numberofbases | bc -l`
# Printing
# Breadth of coverage
  echo "Breadth of coverage" >> $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_MappingReport.txt
  echo "$breadthofcoverage" >> $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_MappingReport.txt
# Depth of coverage
  echo "Depth of coverage" >> $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_MappingReport.txt
  echo "$depthofcoverage" >> $PATHS/Sample-"${NAME}"/$outdir/"${NAME}"_MappingReport.txt

###########################################
# Print all variables to tables
echo -e "${NAME} \t $AnalysisReady \t $TotalMapped \t $TotalQ30 \t $TotalNoDuplicates \t $TotalUnique \t $AvgLength \t $breadthofcoverage \t $depthofcoverage \t" >> $PATHS/../"${species}"_MappingResults.txt

done