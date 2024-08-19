#!/bin/bash
 #PBS -l nodes=1:ppn=1
 #PBS -l walltime=40:00:00
 #PBS -l pmem=26gb
 #PBS -A open
 #PBS -o /gpfs/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping
 #PBS -e /gpfs/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping
 
## NOV 2022 ##
# Based very heavily on and including sections from the mapping portion of the mapDamage script by Sterling Wright, this script performs competitive mapping against Streptococcus species and produces a mapping report file for each sample in the sample folder and a tab-delimited file with results from all samples.
 
## Prepare the concatenated reference ##
# Change the reference name in each FASTA file (including for all contigs) to include species name before it (ex: S.oralisNZ_LR134336.1), put FASTAs to be joined in the same directory, and use: cat Streptococcus* > combinedreference
# Index the reference with BWA

## Define paths ##
# Path to directory with indexed reference
ref=/gpfs/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/JoinedReferences/combinedreference

# Path to directory with Samples (each in their own directory)
PATHS=/gpfs/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples

# Load modules
module load gcc bwa samtools

##############################################################################

# Prepare the table
#echo -e "Sample \t Analysis-ready reads \t Total mapped reads \t Total Q30 mapped reads \t  Total after removing PCR duplicates \t Total unique reads \t Nonunique Streptococcus cristatus \t Nonunique Streptococcus gordonii \t Nonunique Streptococcus mitis \t Nonunique Streptococcus oralis \t Nonunique Streptococcus sanguinis \t Nonunique Streptococcus sp.DD04 \t Unique Streptococcus cristatus \t Unique Streptococcus gordonii \t Unique Streptococcus mitis \t Unique Streptococcus oralis \t Unique Streptococcus sanguinis \t Unique Streptococcus sp.DD04 \t" > $PATHS/CollectedCompetitiveMappingResults.txt

# Start of loop
for Sample in $PATHS/Sample-*;
 do
  cd $Sample

  ID=$(basename $Sample)
  NAME=$(echo $ID |cut -d "-" -f2)
  echo "Running competitve mapping on "${NAME}"..."
  
###########################################

# Alignment with bwa aln
  bwa aln -l 1024 -n 0.01 -o 2 -t 8 $ref $PATHS/Sample-"${NAME}"/*fastq.gz > $PATHS/Sample-"${NAME}"/"${NAME}"_firstalignment.sai

# Alignment in sam file format
  bwa samse $ref $PATHS/Sample-"${NAME}"/"${NAME}"_firstalignment.sai $PATHS/Sample-"${NAME}"/"${NAME}"*fastq.gz > $PATHS/Sample-"${NAME}"/"${NAME}"_firstalignment.sam
  
###########################################

# Filtering, sorting, removing duplicates and non-unique reads using SAMtools

# SAM to BAM, keep header
  samtools view -bSh $PATHS/Sample-"${NAME}"/"${NAME}"_firstalignment.sam > $PATHS/Sample-"${NAME}"/"${NAME}"_firstalignment.bam

# Filter out unmappped and low quality reads
  samtools view -bh -F4 $PATHS/Sample-"${NAME}"/"${NAME}"_firstalignment.bam > $PATHS/Sample-"${NAME}"/"${NAME}"_mapped.bam
  
# Sort mapped alignments by leftmost coordinates
  samtools sort -o $PATHS/Sample-"${NAME}"/"${NAME}"_mapped_sorted.bam $PATHS/Sample-"${NAME}"/"${NAME}"_mapped.bam

  samtools view -bh -q 30  $PATHS/Sample-"${NAME}"/"${NAME}"_mapped_sorted.bam > $PATHS/Sample-"${NAME}"/"${NAME}"_mapped_q30_sorted.bam

# Remove PCR duplicate
  samtools rmdup -s $PATHS/Sample-"${NAME}"/"${NAME}"_mapped_q30_sorted.bam $PATHS/Sample-"${NAME}"/"${NAME}"_mapped_q30_sorted_rmdup.bam
  
## Remove reads with multiple mappings
## grep -v means invert the match, it returns all non matching lines (all the reads that are not multiple mappings)
## [XT:A:U flag in the SAM file denotes unique read, and XT:A:R and XA:Z denote multiple mappings for that read]
## When working with paired end reads, it may also be useful to consider filtering out reads with the flag XT:A:M
## (one-mate recovered) which means that one of the pairs is uniquely mapped and the other isn't]
## Use awk to scan input file for lines that match the pattern: {if($0~/X1:i:0/||$0~/^@/)print $0}
## X1= Number of suboptimal hits found by BWA, if X1=0 then keep that read in file
  samtools view -h $PATHS/Sample-"${NAME}"/"${NAME}"_mapped_q30_sorted_rmdup.bam |grep -v 'XT:A:R'|grep -v 'XA:Z' |grep -v 'XT:A:M' |awk '{if($0~/X1:i:0/||$0~/^@/)print $0}' | samtools view -bS - > $PATHS/Sample-"${NAME}"/"${NAME}"_mapped_q30_sorted_rmdup_uniq.bam

###########################################

# Print general alignment statistics

# Create a report summary file
  echo " " > $PATHS/Sample-"${NAME}"/"${NAME}"_CompetitiveMappingReport.txt
  echo "${NAME}" >> $PATHS/Sample-"${NAME}"/"${NAME}"_CompetitiveMappingReport.txt
  echo "---------------------------------------" >> $PATHS/Sample-"${NAME}"/"${NAME}"_CompetitiveMappingReport.txt

# Analysis-ready reads
 echo "Analysis-ready reads" >> $PATHS/Sample-"${NAME}"/"${NAME}"_CompetitiveMappingReport.txt
  samtools view -c $PATHS/Sample-"${NAME}"/"${NAME}"_firstalignment.bam >> $PATHS/Sample-"${NAME}"/"${NAME}"_CompetitiveMappingReport.txt
# Total mapped reads
  echo "Total mapped reads" >> $PATHS/Sample-"${NAME}"/"${NAME}"_CompetitiveMappingReport.txt
  samtools view -c $PATHS/Sample-"${NAME}"/"${NAME}"_mapped.bam >> $PATHS/Sample-"${NAME}"/"${NAME}"_CompetitiveMappingReport.txt
# Total Q30 mapped reads
  echo "Total Q30 mapped reads" >> $PATHS/Sample-"${NAME}"/"${NAME}"_CompetitiveMappingReport.txt
  samtools view -c $PATHS/Sample-"${NAME}"/"${NAME}"_mapped_q30_sorted.bam >> $PATHS/Sample-"${NAME}"/"${NAME}"_CompetitiveMappingReport.txt
# Total after removing duplicates
  echo "Total after removing PCR duplicates" >> $PATHS/Sample-"${NAME}"/"${NAME}"_CompetitiveMappingReport.txt
  samtools view -c $PATHS/Sample-"${NAME}"/"${NAME}"_mapped_q30_sorted_rmdup.bam >> $PATHS/Sample-"${NAME}"/"${NAME}"_CompetitiveMappingReport.txt
# Total after removing non-unique reads
echo "Total unique reads" >> $PATHS/Sample-"${NAME}"/"${NAME}"_CompetitiveMappingReport.txt
  samtools view -c $PATHS/Sample-"${NAME}"/"${NAME}"_mapped_q30_sorted_rmdup_uniq.bam >> $PATHS/Sample-"${NAME}"/"${NAME}"_CompetitiveMappingReport.txt
  
###########################################

# Calculate and print statistics by species

# Index SAM files to use in SAMtools idxstats
samtools index $PATHS/Sample-"${NAME}"/"${NAME}"_mapped_sorted.bam $PATHS/Sample-"${NAME}"/"${NAME}"_mapped_sorted.bai
samtools index $PATHS/Sample-"${NAME}"/"${NAME}"_mapped_q30_sorted_rmdup.bam $PATHS/Sample-"${NAME}"/"${NAME}"_mapped_q30_sorted_rmdup.bai
samtools index $PATHS/Sample-"${NAME}"/"${NAME}"_mapped_q30_sorted_rmdup_uniq.bam $PATHS/Sample-"${NAME}"/"${NAME}"_mapped_q30_sorted_rmdup_uniq.bai

# Use SAMtools index to print results based on reference. Columns are reference sequence name, sequence length, number of mapped reads, and number of unmapped reads (which would be 0 because unmapped reads were already filtered out). The -f 1,3 prints only the first and third columns.
# An extra step is added for Streptococcus sp. DD04, Streptococcus sanguinis, and Streptococcus mitis because the references are not a complete genome. Reads mapped to each contig are added together for the CompetitiveMappingReport but the un-summed data remains in the secondary file.
# Total reads mapped to species
echo "---------------------------------------" >> $PATHS/Sample-"${NAME}"/"${NAME}"_CompetitiveMappingReport.txt
echo "Total reads mapped to species" >> $PATHS/Sample-"${NAME}"/"${NAME}"_CompetitiveMappingReport.txt
samtools idxstats $PATHS/Sample-"${NAME}"/"${NAME}"_mapped_sorted.bam | cut -f 1,3 >> $PATHS/Sample-"${NAME}"/"${NAME}"_Totalidxstats.txt

echo -e "S.spDD04AllContigs \t `grep "S.sp.DD04" $PATHS/Sample-"${NAME}"/"${NAME}"_Totalidxstats.txt | cut -f2 | paste -sd+ | bc`" >> $PATHS/Sample-"${NAME}"/"${NAME}"_Totalidxstats.txt
echo -e "S.mitisAllContigs \t `grep "S.mitis" $PATHS/Sample-"${NAME}"/"${NAME}"_Totalidxstats.txt | cut -f2 | paste -sd+ | bc`" >> $PATHS/Sample-"${NAME}"/"${NAME}"_Totalidxstats.txt
echo -e "S.sanguinisAllContigs \t `grep "S.sanguinis" $PATHS/Sample-"${NAME}"/"${NAME}"_Totalidxstats.txt | cut -f2 | paste -sd+ | bc`" >> $PATHS/Sample-"${NAME}"/"${NAME}"_Totalidxstats.txt

grep -v "*" $PATHS/Sample-"${NAME}"/"${NAME}"_Totalidxstats.txt | grep -v "S.mitisNZ" | grep -v "S.sanguinisNZ" | grep -v "NZ_KQ9" >> $PATHS/Sample-"${NAME}"/"${NAME}"_CompetitiveMappingReport.txt


# Non-unique reads mapped to species (Q30+, PCR duplicates removed)
echo "---------------------------------------" >> $PATHS/Sample-"${NAME}"/"${NAME}"_CompetitiveMappingReport.txt
echo "Non-unique reads mapped to species (Q30 and PCR duplicates removed)" >> $PATHS/Sample-"${NAME}"/"${NAME}"_CompetitiveMappingReport.txt
samtools idxstats $PATHS/Sample-"${NAME}"/"${NAME}"_mapped_q30_sorted_rmdup.bam | cut -f 1,3 >> $PATHS/Sample-"${NAME}"/"${NAME}"_NUidxstats.txt

echo -e "S.spDD04AllContigs \t `grep "S.sp.DD04" $PATHS/Sample-"${NAME}"/"${NAME}"_NUidxstats.txt | cut -f2 | paste -sd+ | bc`" >> $PATHS/Sample-"${NAME}"/"${NAME}"_NUidxstats.txt
echo -e "S.mitisAllContigs \t `grep "S.mitis" $PATHS/Sample-"${NAME}"/"${NAME}"_NUidxstats.txt | cut -f2 | paste -sd+ | bc`" >> $PATHS/Sample-"${NAME}"/"${NAME}"_NUidxstats.txt
echo -e "S.sanguinisAllContigs \t `grep "S.sanguinis" $PATHS/Sample-"${NAME}"/"${NAME}"_NUidxstats.txt | cut -f2 | paste -sd+ | bc`" >> $PATHS/Sample-"${NAME}"/"${NAME}"_NUidxstats.txt

grep -v "*" $PATHS/Sample-"${NAME}"/"${NAME}"_NUidxstats.txt | grep -v "S.mitisNZ" | grep -v "S.sanguinisNZ" | grep -v "NZ_KQ9" >> $PATHS/Sample-"${NAME}"/"${NAME}"_CompetitiveMappingReport.txt


# Unique reads mapped to species
echo "---------------------------------------" >> $PATHS/Sample-"${NAME}"/"${NAME}"_CompetitiveMappingReport.txt
echo "Unique reads mapped to species" >> $PATHS/Sample-"${NAME}"/"${NAME}"_CompetitiveMappingReport.txt
samtools idxstats $PATHS/Sample-"${NAME}"/"${NAME}"_mapped_q30_sorted_rmdup_uniq.bam | cut -f 1,3 >> $PATHS/Sample-"${NAME}"/"${NAME}"_Uidxstats.txt

echo -e "S.spDD04AllContigs \t `grep "S.sp.DD04" $PATHS/Sample-"${NAME}"/"${NAME}"_Uidxstats.txt | cut -f2 | paste -sd+ | bc`" >> $PATHS/Sample-"${NAME}"/"${NAME}"_Uidxstats.txt
echo -e "S.mitisAllContigs \t `grep "S.mitis" $PATHS/Sample-"${NAME}"/"${NAME}"_Uidxstats.txt | cut -f2 | paste -sd+ | bc`" >> $PATHS/Sample-"${NAME}"/"${NAME}"_Uidxstats.txt
echo -e "S.sanguinisAllContigs \t `grep "S.sanguinis" $PATHS/Sample-"${NAME}"/"${NAME}"_Uidxstats.txt | cut -f2 | paste -sd+ | bc`" >> $PATHS/Sample-"${NAME}"/"${NAME}"_Uidxstats.txt

grep -v "*" $PATHS/Sample-"${NAME}"/"${NAME}"_Uidxstats.txt | grep -v "S.mitisNZ" | grep -v "S.sanguinisNZ" | grep -v "NZ_KQ9" >> $PATHS/Sample-"${NAME}"/"${NAME}"_CompetitiveMappingReport.txt


###########################################

# Make tabbed file for all results

# Collect general data
AnalysisReady=`grep -A 1 "Analysis-ready reads" $PATHS/Sample-"${NAME}"/"${NAME}"_CompetitiveMappingReport.txt | grep -v "Analysis-ready reads"`
TotalMapped=`grep -A 1 "Total mapped reads" $PATHS/Sample-"${NAME}"/"${NAME}"_CompetitiveMappingReport.txt | grep -v "Total mapped reads"`
TotalQ30=`grep -A 1 "Total Q30 mapped reads" $PATHS/Sample-"${NAME}"/"${NAME}"_CompetitiveMappingReport.txt | grep -v "Total Q30 mapped reads"`
TotalNoDuplicates=`grep -A 1 "Total after removing PCR duplicates" $PATHS/Sample-"${NAME}"/"${NAME}"_CompetitiveMappingReport.txt | grep -v "Total after removing PCR duplicates"`
TotalUnique=`grep -A 1 "Total unique reads" $PATHS/Sample-"${NAME}"/"${NAME}"_CompetitiveMappingReport.txt | grep -v "Total unique reads"`

# Prepare variables for species
Scristatus=`grep "S.cristatus" $PATHS/Sample-"${NAME}"/"${NAME}"_CompetitiveMappingReport.txt | cut -f2`
Sgordonii=`grep "S.gordonii" $PATHS/Sample-"${NAME}"/"${NAME}"_CompetitiveMappingReport.txt | cut -f2`
Smitis=`grep "S.mitis" $PATHS/Sample-"${NAME}"/"${NAME}"_CompetitiveMappingReport.txt | cut -f2`
Soralis=`grep "S.oralis" $PATHS/Sample-"${NAME}"/"${NAME}"_CompetitiveMappingReport.txt | cut -f2`
Ssanguinis=`grep "S.sanguinis" $PATHS/Sample-"${NAME}"/"${NAME}"_CompetitiveMappingReport.txt | cut -f2`
SspDD04=`grep "S.spDD04" $PATHS/Sample-"${NAME}"/"${NAME}"_CompetitiveMappingReport.txt | cut -f2`

# Collect non-unique vs unique reads data
ScristatusNU=`echo $Scristatus | cut -d " " -f2`
SgordoniiNU=`echo $Sgordonii | cut -d " " -f2`
SmitisNU=`echo $Smitis | cut -d " " -f2`
SoralisNU=`echo $Soralis | cut -d " " -f2`
SsanguinisNU=`echo $Ssanguinis | cut -d " " -f2`
SspDD04NU=`echo $SspDD04 | cut -d " " -f2`

ScristatusU=`echo $Scristatus | cut -d " " -f3`
SgordoniiU=`echo $Sgordonii | cut -d " " -f3`
SmitisU=`echo $Smitis | cut -d " " -f3`
SoralisU=`echo $Soralis | cut -d " " -f3`
SsanguinisU=`echo $Ssanguinis | cut -d " " -f3`
SspDD04U=`echo $SspDD04 | cut -d " " -f3`

# Print to file
 echo -e "${NAME} \t $AnalysisReady \t $TotalMapped \t $TotalQ30 \t $TotalNoDuplicates \t $TotalUnique \t $ScristatusNU \t $SgordoniiNU \t $SmitisNU \t $SoralisNU \t $SsanguinisNU \t $SspDD04U \t $ScristatusU \t $SgordoniiU \t $SmitisU \t $SoralisU \t $SsanguinisU \t $SspDD04U \t" >> $PATHS/CollectedCompetitiveMappingResults.txt
 
# Remove intermediate files and make directory for SAM and BAM files
rm $PATHS/Sample-"${NAME}"/"${NAME}"_firstalignment.sai
rm $PATHS/Sample-"${NAME}"/"${NAME}"_mapped_sorted.bai
rm $PATHS/Sample-"${NAME}"/"${NAME}"_mapped_q30_sorted_rmdup.bai
rm $PATHS/Sample-"${NAME}"/"${NAME}"_mapped_q30_sorted_rmdup_uniq.bai
mkdir $PATHS/Sample-"${NAME}"/"${NAME}"_SAMandBAMfiles
mv $PATHS/Sample-"${NAME}"/"${NAME}"_*.bam $PATHS/Sample-"${NAME}"/"${NAME}"_SAMandBAMfiles
mv $PATHS/Sample-"${NAME}"/"${NAME}"_*.sam $PATHS/Sample-"${NAME}"/"${NAME}"_SAMandBAMfiles

done
  


