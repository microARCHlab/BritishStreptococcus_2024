## OCT 2023 ##
## Ava Gabrys ##

### Maps to a joined reference and produces a summary table of unfiltered and filtered reads mapped to each reference (totaling reads from all scaffolds/contigs) ###

## Prepare the concatenated reference ##
# Change the reference name in each FASTA file (including for all contigs/scaffolds/chromosomes) to include species name before it (ex: >NZ_LR134336.1 becomes >S.oralisNZ_LR134336.1), put FASTAs to be joined in the same directory, and use: cat *.fna > combinedreference
# Index the reference with samtools and BWA to prepare for mapping
# Create a list of all of these contig/scaffold/chromosome labels (1 for each species/original reference FASTA included in the combined file)

## List of identifying names used in combined FASTA reference ##
list=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/IdentiferList.txt

## Define paths ##
# Path to indexed reference
ref=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/JoinedReferences/NewReference/NewReference.fna

# Path to directory with Samples (each in their own directory)
PATHS=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples

# Out directory name (to be created within sample directory)
OUTDIR=AllSpecies

# Load modules
module load bwa samtools

##############################################################################
# Prepare the summary table
# Get headers based on species names from provided list
touch $PATHS/existingheaders.txt

headers=`cat $PATHS/existingheaders.txt`
while read species
 do
 	headers=`cat $PATHS/existingheaders.txt`
 	echo -e "$headers \t $species" > $PATHS/existingheaders.txt
done < $list

finalheaders=`cat $PATHS/existingheaders.txt`
 	
echo -e "Sample \t Analysis-ready reads \t Total mapped reads \t Total Q30 mapped reads \t  Total after removing PCR duplicates \t Total unique reads $finalheaders \t" > $PATHS/"$OUTDIR"CollectedCompetitiveMappingResults.txt

# Start of loop

for Sample in $PATHS/Sample*;
 do
  cd $Sample

  ID=$(basename $Sample)
  NAME=$(echo $ID |cut -d "-" -f2)
  echo "Running competitve mapping on "${NAME}"..."
  
  mkdir $OUTDIR
  
###########################################

# Alignment with bwa aln
  bwa aln -l 1024 -n 0.01 -o 2 -t 8 $ref $PATHS/Sample-"${NAME}"/*fastq.gz > $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_firstalignment.sai

# Alignment in sam file format
  bwa samse $ref $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_firstalignment.sai $PATHS/Sample-"${NAME}"/"${NAME}"*fastq.gz > $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_firstalignment.sam
  
###########################################
# Filtering, sorting, removing duplicates and non-unique reads using SAMtools

# SAM to BAM, keep header
  samtools view -bSh $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_firstalignment.sam > $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_firstalignment.bam

# Filter out unmappped and low quality reads
  samtools view -bh -F4 $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_firstalignment.bam > $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_mapped.bam
  
# Sort mapped alignments by leftmost coordinates
  samtools sort -o $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_mapped_sorted.bam $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_mapped.bam

  samtools view -bh -q 30  $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_mapped_sorted.bam > $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_mapped_q30_sorted.bam

# Remove PCR duplicate
  samtools rmdup -s $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_mapped_q30_sorted.bam $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_mapped_q30_sorted_rmdup.bam
  
# Remove reads with multiple mappings
samtools view -h $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_mapped_q30_sorted_rmdup.bam |grep -v 'XT:A:R'|grep -v 'XA:Z' |grep -v 'XT:A:M' |awk '{if($0~/X1:i:0/||$0~/^@/)print $0}' | samtools view -bS - > $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_mapped_q30_sorted_rmdup_uniq.bam

###########################################
# Print general alignment statistics

# Create a report summary file for the sample
  echo " " > $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_CompetitiveMappingReport.txt
  echo "${NAME}" >> $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_CompetitiveMappingReport.txt
  echo "---------------------------------------" >> $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_CompetitiveMappingReport.txt

# Analysis-ready reads
 echo "Analysis-ready reads" >> $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_CompetitiveMappingReport.txt
  samtools view -c $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_firstalignment.bam >> $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_CompetitiveMappingReport.txt
# Total mapped reads
  echo "Total mapped reads" >> $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_CompetitiveMappingReport.txt
  samtools view -c $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_mapped.bam >> $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_CompetitiveMappingReport.txt
# Total Q30 mapped reads
  echo "Total Q30 mapped reads" >> $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_CompetitiveMappingReport.txt
  samtools view -c $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_mapped_q30_sorted.bam >> $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_CompetitiveMappingReport.txt
# Total after removing duplicates
  echo "Total after removing PCR duplicates" >> $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_CompetitiveMappingReport.txt
  samtools view -c $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_mapped_q30_sorted_rmdup.bam >> $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_CompetitiveMappingReport.txt
# Total after removing non-unique reads
echo "Total unique reads" >> $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_CompetitiveMappingReport.txt
  samtools view -c $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_mapped_q30_sorted_rmdup_uniq.bam >> $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_CompetitiveMappingReport.txt
  
###########################################
# Calculate and print statistics by species

# Index SAM files to use in SAMtools idxstats
samtools index $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_mapped_sorted.bam $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_mapped_sorted.bai
samtools index $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_mapped_q30_sorted_rmdup_uniq.bam $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_mapped_q30_sorted_rmdup_uniq.bai

# Use SAMtools index to print results based on reference. Columns are reference sequence name, sequence length, number of mapped reads, and number of unmapped reads (which would be 0 because unmapped reads were already filtered out). The -f 1,3 prints only the first and third columns.

########################
# Total reads mapped to species
echo "---------------------------------------" >> $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_CompetitiveMappingReport.txt
echo "Total reads mapped to species" >> $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_CompetitiveMappingReport.txt
samtools idxstats $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_mapped_sorted.bam | cut -f 1,3 >> $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_Totalidxstats.txt

# Add up reads from all contigs for counts by species
while read species
 do
	echo -e "$species all contigs \t `grep "$species" $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_Totalidxstats.txt | cut -f2 | paste -sd+ | bc`" >> $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_Totalidxstats.txt
 done < $list

# Print these total counts to samples summary file
grep "all contigs" $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_Totalidxstats.txt >> $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_CompetitiveMappingReport.txt

#########################
# Unique reads mapped to species (Q30+, PCR duplicates removed)
echo "---------------------------------------" >> $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_CompetitiveMappingReport.txt
echo "Unique reads mapped to species (Q30+, PCR duplicates removed)" >> $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_CompetitiveMappingReport.txt
samtools idxstats $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_mapped_q30_sorted_rmdup_uniq.bam | cut -f 1,3 >> $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_Uidxstats.txt

# Add up reads from all contigs
while read species
 do
echo -e "$species all contigs \t `grep "$species" $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_Uidxstats.txt | cut -f2 | paste -sd+ | bc`" >> $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_Uidxstats.txt
 done < $list

# Print these total counts to sample summary file
grep "all contigs" $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_Uidxstats.txt >> $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_CompetitiveMappingReport.txt

###########################################
# Add to tabbed results file for all samples

# Collect general data
AnalysisReady=`grep -A 1 "Analysis-ready reads" $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_CompetitiveMappingReport.txt | grep -v "Analysis-ready reads"`
TotalMapped=`grep -A 1 "Total mapped reads" $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_CompetitiveMappingReport.txt | grep -v "Total mapped reads"`
TotalQ30=`grep -A 1 "Total Q30 mapped reads" $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_CompetitiveMappingReport.txt | grep -v "Total Q30 mapped reads"`
TotalNoDuplicates=`grep -A 1 "Total after removing PCR duplicates" $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_CompetitiveMappingReport.txt | grep -v "Total after removing PCR duplicates"`
TotalUnique=`grep -A 1 "Total unique reads" $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_CompetitiveMappingReport.txt | grep -v "Total unique reads"`

# Prepare to paste read data by species
touch $PATHS/existingspeciescounts.txt
while read species
 do 
 	UniqueReads=`grep "$species" 	$PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_CompetitiveMappingReport.txt | cut -f2 | tail -1`
 existingspeciescounts=`cat $PATHS/existingspeciescounts.txt`
	# Add onto the tabbed file
	echo -e "$existingspeciescounts \t $UniqueReads" > $PATHS/existingspeciescounts.txt
done < $list

finalspeciescounts=`cat $PATHS/existingspeciescounts.txt`

# Paste all to summary file
echo -e "${NAME} \t $AnalysisReady \t $TotalMapped \t $TotalQ30 \t $TotalNoDuplicates \t $TotalUnique $finalspeciescounts" >> $PATHS/"$OUTDIR"CollectedCompetitiveMappingResults.txt

###########################################
# Remove intermediate files and make directory for SAM and BAM files
rm $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_firstalignment.sai
rm $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_mapped_sorted.bai
rm $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_mapped_q30_sorted_rmdup.bai
rm $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_mapped_q30_sorted_rmdup_uniq.bai
mkdir $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_SAMandBAMfiles
mv $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_*.bam $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_SAMandBAMfiles
mv $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_*.sam $PATHS/Sample-"${NAME}"/$OUTDIR/"${NAME}"_SAMandBAMfiles
rm $PATHS/existingheaders.txt
rm $PATHS/existingspeciescounts.txt

done