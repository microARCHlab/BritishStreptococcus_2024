## Usage: Calculate genome coverage statistics (breadth and depth) for references mapped to in competitive mapping. Isolates reads mapped to a particular specified reference and extracts coverage information from that BAM file.

# PATHS
PATHS=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples

# Output directory for consensus sequences
outdir=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/Phylogenetics/ConsensusSequences/AllSpeciesCompetitiveMapping

# Common identifier in all contigs/chromosomes for species of interest, as was used in competitive mapping
identifier=S.spDD04

# Species name (for file naming purposes)
species=DD04

# Load modules
module load samtools/1.18 

#######################################################################################
## Getting the reads only mapped to the species of interest ##

# Use faidx file of combined reference to get names of chromosome/contigs/scaffolds of interest. These were named with appropriate species identifier
grep "$identifier" $FAIDX | cut -f1 > $PATHS/"species"_regions.txt

# Prepare summary file
echo -e "Sample \t Total mapped filtered reads \t Average read length \t Breadth of coverage \t Depth of coverage" > $PATHS/"$species"_competitivemappingcoveragesummary.txt

for Sample in $PATHS/Sample-*
 do
  cd $Sample

  ID=$(basename $Sample)
  NAME=$(echo $ID |cut -d "-" -f2)

echo "Finding stats for $NAME"
 
# Extract these regions from the sorted filtered bam file for selected samples.
BAM=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples/Sample-"$NAME"/AllSpecies/"$NAME"_SAMandBAMfiles/"$NAME"_mapped_q30_sorted_rmdup_uniq.bam

samtools index $BAM # BAM files need to be indexed
cat $PATHS/"species"_regions.txt | tr "\n" " " | xargs samtools view -bh $BAM > /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples/Sample-"$NAME"/AllSpecies/"$NAME"_SAMandBAMfiles/"$NAME"_filtered_"$species".bam

NEWBAM=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples/Sample-"$NAME"/AllSpecies/"$NAME"_SAMandBAMfiles/"$NAME"_filtered_"$species".bam

#######################################################################################
## Gather and print alignment statistics for this extracted BAM ##

# Create a report summary file
  echo " " > /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples/Sample-"$NAME"/AllSpecies/"$NAME"_"$species"_CMCoverageReport.txt
  echo "${NAME}" >> /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples/Sample-"$NAME"/AllSpecies/"$NAME"_"$species"_CMCoverageReport.txt
  echo "---------------------------------------" >> /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples/Sample-"$NAME"/AllSpecies/"$NAME"_"$species"_CMCoverageReport.txt

# Print general data
# Total after removing non-unique reads
echo "Total unique reads" >> /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples/Sample-"$NAME"/AllSpecies/"$NAME"_"$species"_CMCoverageReport.txt
  samtools view -c $NEWBAM >> /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples/Sample-"$NAME"/AllSpecies/"$NAME"_"$species"_CMCoverageReport.txt
# Average length of reads
  echo "Average length of reads" >> /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples/Sample-"$NAME"/AllSpecies/"$NAME"_"$species"_CMCoverageReport.txt
  samtools view $NEWBAM |awk '{SUM+=length($10);DIV++}END{print SUM/DIV}' >> /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples/Sample-"$NAME"/AllSpecies/"$NAME"_"$species"_CMCoverageReport.txt
  
# Bring to variables for later printing in summary table
TotalUnique=`grep -A 1 "Total unique reads" /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples/Sample-"$NAME"/AllSpecies/"$NAME"_"$species"_CMCoverageReport.txt | grep -v "Total unique reads"`
AvgLength=`grep -A 1 "Average length of reads" /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples/Sample-"$NAME"/AllSpecies/"$NAME"_"$species"_CMCoverageReport.txt | grep -v "Average length of reads"`

# Calculate and print coverage statistics
# Use SAMtools Coverage
samtools coverage $NEWBAM -o /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples/Sample-"$NAME"/AllSpecies/"$NAME"_"$species"_Coverage.txt
# For breadth of coverage
covbases=`grep -v "name" /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples/Sample-"$NAME"/AllSpecies/"$NAME"_"$species"_Coverage.txt | cut -f5 | paste -sd+ | bc `
numberofbases=`grep "$identifier" /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples/Sample-"$NAME"/AllSpecies/"$NAME"_"$species"_Coverage.txt | cut -f3 | paste -sd+ | bc`
breadthofcoverage=`echo $covbases / $numberofbases | bc -l`

# For depth of coverage
totalmappinglength=`echo $AvgLength \* $TotalUnique | bc -l`
depthofcoverage=`echo $totalmappinglength / $numberofbases | bc -l`
# Printing
# Breadth of coverage
  echo "Breadth of coverage" > /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples/Sample-"$NAME"/AllSpecies/"$NAME"_"$species"_CMCoverageReport.txt
  echo "$breadthofcoverage" >> /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples/Sample-"$NAME"/AllSpecies/"$NAME"_"$species"_CMCoverageReport.txt
# Depth of coverage
  echo "Depth of coverage" >> /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples/Sample-"$NAME"/AllSpecies/"$NAME"_"$species"_CMCoverageReport.txt
  echo "$depthofcoverage" >> /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples/Sample-"$NAME"/AllSpecies/"$NAME"_"$species"_CMCoverageReport.txt

# Collect variables
BreadthofCoverage=`grep -A 1 "Breadth of coverage" /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples/Sample-"$NAME"/AllSpecies/"$NAME"_"$species"_CMCoverageReport.txt | grep -v "Breadth of coverage"`
DepthofCoverage=`grep -A 1 "Depth of coverage" /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples/Sample-"$NAME"/AllSpecies/"$NAME"_"$species"_CMCoverageReport.txt | grep -v "Depth of coverage"`

# Print to summary file
echo -e "$NAME \t $TotalUnique \t $AvgLength \t $BreadthofCoverage \t $DepthofCoverage" >> $PATHS/"$species"_competitivemappingcoveragesummary.txt

done

