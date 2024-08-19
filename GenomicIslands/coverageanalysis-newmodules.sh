## Find breadth and depth of coverage for each coding sequence in the annotated genome. Later used to help identify genomic islands/gaps in coverage. ##

# Define paths
PATHS=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/GenomeAnalysis/CoverageAnalysis
# Name of directory within each sample folder that has the relevant mapping data
COVERAGEFOR=Streptococcus_sanguinis
# Specify sample type for directory naming to prevent overwrite when running again for modern, competitive mapping, types of sequences, etc.
TYPE=PerioModernIndividualMapping
# Reference FASTA (to calculate GC content)
REF=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/IndividualAlignments/References/Streptococcus_sanguinis/Streptococcus_sanguinis_SK405

# Needed starting file: *.tsv of CDS (or could really pick any genes of interest) downloaded from annotation details of genome assembly (ex: https://www.ncbi.nlm.nih.gov/datasets/gene/GCF_000767835.1/). 
# Make sure these chromosome names have the additional species identifier if appropriate (if using with references that were used in CM, where names had species added to them)
STARTINGTSV=$PATHS/$COVERAGEFOR/sanguinis-AllSequencesGene-AllGenes.tsv

# load module
module use /storage/group/LiberalArts/default/lsw132_collab/sw8/modules
module load bedtools/2.31.1

#######################################################################################
cd $PATHS/$COVERAGEFOR
# Start the combined summary tables with all sequence information at hand
cat $STARTINGTSV > "$COVERAGEFOR"_BreadthCoverage_"$TYPE"combinedtable.txt
cat $STARTINGTSV > "$COVERAGEFOR"_DepthCoverage_"$TYPE"combinedtable.txt
BASENAME=`basename $STARTINGTSV | cut -d "." -f1`

# Do not need to do the next 3 commands if the same starting information file has already been used.
# Create a bed file with information from the *.tsv with just chromosome, start, and end; get rid of header
cat *tsv | cut -f1,2,3 | grep -v "Accession" > $PATHS/$COVERAGEFOR/"$BASENAME".bed

# Get features of interest for listing in the summary sample summary files (eg. gene name); all features in original *.tsv will be in combined summary files
cat *tsv | cut -f1,2,3,6,11 > $PATHS/$COVERAGEFOR/"$BASENAME"_info.txt

# Calculate GC content with new header
echo "GC content" > $PATHS/$COVERAGEFOR/"$BASENAME"_GCcontent.txt
bedtools nuc -fi $REF -bed $PATHS/$COVERAGEFOR/"$BASENAME".bed | grep -v "#" | cut -f5 >> $PATHS/$COVERAGEFOR/"$BASENAME"_GCcontent.txt

# Add in GC content to the combined tables. Copy combined table as it is and paste from a temporary file to be able to add back to itself with GC content column.
cat "$COVERAGEFOR"_BreadthCoverage_"$TYPE"combinedtable.txt > tempcombined.txt
paste tempcombined.txt $PATHS/$COVERAGEFOR/"$BASENAME"_GCcontent.txt > "$COVERAGEFOR"_BreadthCoverage_"$TYPE"combinedtable.txt

cat "$COVERAGEFOR"_DepthCoverage_"$TYPE"combinedtable.txt > tempcombined.txt
paste tempcombined.txt $PATHS/$COVERAGEFOR/"$BASENAME"_GCcontent.txt > "$COVERAGEFOR"_DepthCoverage_"$TYPE"combinedtable.txt

mkdir $PATHS/$COVERAGEFOR/SampleCoverages_"$TYPE"

#######################################################################################
# Start of loop

for Sample in /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/PerioModernSamples/Sample-*
#$PATHS/../../IndividualAlignments/ModernSamples/Sample-* #Sample in $PATHS/../../IndividualAlignments/Samples/Sample-*
 do
  cd $Sample

  ID=$(basename $Sample)
  NAME=$(echo $ID |cut -d "-" -f2)
  echo "Running breadth of coverage and depth analysis on "${NAME}"..."

####################################
# Get breadth of coverage data (single column with name on top)
echo "$NAME breadth" > $PATHS/$COVERAGEFOR/tempbreadth.txt
bedtools coverage -a $PATHS/$COVERAGEFOR/"$BASENAME".bed -b $Sample/$COVERAGEFOR/"$NAME"_mapped_q30_sorted_rmdup_uniq.bam | cut -f7 >> $PATHS/$COVERAGEFOR/tempbreadth.txt

# Get depth of coverage data (single column with name on top)
echo "$NAME depth" > $PATHS/$COVERAGEFOR/tempdepth.txt
bedtools coverage -mean -a $PATHS/$COVERAGEFOR/"$BASENAME".bed -b $Sample/$COVERAGEFOR/"$NAME"_mapped_q30_sorted_rmdup_uniq.bam | cut -f4 >> $PATHS/$COVERAGEFOR/tempdepth.txt

# Paste coverages and GC content together into sample summary file with selected information
paste $PATHS/$COVERAGEFOR/"$BASENAME"_info.txt $PATHS/$COVERAGEFOR/"$BASENAME"_GCcontent.txt $PATHS/$COVERAGEFOR/tempbreadth.txt $PATHS/$COVERAGEFOR/tempdepth.txt > $PATHS/$COVERAGEFOR/SampleCoverages_"$TYPE"/"$NAME".txt

####################################
# Make breadth of coverage combined summary table
# Copy combined table as it is and paste from a temporary file to be able to add back to itself with breadth data, which already has the sample name on the column
cat $PATHS/$COVERAGEFOR/"$COVERAGEFOR"_BreadthCoverage_"$TYPE"combinedtable.txt > $PATHS/$COVERAGEFOR/tempcombined.txt
paste $PATHS/$COVERAGEFOR/tempcombined.txt $PATHS/$COVERAGEFOR/tempbreadth.txt > $PATHS/$COVERAGEFOR/"$COVERAGEFOR"_BreadthCoverage_"$TYPE"combinedtable.txt

####################################
# Make depth of coverage combined summary table
# Copy combined table as it is and paste from a temporary file to be able to add back to itself with breadth data, which already has the sample name on the column
cat $PATHS/$COVERAGEFOR/"$COVERAGEFOR"_DepthCoverage_"$TYPE"combinedtable.txt > $PATHS/$COVERAGEFOR/tempcombined.txt
paste $PATHS/$COVERAGEFOR/tempcombined.txt $PATHS/$COVERAGEFOR/tempdepth.txt > $PATHS/$COVERAGEFOR/"$COVERAGEFOR"_DepthCoverage_"$TYPE"combinedtable.txt

done

####################################
# Remove temporary files
rm $PATHS/$COVERAGEFOR/tempcombined.txt 
rm $PATHS/$COVERAGEFOR/tempbreadth.txt
rm $PATHS/$COVERAGEFOR/tempdepth.txt
