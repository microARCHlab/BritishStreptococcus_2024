## Usage: Determine presence/absence of previously-identified genomic islands based on overall coverage of that region. Also create attribute tables for the islands and samples (including relative abundance in competitive mapping) for use in heatmap plotting. This script goes through the three S. sanguinis sample sets (ancient, healthy modern, and modern with periodontal disease). Obtain coordinates for plotting in circos as highlights for islands >2000 bp. ##

PATHS=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/GenomeAnalysis/GenomicIslands

# Directory for outputs to be made
outdir=$PATHS/ByWholeIsland_SanguinisIndividualMapping-03062024

# Genomic islands input
genomicislands=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/GenomeAnalysis/GenomicIslands/SanguinisGenomicIslands–03062024.txt

# Output of coverage analysis script. Could also do the tsv file with gene info
coveragetable=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/GenomeAnalysis/CoverageAnalysis/Streptococcus_sanguinis/sanguinis-AllSequencesGene-AllGenes_info.txt

# Reference genome for calculating GC content
REF=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/IndividualAlignments/References/Streptococcus_sanguinis/Streptococcus_sanguinis_SK405

# Breadth of coverage data
ancientmappingsummary=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/IndividualAlignments/ResultsSummaries/Streptococcus_sanguinis_MappingResults.txt

modernmappingsummary=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Streptococcus_sanguinis_HealthyModernMappingResults.txt

periomodernmappingsummary=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/IndividualAlignments/Streptococcus_sanguinis_PerioModernMappingResults.txt

# Competitive mapping data (for calculating relative abundance)
ancientcompetitivemappingsummary=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples/AllSpeciesCollectedCompetitiveMappingResults.txt

moderncompetitivemappingsummary=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/HealthyModernSamples/HealthyModernAllSpeciesCollectedCompetitiveMappingResults.txt

periomoderncompetitivemappingsummary=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/PerioModernSamples/PerioModernAllSpeciesCollectedCompetitiveMappingResults.txt

# Species name in competitive mapping table
speciesname=S.sanguinis

# Formatted dates (for ancient samples)
metadata=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/GenomeAnalysis/BritishFormattedDates-02292024.txt

# Specify minimum overall breadth of coverage for a sample to run the analysis on. As a whole number/percent
minoverallcoverage=50

# For naming/finding the bam files within the sample folders
COVERAGEFOR=Streptococcus_sanguinis

module use /storage/group/LiberalArts/default/lsw132_collab/sw8/modules
module load bedtools/2.31.1
##########################################################################################
# Find which ancient samples have enough coverage to analyze
# While I filtered here, I ultimately filtered further in R when making the plots, so this isn't all that necessary.

mkdir $outdir
cd $outdir

for Sample in $PATHS/../../IndividualAlignments/Samples/Sample-*
 do

  ID=$(basename $Sample)
  NAME=$(echo $ID |cut -d "-" -f2)
  
if [[ "$NAME" == *"EBC"* ]]
  then
  echo "$NAME is an EBC"
 else
 
  breadthcov=`grep "$NAME" $ancientmappingsummary | cut -f8`
  breadthcovpercent=`echo "$breadthcov * 1000000" | bc`
  # Multiplying is necessary to be able to work with (and remove) really low coverage values when shell does not work with decimals.
  breadthcovpercentremovedecimal=`echo $breadthcovpercent | cut -d "." -f1`
  breadthcovpercentremoveddecimaladdone=`echo $breadthcovpercentremovedecimal + 1 | bc`

  mincovaddone=`echo "$minoverallcoverage + 1" | bc`
  mincovmultiplied=`echo "$mincovaddone * 10000" | bc`

if [ "$breadthcovpercentremoveddecimaladdone" -lt "$mincovmultiplied" ]
 then
  echo "Not enough breadth for $NAME"
 else
 echo "Enough breadth for $NAME"
 echo "$NAME" >> $outdir/"$COVERAGEFOR"_AnalyzableSampleList_Ancient.txt
fi
fi

done

############################################
# Find which healthy modern samples have enough coverage to analyze

for Sample in /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/HealthyModernSamples/Sample-*
 do

  ID=$(basename $Sample)
  NAME=$(echo $ID |cut -d "-" -f2)
  
if [[ "$NAME" == *"EBC"* ]]
  then
  echo "$NAME is an EBC"
 else
 
  breadthcov=`grep "$NAME" $modernmappingsummary | cut -f8`
  breadthcovpercent=`echo "$breadthcov * 1000000" | bc`
  # Multiplying is necessary to be able to work with (and remove) really low coverage values when shell does not work with decimals.
  breadthcovpercentremovedecimal=`echo $breadthcovpercent | cut -d "." -f1`
  breadthcovpercentremoveddecimaladdone=`echo $breadthcovpercentremovedecimal + 1 | bc`

  mincovaddone=`echo "$minoverallcoverage + 1" | bc`
  mincovmultiplied=`echo "$mincovaddone * 10000" | bc`

if [ "$breadthcovpercentremoveddecimaladdone" -lt "$mincovmultiplied" ]
 then
  echo "Not enough breadth for $NAME"
 else
 echo "Enough breadth for $NAME"
 echo "$NAME" >> $outdir/"$COVERAGEFOR"_AnalyzableSampleList_HealthyModern.txt
fi
fi

done

############################################
# Find which perio modern samples have enough coverage to analyze

for Sample in /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/PerioModernSamples/Sample-*
 do

  ID=$(basename $Sample)
  NAME=$(echo $ID |cut -d "-" -f2)
  
if [[ "$NAME" == *"EBC"* ]]
  then
  echo "$NAME is an EBC"
 else
 
  breadthcov=`grep "$NAME" $periomodernmappingsummary | cut -f8`
  breadthcovpercent=`echo "$breadthcov * 1000000" | bc`
  # Multiplying is necessary to be able to work with (and remove) really low coverage values when shell does not work with decimals.
  breadthcovpercentremovedecimal=`echo $breadthcovpercent | cut -d "." -f1`
  breadthcovpercentremoveddecimaladdone=`echo $breadthcovpercentremovedecimal + 1 | bc`

  mincovaddone=`echo "$minoverallcoverage + 1" | bc`
  mincovmultiplied=`echo "$mincovaddone * 10000" | bc`

if [ "$breadthcovpercentremoveddecimaladdone" -lt "$mincovmultiplied" ]
 then
  echo "Not enough breadth for $NAME"
 else
 echo "Enough breadth for $NAME"
 echo "$NAME" >> $outdir/"$COVERAGEFOR"_AnalyzableSampleList_PerioModern.txt
fi
fi

done

##########################################################################################
# Get island information

echo -e "Genomic Island" > $outdir/"$COVERAGEFOR"_GenomicIslandCoverageSummary_Ancient.txt
echo -e "Genomic Island" > $outdir/"$COVERAGEFOR"_GenomicIslandCoverageSummary_HealthyModern.txt
echo -e "Genomic Island" > $outdir/"$COVERAGEFOR"_GenomicIslandCoverageSummary_PerioModern.txt

echo -e "Genomic Island \t Island Length \t Number of Genes" > $outdir/"$COVERAGEFOR"_GenomicIslandAttributes.txt

# Get a non-redundant list of island names
# Output a single copy of each repeated line
cat $genomicislands | cut -f1 | uniq -d > islandlist.txt
# Output lines that were not repeated as well, without including header
cat $genomicislands | cut -f1 | uniq -u >> islandlist.txt

# Get all genes part of that genomic island into the same file
mkdir GenomicIslands
while read IslandName
 do 
 	touch $outdir/GenomicIslands/"$IslandName".txt 
 	grep "$IslandName" $genomicislands >> $outdir/GenomicIslands/"$IslandName".txt 
done < islandlist.txt

rm islandlist.txt

for GIfile in $outdir/GenomicIslands/GI*.txt
 do
 	touch $outdir/tempcoordinates.txt
	touch $outdir/tempcoordinatessorted.txt
 	cat $GIfile | cut -f4 > tempgenefile.txt
	ID=$(basename $GIfile)
	NAME=`echo $ID | cut -d "." -f1`
	
# Find the largest and smallest coordinates of genes in the genomic islands (in case they are not list in order or reversed). Use these as the starting and ending positions of the genomic island.
# ALL GENES HAVE TO BE ON SAME CONTIG FOR THIS TO WORK PROPERLY
	while read genomicislandline
	 do
	 locus=`echo "$genomicislandline" | cut -f4`
	 genestart=`grep "$locus" $coveragetable | cut -f2`
	 geneend=`grep "$locus" $coveragetable | cut -f3`
	 echo "Finding coordinates for $locus"
	 
	 echo "$genestart" >> $outdir/tempcoordinates.txt
	 echo "$geneend" >> $outdir/tempcoordinates.txt
	done < $GIfile
	
  cat $outdir/tempcoordinates.txt | sort -n > $outdir/tempcoordinatessorted.txt
  startingpos=`head -n 1 $outdir/tempcoordinatessorted.txt`
  endingpos=`tail -n 1 $outdir/tempcoordinatessorted.txt`
  anylocus=`head -n 1 $GIfile | cut -f4`
  contig=`grep "$anylocus" $coveragetable | cut -f1`

# Add to combined bed file
echo "Adding coordinates for $NAME to combined bed file"
# Make sure for a single gene the coordinates are still in order
if [ $startingpos -lt $endingpos ]
then
echo -e "$contig\t$startingpos\t$endingpos" >> $outdir/combinedgenomicislandcoordinates.bed
else
echo -e "$contig\t$endingpos\t$startingpos" >> $outdir/combinedgenomicislandcoordinates.bed
fi

# Total length (bp) as absolute value
totallength=$(echo "$endingpos - $startingpos" | bc)
totallength="${totallength/#-}"

# Number of genes (any annotated, protein-coding or not)
numberofgenes=`cat $GIfile | wc -l`

# Add list to what will be the first column of the summary file
echo -e "$NAME" >> $outdir/"$COVERAGEFOR"_GenomicIslandCoverageSummary_Ancient.txt
echo -e "$NAME" >> $outdir/"$COVERAGEFOR"_GenomicIslandCoverageSummary_HealthyModern.txt
echo -e "$NAME" >> $outdir/"$COVERAGEFOR"_GenomicIslandCoverageSummary_PerioModern.txt

# Also add to a second attributes table
echo -e "$NAME \t $totallength \t $numberofgenes" >> $outdir/"$COVERAGEFOR"_GenomicIslandAttributes.txt

rm $outdir/tempcoordinates.txt
rm $outdir/tempcoordinatessorted.txt

done

# Calculate GC content
echo "GC content" > $outdir/"$COVERAGEFOR"_GCcontent.txt
bedtools nuc -fi $REF -bed $outdir/combinedgenomicislandcoordinates.bed | grep -v "#" | cut -f5 >> $outdir/"$COVERAGEFOR"_GCcontent.txt

# Add GC content to the attribute file
cat $outdir/"$COVERAGEFOR"_GenomicIslandAttributes.txt > $outdir/tempcombined.txt
paste $outdir/tempcombined.txt $outdir/"$COVERAGEFOR"_GCcontent.txt > $outdir/"$COVERAGEFOR"_GenomicIslandAttributes.txt

# Make a separate attributes table for ancient and modern. Keep same formatting to make it easy to combine.
echo -e "Sample \t Whole Genome Coverage \t Relative Abundance \t Date Range \t Early Date \t Late Date \t Date" > $outdir/"$COVERAGEFOR"_SampleAttributes_Ancient.txt

echo -e "Sample \t Whole Genome Coverage \t Relative Abundance \t Date Range \t Early Date \t Late Date \t Date" > $outdir/"$COVERAGEFOR"_SampleAttributes_HealthyModern.txt

echo -e "Sample \t Whole Genome Coverage \t Relative Abundance \t Date Range \t Early Date \t Late Date \t Date" > $outdir/"$COVERAGEFOR"_SampleAttributes_PerioModern.txt

##########################################################################################
# Test coverage of these genomic islands in each ancient sample
head -n 1 $ancientcompetitivemappingsummary | tr '\t' '\n' > $outdir/headersverticallist.txt
while read NAME
 do
  echo "Running genomic island breadth analysis on "${NAME}"..."

# Get associated data  
breadthcov=`grep "$NAME" $ancientmappingsummary | cut -f8`
Date=`grep "$NAME" $metadata | cut -f 2`
EarlyDate=`grep "$NAME" $metadata | cut -f 4`
LateDate=`grep "$NAME" $metadata | cut -f 6`
BlockPeriod=`grep "$NAME" $metadata | cut -f 9`

totalstreptococcimappedreadsCM=`grep "$NAME" $ancientcompetitivemappingsummary | cut -f 6`
col=`grep -n -m 1 "$speciesname" $outdir/headersverticallist.txt | cut -d ":" -f1`
mappedtospecies=`grep "$NAME" $ancientcompetitivemappingsummary | cut -f$col`
relativeabundance=`echo "scale=3; $mappedtospecies / $totalstreptococcimappedreadsCM" | bc`

    
# Get breadth of coverage data for the islands (single column with name on top)
echo "$NAME" > $outdir/tempbreadth.txt
bedtools coverage -a $outdir/combinedgenomicislandcoordinates.bed -b $PATHS/../../IndividualAlignments/Samples/Sample-"$NAME"/$COVERAGEFOR/"$NAME"_mapped_q30_sorted_rmdup_uniq.bam | cut -f7 >> $outdir/tempbreadth.txt

# Paste it into the summary file. Need to use temporary files to be able to add to themselves.
cat $outdir/"$COVERAGEFOR"_GenomicIslandCoverageSummary_Ancient.txt > $outdir/tempcombined.txt
paste $outdir/tempcombined.txt $outdir/tempbreadth.txt > $outdir/"$COVERAGEFOR"_GenomicIslandCoverageSummary_Ancient.txt

# Also add to separate attribute table
echo -e "$NAME \t $breadthcov \t $relativeabundance \t $Date \t $EarlyDate \t $LateDate \t $BlockPeriod" >> $outdir/"$COVERAGEFOR"_SampleAttributes_Ancient.txt

done < $outdir/"$COVERAGEFOR"_AnalyzableSampleList_Ancient.txt

rm $outdir/headersverticallist.txt
rm $outdir/tempbreadth.txt
rm $outdir/tempcombined.txt
rm $outdir/tempgenefile.txt

###################################
# Test coverage of these genomic islands in each healthy modern sample

head -n 1 $moderncompetitivemappingsummary | tr '\t' '\n' > $outdir/headersverticallist.txt
while read NAME
 do
  echo "Running genomic island breadth analysis on "${NAME}"..."

# Get associated data  
breadthcov=`grep "$NAME" $modernmappingsummary | cut -f8`

totalstreptococcimappedreadsCM=`grep "$NAME" $moderncompetitivemappingsummary | cut -f 6`
col=`grep -n -m 1 "$speciesname" $outdir/headersverticallist.txt | cut -d ":" -f1`
mappedtospecies=`grep "$NAME" $moderncompetitivemappingsummary | cut -f$col`
relativeabundance=`echo "scale=3; $mappedtospecies / $totalstreptococcimappedreadsCM" | bc`

# Get breadth of coverage data for the islands (single column with name on top)
echo "$NAME" > $outdir/tempbreadth.txt
bedtools coverage -a $outdir/combinedgenomicislandcoordinates.bed -b /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/HealthyModernSamples/Sample-"$NAME"/$COVERAGEFOR/"$NAME"_mapped_q30_sorted_rmdup_uniq.bam | cut -f7 >> $outdir/tempbreadth.txt

# Paste it into the summary file. Need to use temporary files to be able to add to themselves.
cat $outdir/"$COVERAGEFOR"_GenomicIslandCoverageSummary_HealthyModern.txt > $outdir/tempcombined.txt
paste $outdir/tempcombined.txt $outdir/tempbreadth.txt > $outdir/"$COVERAGEFOR"_GenomicIslandCoverageSummary_HealthyModern.txt

# Also add to separate attribute table
# Put in arbitrary "Early Date" for ordering purposes
echo -e "$NAME \t $breadthcov \t $relativeabundance \t Healthy Modern \t 2016 \t Healthy Modern \t Healthy Modern" >> $outdir/"$COVERAGEFOR"_SampleAttributes_HealthyModern.txt

done < $outdir/"$COVERAGEFOR"_AnalyzableSampleList_HealthyModern.txt

rm $outdir/headersverticallist.txt
rm $outdir/tempbreadth.txt
rm $outdir/tempcombined.txt

###################################
# Test coverage of these genomic islands in each perio modern sample

head -n 1 $periomoderncompetitivemappingsummary | tr '\t' '\n' > $outdir/headersverticallist.txt
while read NAME
 do
  echo "Running genomic island breadth analysis on "${NAME}"..."

# Get associated data  
breadthcov=`grep "$NAME" $periomodernmappingsummary | cut -f8`

totalstreptococcimappedreadsCM=`grep "$NAME" $periomoderncompetitivemappingsummary | cut -f 6`
col=`grep -n -m 1 "$speciesname" $outdir/headersverticallist.txt | cut -d ":" -f1`
mappedtospecies=`grep "$NAME" $periomoderncompetitivemappingsummary | cut -f$col`
relativeabundance=`echo "scale=3; $mappedtospecies / $totalstreptococcimappedreadsCM" | bc`

# Get breadth of coverage data for the islands (single column with name on top)
echo "$NAME" > $outdir/tempbreadth.txt
bedtools coverage -a $outdir/combinedgenomicislandcoordinates.bed -b /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/PerioModernSamples/Sample-"$NAME"/$COVERAGEFOR/"$NAME"_mapped_q30_sorted_rmdup_uniq.bam | cut -f7 >> $outdir/tempbreadth.txt

# Paste it into the summary file. Need to use temporary files to be able to add to themselves
cat $outdir/"$COVERAGEFOR"_GenomicIslandCoverageSummary_PerioModern.txt > $outdir/tempcombined.txt
paste $outdir/tempcombined.txt $outdir/tempbreadth.txt > $outdir/"$COVERAGEFOR"_GenomicIslandCoverageSummary_PerioModern.txt

# Also add to separate attribute table
# Put in arbitrary "Early Date" for ordering purposes
echo -e "$NAME \t $breadthcov \t $relativeabundance \t Perio Modern \t 2017 \t Perio Modern \t Perio Modern" >> $outdir/"$COVERAGEFOR"_SampleAttributes_PerioModern.txt

done < $outdir/"$COVERAGEFOR"_AnalyzableSampleList_PerioModern.txt

rm $outdir/headersverticallist.txt
rm $outdir/tempbreadth.txt
rm $outdir/tempcombined.txt

##########################################################################################
# Make combined attribute and summary tables

cat $outdir/"$COVERAGEFOR"_SampleAttributes_Ancient.txt > $outdir/"$COVERAGEFOR"_SampleAttributes_Combined.txt
cat $outdir/"$COVERAGEFOR"_SampleAttributes_HealthyModern.txt | grep -v "Sample" >> $outdir/"$COVERAGEFOR"_SampleAttributes_Combined.txt
cat $outdir/"$COVERAGEFOR"_SampleAttributes_PerioModern.txt | grep -v "Sample" >> $outdir/"$COVERAGEFOR"_SampleAttributes_Combined.txt

cat $outdir/"$COVERAGEFOR"_GenomicIslandCoverageSummary_HealthyModern.txt | cut -f2- > $outdir/"$COVERAGEFOR"_GenomicIslandCoverageSummary_HealthyModerntemp.txt
cat $outdir/"$COVERAGEFOR"_GenomicIslandCoverageSummary_PerioModern.txt  | cut -f2- > $outdir/"$COVERAGEFOR"_GenomicIslandCoverageSummary_PerioModerntemp.txt 

paste $outdir/"$COVERAGEFOR"_GenomicIslandCoverageSummary_Ancient.txt $outdir/"$COVERAGEFOR"_GenomicIslandCoverageSummary_HealthyModerntemp.txt > tempcombined.txt
paste tempcombined.txt $outdir/"$COVERAGEFOR"_GenomicIslandCoverageSummary_PerioModerntemp.txt > $outdir/"$COVERAGEFOR"_GenomicIslandCoverageSummary_Combined.txt

rm $outdir/"$COVERAGEFOR"_GenomicIslandCoverageSummary_HealthyModerntemp.txt
rm $outdir/"$COVERAGEFOR"_GenomicIslandCoverageSummary_PerioModerntemp.txt 
rm tempcombined.txt

##########################################################################################
# Get coordinates for plotting in circos

name=`basename $outdir | cut -d "_" -f2`
touch $outdir/"$name"_islandcircoscoordinates.txt

for GIfile in $outdir/GenomicIslands/GI*.txt
 do
	ID=$(basename $GIfile)
	NAME=`echo $ID | cut -d "." -f1`
    length=`grep "$NAME" "$COVERAGEFOR"_GenomicIslandAttributes.txt | cut -f2 | tr -d " "`
    if [[ $length -gt 2000 ]]
    then
	echo "Finding circos coordinates for $NAME"
	while read genomicislandline
	 do
	 locus=`echo "$genomicislandline" | cut -f4`
	 genestart=`grep "$locus" $coveragetable | cut -f2`
	 geneend=`grep "$locus" $coveragetable | cut -f3`
	 
	 echo "$genestart" >> $outdir/tempcoordinates.txt
	 echo "$geneend" >> $outdir/tempcoordinates.txt
	done < $GIfile
	
  cat $outdir/tempcoordinates.txt | sort -n > $outdir/tempcoordinatessorted.txt
  startingpos=`head -n 1 $outdir/tempcoordinatessorted.txt`
  endingpos=`tail -n 1 $outdir/tempcoordinatessorted.txt`
  anylocus=`grep "$NAME" $GIfile | head -n 1 | cut -f4`
  contig=`grep "$anylocus" $coveragetable | cut -f1`

echo -e "$contig \t $startingpos \t $endingpos" >> $outdir/"$name"_islandcircoscoordinates.txt

rm $outdir/tempcoordinates.txt
rm $outdir/tempcoordinatessorted.txt

else 
:
fi
done

