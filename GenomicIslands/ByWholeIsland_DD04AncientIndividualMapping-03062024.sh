## Usage: Determine presence/absence of previously-identified genomic islands based on overall coverage of that region. Also create attribute tables for the islands and samples for use in heatmap plotting and coordinates for plotting in circos as highlights for islands >2000 bp. ##

# Define paths
PATHS=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/GenomeAnalysis/GenomicIslands

# Directory for outputs to be made
outdir=$PATHS/ByWholeIsland_DD04AncientIndividualMapping-03062024

# Genomic island info
genomicislands=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/GenomeAnalysis/GenomicIslands/DD04GenomicIslands-03052024.txt


# Output of coverage analysis script. Could also do the tsv file with gene info
coveragetable=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/GenomeAnalysis/CoverageAnalysis/Streptococcus_spDD04/DD04-AllSequencesGene-AllGenes_info.txt

# Reference genome for calculating GC content
REF=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/IndividualAlignments/References/Streptococcus_spDD04/Streptococcus_sp._DD04

# Breadth of coverage data
mappingsummary=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/IndividualAlignments/ResultsSummaries/Streptococcus_spDD04_MappingResults.txt

# Competitive mapping data (for calculating relative abundance)
competitivemappingsummary=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples/AllSpeciesCollectedCompetitiveMappingResults.txt

# Species name in competitive mapping table (for calculating relative abundance)
speciesname=S.spDD04

# Formatted dates
metadata=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/GenomeAnalysis/BritishFormattedDates-02292024.txt

# Specify minimum overall breadth of coverage for a sample to run the analysis on, as a whole number/percent
minoverallcoverage=50

# For naming/finding the bam files within the sample folders
COVERAGEFOR=Streptococcus_spDD04

module use /storage/group/LiberalArts/default/lsw132_collab/sw8/modules
module load bedtools/2.31.1

#######################################################################################
# Find which samples have enough coverage to analyze

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
 
  breadthcov=`grep "$NAME" $mappingsummary | cut -f8`
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
 echo "$NAME" >> $outdir/"$COVERAGEFOR"_AnalyzableSampleList.txt
fi
fi

done

#######################################################################################
# Get island information

echo -e "Genomic Island" > $outdir/"$COVERAGEFOR"_GenomicIslandCoverageSummary.txt

echo -e "Genomic Island \t Island Length \t Number of Genes" > $outdir/"$COVERAGEFOR"_GenomicIslandAttributes.txt

# Get a non-redundant list of island names
# Output a single copy of each repeated line. For DD04 need to remove GI-3- bc extends beyond contig lines
cat $genomicislands | cut -f1 | uniq -d | grep -v "GI-3-" > islandlist.txt
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
	
# Find the largest and smallest coordinates of genes in the genomic islands (in case they are not listed in order or reversed). Use these as the starting and ending positions of the genomic island.
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
echo -e "$NAME" >> $outdir/"$COVERAGEFOR"_GenomicIslandCoverageSummary.txt

# Also add to a second attributes table
echo -e "$NAME \t $totallength \t $numberofgenes" >> $outdir/"$COVERAGEFOR"_GenomicIslandAttributes.txt

rm $outdir/tempcoordinates.txt
rm $outdir/tempcoordinatessorted.txt

done

# Calculate GC content
echo "GC content" > $outdir/"$COVERAGEFOR"_GCcontent.txt
bedtools nuc -fi $REF -bed $outdir/combinedgenomicislandcoordinates.bed | grep -v "#" | cut -f5 >> $outdir/"$COVERAGEFOR"_GCcontent.txt

# Add GC content to the attribute file
#cat $outdir/"$COVERAGEFOR"_GenomicIslandCoverageSummary.txt > $outdir/tempcombined.txt
#paste $outdir/tempcombined.txt $outdir/"$COVERAGEFOR"_GCcontent.txt > "$COVERAGEFOR"_GenomicIslandCoverageSummary.txt

cat $outdir/"$COVERAGEFOR"_GenomicIslandAttributes.txt > $outdir/tempcombined.txt
paste $outdir/tempcombined.txt $outdir/"$COVERAGEFOR"_GCcontent.txt > $outdir/"$COVERAGEFOR"_GenomicIslandAttributes.txt

# Make a separate attributes table
echo -e "Sample \t Whole Genome Coverage \t Relative Abundance \t Date Range \t Early Date \t Late Date \t Date" > $outdir/"$COVERAGEFOR"_SampleAttributes.txt

#######################################################################################
# Test coverage of these genomic islands in each sample
head -n 1 $competitivemappingsummary | tr '\t' '\n' > $outdir/headersverticallist.txt
while read NAME
 do
  echo "Running genomic island breadth analysis on "${NAME}"..."

# Get associated data  
breadthcov=`grep "$NAME" $mappingsummary | cut -f8`
Date=`grep "$NAME" $metadata | cut -f 2`
EarlyDate=`grep "$NAME" $metadata | cut -f 4`
LateDate=`grep "$NAME" $metadata | cut -f 6`
BlockPeriod=`grep "$NAME" $metadata | cut -f 9`

totalstreptococcimappedreadsCM=`grep "$NAME" $competitivemappingsummary | cut -f 6`
col=`grep -n -m 1 "$speciesname" $outdir/headersverticallist.txt | cut -d ":" -f1`
mappedtospecies=`grep "$NAME" $competitivemappingsummary | cut -f$col`
relativeabundance=`echo "scale=3; $mappedtospecies / $totalstreptococcimappedreadsCM" | bc`
  
# Get breadth of coverage data for the islands (single column with name on top)
echo "$NAME" > $outdir/tempbreadth.txt
bedtools coverage -a $outdir/combinedgenomicislandcoordinates.bed -b $PATHS/../../IndividualAlignments/Samples/Sample-"$NAME"/$COVERAGEFOR/"$NAME"_mapped_q30_sorted_rmdup_uniq.bam | cut -f7 >> $outdir/tempbreadth.txt

# Paste it into the summary file. Need to use temporary files to be able to add to themselves.
cat $outdir/"$COVERAGEFOR"_GenomicIslandCoverageSummary.txt > $outdir/tempcombined.txt
paste $outdir/tempcombined.txt $outdir/tempbreadth.txt > $outdir/"$COVERAGEFOR"_GenomicIslandCoverageSummary.txt

# Also add to separate attribute table
echo -e "$NAME \t $breadthcov \t $relativeabundance \t $Date \t $EarlyDate \t $LateDate \t $BlockPeriod" >> $outdir/"$COVERAGEFOR"_SampleAttributes.txt

done < $outdir/"$COVERAGEFOR"_AnalyzableSampleList.txt

rm $outdir/headersverticallist.txt
rm $outdir/tempbreadth.txt
rm $outdir/tempcombined.txt
rm $outdir/tempgenefile.txt

#######################################################################################
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
