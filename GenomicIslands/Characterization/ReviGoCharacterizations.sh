## Usage: try to reduce gene ontology lists for the islands produced by UniProt using Revigo (downloaded R scripts). Input needs the Uniprot raw results with island attribute added to each result (see past script). ##

# PATHS=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/GenomeAnalysis/GenomicIslands/ByWholeIsland_SanguinisIndividualMapping-03062024
PATHS=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/GenomeAnalysis/GenomicIslands/ByWholeIsland_DD04AncientIndividualMapping-03062024

# Existing output directory for formatted ontology lists
outdir=$PATHS/CharacterizingIslands/OntologyLists

islandtable=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/GenomeAnalysis/GenomicIslands/DD04GenomicIslands-03052024.txt
#islandtable=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/GenomeAnalysis/GenomicIslands/SanguinisGenomicIslands–03062024.txt

# for naming
NAME=DD04

# Directory with island files (from coverage script)
#islandfiles=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/GenomeAnalysis/GenomicIslands/ByWholeIsland_SanguinisIndividualMapping-03062024/GenomicIslands
islandfiles=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/GenomeAnalysis/GenomicIslands/ByWholeIsland_DD04AncientIndividualMapping-03062024/GenomicIslands

# List of filtered islands (above 2000 bp cut off)
listofislands=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/GenomeAnalysis/GenomicIslands/ByWholeIsland_DD04AncientIndividualMapping-03062024/Filtered_islands_2000bp.txt
#listofislands=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/GenomeAnalysis/GenomicIslands/ByWholeIsland_SanguinisIndividualMapping-03062024/Filtered_islands_2000bp.txt

# UniProt results with island names
#uniprot=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/GenomeAnalysis/GenomicIslands/ByWholeIsland_SanguinisIndividualMapping-03062024/CharacterizingIslands/uniprotoutput_withislandnames.txt
uniprot=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/GenomeAnalysis/GenomicIslands/ByWholeIsland_DD04AncientIndividualMapping-03062024/CharacterizingIslands/uniprotoutput_withislandnames.txt

module load r

# Start summary file
echo -e "Genomic Island \t Total Genes \t Total Genes with Protein Accessions \t Total Genes with Assigned GO Terms \t Genes with Multiple UniProt Hits \t Genes Where Multiple Uniprot Hits Added GO Terms \t Fraction of CDS Genes with GO Terms \t All GO Terms" > $PATHS/"$NAME"_ontologysummarywithReviGoinput.txt

# Remove genes from the UniProt results that do not have ontologies
touch $PATHS/"$NAME"_allwithontologies.txt
while read line
 do
 ontologycolumn=`echo "$line" | cut -f14`
 if [[ -z $ontologycolumn ]]
 then 
 echo "Ontology column is empty"
 else 
	 echo "$line" >> $PATHS/"$NAME"_allwithontologies.txt
 fi
done < $uniprot

# Count how many genes were in each of these islands to start with
while read Island
do
island=`echo "$Island" | tr -d '\r'`
mkdir $outdir/"$island"
echo "Summarizing for $island"
islandfile="$islandfiles"/"$island".txt
totalgenes=`cat $islandfile | wc -l`
# Genes with protein accessions were input in uniprot
totalgeneswithaccessions=`cat $islandfile | cut -f3 | grep -v '^$' | wc -l`

# Count how many genes from the island have GO terms in the UniProt results (make non-redundant list)
grep "$island" $PATHS/"$NAME"_allwithontologies.txt | cut -f4 | sort > $PATHS/tempaccessionlist.txt
cat $PATHS/tempaccessionlist.txt | sort | uniq -u > $PATHS/tempaccessionlist2.txt
cat $PATHS/tempaccessionlist.txt | sort | uniq -d >> $PATHS/tempaccessionlist2.txt
totalgenesgoterms=`cat $PATHS/tempaccessionlist2.txt | wc -l`

# Count how many genes from the island have multiple UniProt results
cat $PATHS/tempaccessionlist.txt | sort | uniq -d > $PATHS/multipleuniprotresults.txt
multipleuniprotresults=`cat $PATHS/multipleuniprotresults.txt | wc -l`

# Count how many of the genes with multiple UniProt results have new gene ontologies as a result (basically, do these multiple results have duplicate ontologies?)
touch $PATHS/countingfile.txt
while read accession
 do
 grep "$accession" $PATHS/"$NAME"_allwithontologies.txt > $PATHS/linesforgene.txt
 grep "$accession" $PATHS/"$NAME"_allwithontologies.txt | cut -f14 | tr ';' '\n' | tr -d " " > $PATHS/tempaccessionlist3.txt
# If there are uniq values, are they all from the same line (UniProt hit)?
cat $PATHS/tempaccessionlist3.txt | sort | uniq -u > $PATHS/uniqueGOs.txt
while read GO
do
grep "$GO" $PATHS/linesforgene.txt | cut -f1 > $PATHS/entriesforgene.txt
done < $PATHS/uniqueGOs.txt
unique=`cat $PATHS/entriesforgene.txt | sort | uniq -u`
# If there were multiple entries contributing to this list of unique genes, that means unique GO terms came from having additional entries.
if [[ -z $unique ]]
 then 
 echo "$accession results do not have multiple hits contributing to multiple ontologies"
 else
 echo "$accession" >> $PATHS/countingfile.txt
fi
done < $PATHS/multipleuniprotresults.txt
numberofgeneswithmultipleontologiesduetomultiplehits=`cat $PATHS/countingfile.txt | wc -l`
numberofgeneswithmultipleontologiesduetomultiplehitsnumber=`echo $numberofgeneswithmultipleontologiesduetomultiplehits`
rm $PATHS/countingfile.txt

# Get a list of GO terms for the island formatted as a list
grep "$island" $PATHS/"$NAME"_allwithontologies.txt | cut -f14 | tr ';' '\n' | sort | uniq -u | tr -d " " > $outdir/"$island"/"$island"_ontologieslist.txt
grep "$island" $PATHS/"$NAME"_allwithontologies.txt | cut -f14 | tr ';' '\n' | sort | uniq -d | tr -d " " >> $outdir/"$island"/"$island"_ontologieslist.txt

GOterms=`cat $outdir/"$island"/"$island"_ontologieslist.txt | tr '\n' ';'`

# Add to summary file
echo -e "$island \t $totalgenes \t $totalgeneswithaccessions \t $totalgenesgoterms \t $multipleuniprotresults \t $numberofgeneswithmultipleontologiesduetomultiplehitsnumber \t $totalgenesgoterms/$totalgeneswithaccessions \t $GOterms" >> $PATHS/"$NAME"_ontologysummarywithReviGoinput.txt

#rm "$NAME"_allwithontologies.txt

## Running revigo

# Run for biological processes
echo "Running ReviGo for biological processes for $island"
Rscript /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/GenomeAnalysis/GenomicIslands/ReviGoScripts/revigo-bioprocesses.R $outdir/"$island"/"$island"_ontologieslist.txt $outdir/"$island"/"$island"_tempprocesses.txt
touch $outdir/"$island"/"$island"_biologicalprocesses.txt
while read line
 do
 echo -e "$island \t Biological Process \t $line" >> $outdir/"$island"/"$island"_biologicalprocesses.txt
done < $outdir/"$island"/"$island"_tempprocesses.txt

# Run for molecular function
echo "Running ReviGo for molecular function for $island"
Rscript /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/GenomeAnalysis/GenomicIslands/ReviGoScripts/revigo-molecularfunction.R $outdir/"$island"/"$island"_ontologieslist.txt $outdir/"$island"/"$island"_tempprocesses.txt
touch $outdir/"$island"/"$island"_molecularfunction.txt
while read line
 do
 echo -e "$island \t Molecular Function \t $line" >> $outdir/"$island"/"$island"_molecularfunction.txt
done < $outdir/"$island"/"$island"_tempprocesses.txt

# Run for cellular component
echo "Running ReviGo for cellular component for $island"
Rscript /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/GenomeAnalysis/GenomicIslands/ReviGoScripts/revigo-cellularcomponent.R $outdir/"$island"/"$island"_ontologieslist.txt $outdir/"$island"/"$island"_tempprocesses.txt
touch $outdir/"$island"/"$island"_cellularcomponent.txt
while read line
 do
 echo -e "$island \t Cellular Component \t $line" >> $outdir/"$island"/"$island"_cellularcomponent.txt
done < $outdir/"$island"/"$island"_tempprocesses.txt

rm $outdir/"$island"/"$island"_tempprocesses.txt

# Make island summary file
echo -e "TermID \t Name \t Value \t LogSize \t Frequency \t Uniqueness \t Dispensability \t Representative" > $outdir/"$island"/"$island"_revigocombinedresults.txt 
for file in $outdir/"$island"/"$island"_biologicalprocesses.txt $outdir/"$island"/"$island"_molecularfunction.txt $outdir/"$island"/"$island"_cellularcomponent.txt
do
grep -v "TermID" $file >> $outdir/"$island"/"$island"_revigocombinedresults.txt 
done

# Make summary line (will be used to make summary file with all islands)
biologicalprocesses=`grep -v "TermID" $outdir/"$island"/"$island"_biologicalprocesses.txt | cut -f3,4 | tr -d '"' | tr '\t' '-' | tr '\n' ';'`
molecularfunction=`grep -v "TermID" $outdir/"$island"/"$island"_molecularfunction.txt | cut -f3,4 | tr -d '"' | tr '\t' '-' | tr '\n' ';'`
cellularcomponent=`grep -v "TermID" $outdir/"$island"/"$island"_cellularcomponent.txt | cut -f3,4 | tr -d '"' | tr '\t' '-' | tr '\n' ';'`
echo -e "Island \t Biological Processes \t Molecular Function \t Cellular Component" > $outdir/"$island"/"$island"_summaryline.txt
echo -e "$island \t $biologicalprocesses \t $molecularfunction \t $cellularcomponent" >> $outdir/"$island"/"$island"_summaryline.txt

done < $listofislands

# Join results from all islands
echo -e "TermID \t Name \t Value \t LogSize \t Frequency \t Uniqueness \t Dispensability \t Representative" > $outdir/"$NAME"_revigocombinedresults.txt 
for file in $outdir/GI*/*_revigocombinedresults.txt
do
grep -v "TermID" $file >> $outdir/"$NAME"_revigocombinedresults.txt 
done

# Join summary lines from all islands
echo -e "Island \t Biological Processes \t Molecular Function \t Cellular Component" > $outdir/"$NAME"_revigocombinedresults_summarized.txt
for file in $outdir/GI*/*_summaryline.txt
do
grep -v "Island" $file >> $outdir/"$NAME"_revigocombinedresults_summarized.txt
done

# Join this summary table with the one from before
header1=`head -n 1 "$NAME"_ontologysummarywithReviGoinput.txt`
header2=`head -n 1 $outdir/"$NAME"_revigocombinedresults_summarized.txt | cut -f2-`
echo -e "$header1 \t Accessions \t Gene Names \t UniProt Names \t $header2" > $outdir/"$NAME"_ALLontologysummarizedrevigo.txt

while read Island
do
island=`echo "$Island" | tr -d '\r'`
revigo=`grep "$island" $outdir/"$NAME"_revigocombinedresults_summarized.txt | cut -f2-`
ontologysummary=`grep $island "$NAME"_ontologysummarywithReviGoinput.txt`
islandfile="$islandfiles"/"$island".txt
listofaccessions=`cat $islandfile | cut -f3 | tr '\n' ';'`
listofgenes=`cat -A $islandtable | grep "$island" | cut -f5 -d "^" | tr '\n' ';' | sed 's/;I/; /g' | sed 's/I//' | rev | sed 's/;//' | rev`
cat $islandfile | cut -f5 | tr '\r' ';' | tr -d '/n' > temp.txt
listofunprotgenenames=`grep "$island" $uniprot | cut -f8 | tr '\n' ';' | rev | sed 's/;//' | rev | sed 's/;/; /g'`

echo -e "$ontologysummary \t $listofaccessions \t $listofgenes \t $listofunprotgenenames \t $revigo" >> $outdir/"$NAME"_ALLontologysummarizedrevigo.txt
done < $listofislands

rm $PATHS/entriesforgene.txt
rm $PATHS/tempaccessionlist.txt
rm $PATHS/tempaccessionlist2.txt
rm $PATHS/tempaccessionlist3.txt
rm $PATHS/uniqueGOs.txt
rm $PATHS/linesforgene.txt
rm $PATHS/multipleuniprotresuts.txt


