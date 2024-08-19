## Usage: parse the output from UniProt ID mapping into 2 tables. The first is the UniProt output table with additional columns about the island the gene with UniProt results is a part of. The second table is a summary table where each island is a row. The gene ontologies for the island are made into non-redundant lists. 

PATHS=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/GenomeAnalysis/GenomicIslands

# Existing output directory
outdir=$PATHS/ByWholeIsland_DD04AncientIndividualMapping-03062024/CharacterizingIslands

islandtable=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/GenomeAnalysis/GenomicIslands/DD04GenomicIslands-03052024.txt

uniprotoutputtable=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/GenomeAnalysis/GenomicIslands/ByWholeIsland_DD04AncientIndividualMapping-03062024/CharacterizingIslands/DD04_idmapping_2024_03_06.tsv

#####################################################
# Add in what island each gene corresponds to in the uniprot results table

header=`cat $uniprotoutputtable | head -n 1`
cat $uniprotoutputtable | grep -v "Entry Name" > $outdir/uniprotoutput_noheader.txt
echo -e "Island \t Identifier \t Gene Locus \t $header " > $outdir/uniprotoutput_withislandnames.txt

while read line
 do
 proteinaccession=`echo "$line" | cut -f 1`
 echo "Adding island for $proteinaccession"
 geneinfo=`grep "$proteinaccession" $islandtable | cut -f1,2,4`
 echo -e "$geneinfo \t $line" >> $outdir/uniprotoutput_withislandnames.txt

done < $outdir/uniprotoutput_noheader.txt

#####################################################
# Create a summery table where each island is a row

# Non-redundant list of islands
cat $islandtable | cut -f 1 | sort | uniq -d > $outdir/listofislands.txt
cat $islandtable | cut -f 1 | sort | uniq -u >> $outdir/listofislands.txt

echo -e "Island \t Overall Number of Protein Coding Genes \t List of Protein Accessions \t List of Gene Names \t List of Gene Names from Uniprot \t Number of Genes with Overall Ontology \t Number of Genes with Biological Process \t Number of Genes with Molecular Function \t Number of Genes with Cellular Component \t Gene Ontology \t Biological Process \t Molecular Function \t Cellular Component" > $outdir/islandontologies_summarytable.txt

# Go by island
while read line
do
numberofgenes=`grep "$line" $islandtable | cut -f 3 | grep -v '^$' | wc -l`
# Counts all, including those without protein accession and pseudogenes: numberofgenes=`grep "$line" $islandtable | wc -l`
listofgenes=`cat -A $islandtable | grep "$line" | cut -f5 -d "^" | tr '\n' ';' | sed 's/;I/; /g' | sed 's/I//' | rev | sed 's/;//' | rev`
listofaccessions=`cat -A $islandtable | grep "$line" | cut -f3 -d "^" | tr '\n' ';' | sed 's/;I/; /g' | sed 's/I//' | rev | sed 's/;//' | rev`
grep "$line" $islandtable | cut -f2 | sort | uniq -d > $outdir/tempidentificationlist.txt
grep "$line" $islandtable | cut -f2 | sort | uniq -u >> $outdir/tempidentificationlist.txt
identification=`cat $outdir/tempidentificationlist.txt | tr '\n' ';' | rev | sed 's/;//' | rev | sed 's/;/; /g'`
grep "$line" $outdir/uniprotoutput_withislandnames.txt | cut -f8 | sort | uniq -d > $outdir/tempidentificationlist.txt
grep "$line" $outdir/uniprotoutput_withislandnames.txt | cut -f8 | sort | uniq -u >> $outdir/tempidentificationlist.txt
listofuniprotgenenames=`cat $outdir/tempidentificationlist.txt | tr '\n' ';' | rev | sed 's/;//' | rev | sed 's/;/; /g'`
#listofuniprotgenenames=`cat -A $outdir/uniprotoutput_withislandnames.txt | grep "$line" | cut -f8 -d "^" | tr '\n' ';' | sed 's/;I/; /g' | sed 's/I//' | rev | sed 's/;//' | rev`

grep "$line" $islandtable > $outdir/tempislandfile.txt

touch $outdir/countingtable.txt
while read wholegeneline
do
gene=`echo "$wholegeneline" | cut -f3`
singlefeature=`grep "$gene" $outdir/uniprotoutput_withislandnames.txt | cut -f12`
 if [ -z "$singlefeature" ]
	then
	echo "No overall ontology for $gene"
	else
	echo "$gene" >> $outdir/countingtable.txt
 fi
numberofgeneswithontology=`cat $outdir/countingtable.txt | grep -v '^$' | wc -l`
done < $outdir/tempislandfile.txt
rm $outdir/countingtable.txt
geneontology=`grep "$line" $outdir/uniprotoutput_withislandnames.txt | cut -f 12 | tr ';' '\n' | sort`
echo "$geneontology" | uniq -d > $outdir/tempontologylist.txt
echo "$geneontology" | uniq -u >> $outdir/tempontologylist.txt
geneontoloynonredundant=`cat $outdir/tempontologylist.txt | sed 's/]/];/g' | tr '\n' ' ' | tr -s " " | rev | sed 's/;//' | rev`

touch $outdir/countingtable.txt
while read wholegeneline
do
gene=`echo "$wholegeneline" | cut -f3`
singlefeature=`grep "$gene" $outdir/uniprotoutput_withislandnames.txt | cut -f13`
 if [ -z "$singlefeature" ]
	then
	echo "No biological process for $gene"
	else
	echo "$gene" >> $outdir/countingtable.txt
 fi
numberofgeneswithbiologicalprocess=`cat $outdir/countingtable.txt | grep -v '^$' | wc -l`
done < $outdir/tempislandfile.txt
rm $outdir/countingtable.txt
biologicalprocess=`grep "$line" $outdir/uniprotoutput_withislandnames.txt | cut -f 13 | tr ';' '\n' | sort`
echo "$biologicalprocess" | uniq -d > $outdir/tempontologylist.txt
echo "$biologicalprocess" | uniq -u >> $outdir/tempontologylist.txt
biologicalprocessnonredundant=`cat $outdir/tempontologylist.txt | sed 's/]/];/g' | tr '\n' ' ' | tr -s " " | rev | sed 's/;//' | rev`

touch $outdir/countingtable.txt
while read wholegeneline
do
gene=`echo "$wholegeneline" | cut -f3`
singlefeature=`grep "$gene" $outdir/uniprotoutput_withislandnames.txt | cut -f15`
 if [ -z "$singlefeature" ]
	then
	echo "No molecular function for $gene"
	else
	echo "$gene" >> $outdir/countingtable.txt
 fi
numberofgeneswithmolecularfunction=`cat $outdir/countingtable.txt | grep -v '^$' | wc -l`
done < $outdir/tempislandfile.txt
rm $outdir/countingtable.txt
molecularfunction=`grep "$line" $outdir/uniprotoutput_withislandnames.txt | cut -f 15 | tr ';' '\n' | sort`
echo "$molecularfunction" | uniq -d > $outdir/tempontologylist.txt
echo "$molecularfunction" | uniq -u >> $outdir/tempontologylist.txt
molecularfunctionnonredundant=`cat $outdir/tempontologylist.txt | sed 's/]/];/g' | tr '\n' ' ' | tr -s " " | rev | sed 's/;//' | rev`

touch $outdir/countingtable.txt
while read wholegeneline
do
gene=`echo "$wholegeneline" | cut -f3`
singlefeature=`grep "$gene" $outdir/uniprotoutput_withislandnames.txt | cut -f16`
 if [ -z "$singlefeature" ]
	then
	echo "No cellular function for $gene"
	else
	echo "$gene" >> $outdir/countingtable.txt
 fi
numberofgeneswithcellularcomponent=`cat $outdir/countingtable.txt | grep -v '^$' | wc -l`
done < $outdir/tempislandfile.txt
rm $outdir/countingtable.txt
cellularcomponent=`grep "$line" $outdir/uniprotoutput_withislandnames.txt | cut -f 16 | tr ';' '\n' | sort`
echo "$cellularcomponent" | uniq -d > $outdir/tempontologylist.txt
echo "$cellularcomponent" | uniq -u >> $outdir/tempontologylist.txt
cellularcomponentnonredundant=`cat $outdir/tempontologylist.txt | sed 's/]/];/g' | tr '\n' ' ' | tr -s " " | rev | sed 's/;//' | rev`

echo -e "$line \t $numberofgenes \t $listofaccessions \t $listofgenes \t $listofuniprotgenenames \t $numberofgeneswithontology \t $numberofgeneswithbiologicalprocess \t $numberofgeneswithmolecularfunction \t $numberofgeneswithcellularcomponent \t $geneontoloynonredundant \t $biologicalprocessnonredundant \t $molecularfunctionnonredundant \t $cellularcomponentnonredundant" >> $outdir/islandontologies_summarytable.txt

done < $outdir/listofislands.txt

rm $outdir/tempidentificationlist.txt
rm $outdir/tempislandfile.txt
rm $outdir/tempontologylist.txt
rm $outdir/uniprotoutput_noheader.txt
rm $outdir/listofislands.txt


