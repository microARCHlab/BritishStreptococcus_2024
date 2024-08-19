## Usage: Condense circos read (highlight) coordinates to ultimately limit file size by limiting amount plotted. Overlapping highlights are joined into one single highlight. ##

# For S. sanguinis
# Path where to put output files
outdir=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/GenomeAnalysis/CircosPlot/Streptococcus_sanguinis_SK405_plot-03072024/SampleCoverages/CompressedHighlights

mkdir $outdir

# Assumes data is already sorted by contig/chromosome and then in position order, as it would be following prior script
for file in /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/GenomeAnalysis/CircosPlot/Streptococcus_sanguinis_SK405_plot-03072024/SampleCoverages/highlights.reads.*.txt
 do
 outname=`echo "$file" | rev | cut -d "/" -f1 | cut -d "." -f2 | rev`
 #outfilepath=`echo "$file" | rev | cut -d "/" -f2- | rev`
 outfile="$outdir"/"$outname".txt
 echo "Condensing for $file"
 touch $outfile
while read line
 do
 	contig=`echo "$line" | cut -f1 | tr -d " "`
 	contigbefore=`tail -n 1 $outfile | cut -f1 | tr -d " "`
 	startpos=`echo "$line" | cut -f2 | tr -d " "`
 	startposbefore=`tail -n 1 $outfile | cut -f2 | tr -d " "`
 	endpos=`echo "$line" | cut -f3 | tr -d " "`
 	endposbefore=`tail -n 1 $outfile | cut -f3 | tr -d " "`
 	if [[ "$contig" == "$contigbefore" ]] 
 	 then
 	 	if [[ "$startpos" == "$endposbefore" ]]
 	 	 then
 	 	 # If new "read" starts at the same position as the one before it, make the new "read" have the new end position and but starting position from before; remove the "incomplete" "read" that was already output to the file before to be able to add the full "read" to the output file
# Else, add the "read" as is to the new file
 	 	  finalendpos=`echo "$endpos"`
 	 	  finalstartpos=`echo $startposbefore`
 	 	  sed '$d' $outfile > $outdir/tempfile.txt
 	 	  mv $outdir/tempfile.txt $outfile
 	 	  else
 	 	  finalendpos=`echo "$endpos"`
 	 	  finalstartpos=`echo "$startpos"`
 	 	  fi
 	 	 echo -e "$contig \t $finalstartpos \t $finalendpos" >> $outfile
 	 # If this is the start of the new contig, then there is nothing to compare to and just add the first "read" to the outfile
 	 else
 	 	echo -e "$contig \t $startpos \t $endpos" >> $outfile
 	 	fi
 done < $file
 done
 
rm $outdir/tempfile.txt

 