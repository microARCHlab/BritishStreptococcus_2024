# Check for heterozygosity (indicating potential mixed strain mapping) in the variants called against the S.sp.DD04 (S.sinensis) genome that were used to create consensus sequences and create a phylogeny.
# Create summary tables of the number of mixed variants (allele frequency between 20 and 80%) compared to the number of homozygous variants (allele frequency ≥80%), including average depth.
# I also repeated the analysis comparing alelle frequencies between 40 and 60% and ≥80%.
# Calculate separately for variants with higher (≥x3) or any coverage. Lower coverages are more likely to give heterozygous variants with sequencing error or deamination.

PATHS=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/Phylogenetics
listofsamples=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/Phylogenetics/SspDD04Samples.txt

####################################################################################
# Start summary files, one for each where variants are filtered by depth (≥3) and not
echo -e "Sample \t Count Homozygous \t Count Heterozygous \t Average Depth Homozygous \t Average Depth Heterozygous" > $PATHS/DD04variantscheckheterozygous_3cov.txt
echo -e "Sample \t Count Homozygous \t Count Heterozygous \t Average Depth Homozygous \t Average Depth Heterozygous" > $PATHS/DD04variantscheckheterozygous.txt

while read Sample
 do
echo "Working on $Sample"
 cd /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/IndividualAlignments/Samples/Sample-"$Sample"/Streptococcus_spDD04_Strict/Streptococcus_spDD04_Strict_VariantCallingandConsensus

gzip -d "$Sample"_calls.vcf.gz
 
# Be able to go line-by-line; remove vcf header
grep -v "#" "$Sample"_calls.vcf > "$Sample"_calls_noheader.txt

# Prepare temporary files
touch $PATHS/tempheterozygous3cov.txt
touch $PATHS/temphomzygous3cov.txt
touch $PATHS/tempheterozygous.txt
touch $PATHS/temphomozygous.txt
 
 while read line # Unfiltered VCF
 do
# For each line extract VAF and AD
# AD = number of high quality bases, given for each allele (so add together)
# Bash does not work well with decimals, so will need to convert them. But in the case of 100% allele frequency, it is given as 1, which will not convert the same way. All others covert to tenths place (for example, .90 or 90% –> 9; .92 also –> 9)
VAFcheck=`echo $line | rev | cut -f1 -d ":" | rev`
if [ "$VAFcheck" == 1 ]
	then
	VAF=10
	else
	VAF=`echo $line | rev | cut -f1 -d ":" | rev | cut -d "." -f2 | grep -Po "^."`
fi
 AD1=`echo $line | rev | cut -f2 -d ":" | rev | cut -d "," -f1`
 AD2=`echo $line | rev | cut -f2 -d ":" | rev | cut -d "," -f2`
 AD=`echo $AD1 + $AD2 | bc`
 
# Work with if AD ≥ 3 (at least 3 allele depth) separately
if [ "$AD" -ge 3 ]
 then
 echo "Has ≥ 3 coverage"
# If greater than 20% or less than 80%, record as heterozygous
if [ "$VAF" -gt 2 ] && [ "$VAF" -lt 8 ]
 then
 echo "$AD" >> $PATHS/tempheterozygous3cov.txt
 echo "$AD" >> $PATHS/tempheterozygous.txt
fi
# If greater than 80%, record as homozygous.
if [ "$VAF" -ge 8 ]
then
 echo "$AD" >> $PATHS/temphomzygous3cov.txt
 echo "$AD" >> $PATHS/temphomozygous.txt
fi

else # Now also calculate everything else
echo "Does not have ≥ 3 coverage"
# If greater than 20% or less than 80%, record as heterozygous
if [ "$VAF" -gt 2 ] && [ "$VAF" -lt 8 ]
 then
 echo "$AD" >> $PATHS/tempheterozygous.txt
fi
# If greater than 80%, record as homozygous.
if [ "$VAF" -ge 8 ]
then
 echo "$AD" >> $PATHS/temphomozygous.txt
fi

fi

done < "$Sample"_calls_noheader.txt

#####################
# Calculate the total numbers of each type for at least 3x coverage
counthet3cov=`cat $PATHS/tempheterozygous3cov.txt | wc -l`
counthom3cov=`cat $PATHS/temphomzygous3cov.txt | wc -l`
# Calculate average depth
totdepthhet3cov=`cat $PATHS/tempheterozygous3cov.txt | tr '\n' '+' | tr -d " " | rev | cut -d "+" -f 2- | rev | bc`
avdepthhet3cov=`echo "scale =3; $totdepthhet3cov/$counthet3cov" | bc -l`
totdepthhom3cov=`cat $PATHS/temphomzygous3cov.txt | tr '\n' '+' | tr -d " " | rev | cut -d "+" -f 2- | rev | bc`
avdepthhom3cov=`echo "scale =3; $totdepthhom3cov/$counthom3cov" | bc -l`
# Export to summary file
echo -e "$Sample \t $counthom3cov \t $counthet3cov \t $avdepthhom3cov \t $avdepthhet3cov" >> $PATHS/DD04variantscheckheterozygous_3cov.txt

# Calculate the total numbers of each type (for any coverage)
counthet=`cat $PATHS/tempheterozygous.txt | wc -l`
counthom=`cat $PATHS/temphomozygous.txt | wc -l`
# Calculate average depth
totdepthhet=`cat $PATHS/tempheterozygous.txt | tr '\n' '+' | tr -d " " | rev | cut -d "+" -f 2- | rev | bc`
avdepthhet=`echo "scale =3; $totdepthhet/$counthet" | bc -l`
totdepthhom=`cat $PATHS/temphomozygous.txt | tr '\n' '+' | tr -d " " | rev | cut -d "+" -f 2- | rev | bc`
avdepthhom=`echo "scale =3; $totdepthhom/$counthom" | bc -l`
# Export to summary file
echo -e "$Sample \t $counthom \t $counthet \t $avdepthhom \t $avdepthhet" >> $PATHS/DD04variantscheckheterozygous.txt

rm $PATHS/tempheterozygous.txt
rm $PATHS/temphomozygous.txt
rm $PATHS/tempheterozygous3cov.txt
rm $PATHS/temphomzygous3cov.txt

done < $listofsamples
