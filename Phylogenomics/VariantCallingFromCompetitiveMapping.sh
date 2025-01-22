# Usage: make a phylogeny with sequences competitively mapped. Compare this tree to the same samples individually mapped.

PATHS=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/Phylogenetics/Mugsy/SspDD04/SanguinisOutgroupAdditionalReferences/CompetitiveMappingSequences

# List of samples to run on
listofsamples=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/Phylogenetics/SspDD04Samples.txt

# Output directory for consensus sequences
outdir=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/Phylogenetics/ConsensusSequences/AllSpeciesCompetitiveMapping

# Index file for combined reference used in competitive mapping
FAIDX=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/JoinedReferences/NewReference/NewReference.fna.fai

# Common identifier in all contigs/chromosomes for species of interest
identifer=S.spDD04

# Species name (for file naming purposes)
species=DD04

# INDIVIDUAL reference for species of interest. The contig/chromosome names must match those in the combined reference exactly
REF=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/JoinedReferences/NewReference/GCF_001578725.1_ASM157872v1_genomic.fna

# bcftools and setGTplugin (allows filtering of variants)
bcftools=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/Phylogenetics/snpToolkit/bcftools/installed/bin/bcftools
setGTplugin=+/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/Phylogenetics/snpToolkit/bcftools/installed/libexec/bcftools/setGT.so
filltagsplugin=+/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/Phylogenetics/snpToolkit/bcftools/installed/libexec/bcftools/fill-tags.so

# Load modules; conda environment needed for bcftools here
module load samtools/1.18 anaconda bedtools/2.31.0
conda activate snptoolkit

#######################################################################################
## Getting the reads only mapped to the species of interest ##

# Use faidx file of combined reference to get names of chromosome/contigs/scaffolds of interest. These were named with appropriate species identifier
grep "S.spDD04" $FAIDX | cut -f1 > $PATHS/SspDD04Regions.txt

while read Sample
 do
 echo "Beginning process for $Sample"
 
# Extract these regions from the sorted filtered bam file for selected samples.
BAM=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples/Sample-"$Sample"/AllSpecies/"$Sample"_SAMandBAMfiles/"$Sample"_mapped_q30_sorted_rmdup_uniq.bam

samtools index $BAM # BAM files need to be indexed
cat $PATHS/SspDD04Regions.txt | tr "\n" " " | xargs samtools view -bh $BAM > /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples/Sample-"$Sample"/AllSpecies/"$Sample"_SAMandBAMfiles/"$Sample"_filtered_"$species".bam

NEWBAM=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples/Sample-"$Sample"/AllSpecies/"$Sample"_SAMandBAMfiles/"$Sample"_filtered_"$species".bam

##############################################
## Call variants ##

# Prepare summary table
echo -e "Sample \t Total called variants \t Variants after filtering \t" > $PATHS/"$species"_CompetitiveMappingVariantCallingSummary.txt

# Use bedtools genomecov to calculate coverage at each base (-d) with extracted sorted bam input (ibam).
# Remove lines where coverage (column 3) is greater than or equal to 1 (covered bases). The remaining bases that do not meet this criteria (not covered) will be masked in the consensus sequence.
 echo "Finding regions to mask for $Sample"
bedtools genomecov -d -ibam $NEWBAM | awk '$3 < 1' | cut -f1,2 > /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples/Sample-"$Sample"/AllSpecies/"$Sample"_SAMandBAMfiles/"$Sample"_"$species"_regionsttomask.txt

# Call variants using bcftools.
# The annotation list in mpileup (-a) annotates variants with allelic depth (AD) and number of high-quality bases/overall depth (DP), allowing for filtering.
# Output only variant sites (-v), use the multiallelic-caller (-m), and output compressed BCF (-Ob).
# The fill-in tags plugin calculates and annotates fraction of reads with alternate allele (VAF) for each allele; Oz for compressed VCF output. 
echo "Calling variants for $Sample"
$bcftools mpileup -a AD,DP,SP -Ou -f $REF $NEWBAM | $bcftools call --ploidy 1 -mv |  $bcftools $filltagsplugin -Oz -o /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples/Sample-"$Sample"/AllSpecies/"$Sample"_SAMandBAMfiles/"$Sample"_"$species"_calls.vcf.gz -- -t VAF

# Check how many variants there are before filtering (-H does not include header)
VariantsBeforeFiltering=`$bcftools view -H /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples/Sample-"$Sample"/AllSpecies/"$Sample"_SAMandBAMfiles/"$Sample"_"$species"_calls.vcf.gz | wc -l`

# Filter variants
# Here require a total allelic depth greater than or equal to 3 (FORMAT/AD>=3), quality greater than or equal to 30 (QUAL>30), and variant allele frequency greater than or equal to 0.90 (FORMAT/VAF>=0.90).
$bcftools view -i 'FORMAT/AD>=3 && QUAL>30 && FORMAT/VAF>=0.90' /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples/Sample-"$Sample"/AllSpecies/"$Sample"_SAMandBAMfiles/"$Sample"_"$species"_calls.vcf.gz | $bcftools view -Oz -o /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples/Sample-"$Sample"/AllSpecies/"$Sample"_SAMandBAMfiles/"$Sample"_"$species"_filtered.vcf.gz

# Check how many variants there are after filtering
VariantsAfterFiltering=`$bcftools view -H /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples/Sample-"$Sample"/AllSpecies/"$Sample"_SAMandBAMfiles/"$Sample"_"$species"_filtered.vcf.gz | wc -l`

# Index VCF
$bcftools index /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples/Sample-"$Sample"/AllSpecies/"$Sample"_SAMandBAMfiles/"$Sample"_"$species"_filtered.vcf.gz

# Generate consensus sequence with masked regions
echo "Creating consensus sequence for $Sample"
cat $REF | $bcftools consensus -m /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples/Sample-"$Sample"/AllSpecies/"$Sample"_SAMandBAMfiles/"$Sample"_"$species"_regionsttomask.txt /storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/Samples/Sample-"$Sample"/AllSpecies/"$Sample"_SAMandBAMfiles/"$Sample"_"$species"_filtered.vcf.gz > $outdir/"$Sample"_"$species"_competitivemappingconsensus.fa

##############################################
# Print results to summary file
echo -e "$Sample \t $VariantsBeforeFiltering \t $VariantsAfterFiltering \t" >> $PATHS/"$species"_CompetitiveMappingVariantCallingSummary.txt

done < $listofsamples
