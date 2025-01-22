### FEB 2023 ###

## This script performs variant calling and filtering before producing a consensus sequence from a sorted and filtered BAM file. Regions not covered are masked with an N in the consensus sequence. A summary table of variants before and after filtering is also provided. ##

# Reference for variant calling must be same used in generation of bam files during alignment
REF=/gpfs/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/IndividualAlignments/References/Streptococcus_spDD04/Streptococcus_sp._DD04

# bcftools program and plugins files (plugin files downloaded automatically as part of bcftools)
bcftools=/gpfs/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/Phylogenetics/bcftools/bin/bcftools
setGTplugin=+/gpfs/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/Phylogenetics/bcftools/libexec/bcftools/setGT.so
filltagsplugin=+/gpfs/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/Phylogenetics/bcftools/libexec/bcftools/fill-tags.so

# Path to directories with samples
PATHS=/gpfs/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/IndividualAlignments/Samples

# Species (subdirectory name in sample folder)
SPECIES=Streptococcus_spDD04

# Specify name of out directory (within each sample folder)
OUTDIR="${SPECIES}"_VariantCallingandConsensus

# Specify name of directory for summary table output (only)
SUMOUT=/gpfs/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/Phylogenetics

# File with list of sample names to create consensus sequences for
LIST=/gpfs/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/Phylogenetics/SspDD04Samples.txt

module load gcc samtools bedtools

# Directory set up:
	# Samples directory
		# Individual sample directories
			# Alignment files (sorted and filtered bam)

##################################

# Prepare summary table
echo -e "Sample \t Total called variants \t Variants after filtering \t" > $SUMOUT/"${SPECIES}"_VariantCallingSummary.txt
  
##################################

while read line;
 do
  NAME=`echo $line`
  echo "Creating consensus sequence for "${NAME}""
  cd $PATHS/Sample-"${NAME}"/"${SPECIES}"
  BAM="${NAME}"_mapped_q30_sorted_rmdup_uniq.bam
  mkdir $OUTDIR
  
##################################
# Use bedtools genomecov to calculate coverage at each base (-d) with sorted bam input (ibam).
# Remove lines where coverage (column 3) is greater than or equal to 1 (covered bases). The remaining bases that do not meet this criteria (not covered) will be masked in the consensus sequence.
bedtools genomecov -d -ibam $BAM | awk '$3 < 1' | cut -f1,2 > $OUTDIR/"${NAME}"_regionstomask.txt

# Call variants using bcftools.
# The annotation list in mpileup (-a) annotates variants with allelic depth (AD) and number of high-quality bases/overall depth (DP), allowing for filtering.
# Output only variant sites (-v), use the multiallelic-caller (-m), and output compressed BCF (-Ob).
# The fill-in tags plugin calculates and annotates fraction of reads with alternate allele (VAF) for each allele; Oz for compressed VCF output. 
$bcftools mpileup -a AD,DP,SP -Ou -f $REF $BAM | $bcftools call --ploidy 1 -mv |  $bcftools $filltagsplugin -Oz -o $OUTDIR/"${NAME}"_calls.vcf.gz -- -t VAF

# Check how many variants there are before filtering (-H does not include header)
VariantsBeforeFiltering=`$bcftools view -H $OUTDIR/"${NAME}"_calls.vcf.gz | wc -l`

# Filter variants
# Here require a total allelic depth greater than or equal to 3 (FORMAT/AD>=3), quality greater than or equal to 30 (QUAL>30), and variant allele frequency greater than or equal to 0.90 (FORMAT/VAF>=0.90).
$bcftools view -i 'FORMAT/AD>=3 && QUAL>30 && FORMAT/VAF>=0.90' $OUTDIR/"${NAME}"_calls.vcf.gz | $bcftools view -Oz -o $OUTDIR/"${NAME}"_filtered.vcf.gz

# Check how many variants there are after filtering
VariantsAfterFiltering=`$bcftools view -H $OUTDIR/"${NAME}"_filtered.vcf.gz | wc -l`

# Index VCF
$bcftools index $OUTDIR/"${NAME}"_filtered.vcf.gz

# Generate consensus sequence with masked regions
cat $REF | $bcftools consensus -m $OUTDIR/"${NAME}"_regionstomask.txt $OUTDIR/"${NAME}"_filtered.vcf.gz > $OUTDIR/"${NAME}"_"${SPECIES}"_consensus.fa

##################################

# Print results to summary file
echo -e "${NAME} \t $VariantsBeforeFiltering \t $VariantsAfterFiltering \t" >> $SUMOUT/"${SPECIES}"_VariantCallingSummary.txt

bgzip -d $OUTDIR/"${NAME}"_filtered.vcf.gz

done < $LIST