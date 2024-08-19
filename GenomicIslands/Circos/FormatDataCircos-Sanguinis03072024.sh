## Usage: Format read data for plotting genome coverage as a Circos plot. Also plot GC content and CDS of reference. ##

## Input data ##
# FASTA file
# Indexed FASTA file
# BAM files

# Working directory
PATHS=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/GenomeAnalysis/CircosPlot/Streptococcus_sanguinis_SK405_plot-03072024

# faidx file of reference genome with path. Made with reference FASTA and samtools
FAIDX=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/IndividualAlignments/References/Streptococcus_sanguinis/Streptococcus_sanguinis_SK405.fai

# Specify FASTA file with path to calculate GC coverage
FASTAFILE=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/IndividualAlignments/References/Streptococcus_sanguinis/Streptococcus_sanguinis_SK405

# Name of FASTA (or could put anything to describe the genome, just for file naming purposes)
FASTA=Streptococcus_sanguinis_SK405

# Window size to calculate GC content
WIDTH=2500

# Directory with bam file alignments to display
# BAMS=$PATHS/BAMfiles

# Output directory to be made
OUTDIR=$PATHS/SampleCoverages

module use /storage/group/LiberalArts/default/lsw132_collab/sw8/modules
module load bedtools/2.31.1

mkdir $OUTDIR

#################################################
## Format karyotype using information from .fai file ##
# Need columns: chr - ID LABEL START END COLOR
# ID = identifier used in data files
# LABEL = text that will appear next to ideogram on image
# http://www.htslib.org/doc/faidx.html
# Used Excel to get in correct format. All chromosomes/contigs start at 0, and use size as end value.
#################################################
## Histogram to show GC content of reference ##
# Calculated in 2500 bp (or whatever specified $WIDTH) bins when applicable (when contig is >2500 bp; remainder after diving by this is in a smaller averaging window)
# First use faidx file
cat $FAIDX | cut -f1,2 > $PATHS/"$FASTA"_bedtoolsgenome.txt

# Make the 2500 bp windows
bedtools makewindows -g $PATHS/"$FASTA"_bedtoolsgenome.txt -w $WIDTH > $PATHS/"$FASTA"_"$WIDTH"bp.bed

# Use bedtools nuc to calculate GC content
bedtools nuc -fi $FASTAFILE -bed $PATHS/"$FASTA"_"$WIDTH"bp.bed | cut -f1,2,3,5 | grep -v "#" > $PATHS/histogram.GC."$FASTA".txt
#########################################################################################
# Get read data from bam files
for NAME in 16857_Medieval_MertonPriory_MPY86 Yorkshire8893
 do
 
BAM=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/IndividualAlignments/Samples/Sample-"$NAME"/Streptococcus_sanguinis/"$NAME"_mapped_q30_sorted_rmdup_uniq.bam

## Histogram to show coverage depth ##
# Data format: chr start end value [options]
bedtools genomecov -ibam $BAM -bg > $OUTDIR/histogram.coverage."$NAME".txt
###################
## Highlights to show reads ##
# Just remove column with coverage values and have all positions with any coverage
# Data format: chr start end
cat $OUTDIR/histogram.coverage."$NAME".txt | cut -f1,2,3 > $OUTDIR/highlights.reads."$NAME".txt
###################

done

# Do again for healthy and perio modern samples
NAME=19778_ModernCalculus_PerioDisease_12_kneaddata
BAM=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/PerioModernSamples/Sample-19778_ModernCalculus_PerioDisease_12_kneaddata/Streptococcus_sanguinis/19778_ModernCalculus_PerioDisease_12_kneaddata_mapped_q30_sorted_rmdup_uniq.bam

bedtools genomecov -ibam $BAM -bg > $OUTDIR/histogram.coverage."$NAME".txt
cat $OUTDIR/histogram.coverage."$NAME".txt | cut -f1,2,3 > $OUTDIR/highlights.reads."$NAME".txt

NAME=18760_HealthyModern_5_Calculus_kneaddata
BAM=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/CompetitiveMapping/HealthyModernSamples/Sample-18760_HealthyModern_5_Calculus_kneaddata/Streptococcus_sanguinis/18760_HealthyModern_5_Calculus_kneaddata_mapped_q30_sorted_rmdup_uniq.bam

bedtools genomecov -ibam $BAM -bg > $OUTDIR/histogram.coverage."$NAME".txt
cat $OUTDIR/histogram.coverage."$NAME".txt | cut -f1,2,3 > $OUTDIR/highlights.reads."$NAME".txt


########################################################################################
# "Normalize" as needed by changing chromosome names so all uniform (because were changed to include species name as well as accession in competitive mapping) and only including chromosomes of interest (those belonging to the particular species and not all chromosomes for competitive mapping). Thus, this is not necessary for files used in individual mapping.

#cd $OUTDIR

#for file in highlights.reads.SinensisDD04CM.txt highlights.reads.SinensisOnlyCM.txt histogram.coverage.SinensisDD04CM.txt histogram.coverage.SinensisOnlyCM.txt;

 #do
 #NAME=$(basename $file)
 #grep "S.sinensis" $file > tempfile.txt
 #sed -i "s/S.sinensisN/N/g" tempfile.txt
# cp tempfile.txt "$NAME"_normalizedlabels.txt
 #done
# rm tempfile.txt

