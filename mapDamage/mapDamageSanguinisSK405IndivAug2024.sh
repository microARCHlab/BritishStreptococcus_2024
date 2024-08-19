#!/bin/bash
 #SBATCH --nodes=1
 #SBATCH --time=60:00:00
 #SBATCH --mem=200gb
 #SBATCH --account=one
 #SBATCH --output=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/MapDamage/SanguinisOut.txt
 #SBATCH --error=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/MapDamage/SanguinisError.txt

module load anaconda r
# mapDamage2.0
conda activate mapDamage

# Path to alignments
PATHS=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/IndividualAlignments/Samples

# Path to put mapDamage outputs
MAPDAMAGEPATHS=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/MapDamage/Samples

#species=Streptococcus_spDD04
species=Streptococcus_sanguinis

#ref=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/IndividualAlignments/References/Streptococcus_spDD04/Streptococcus_sp._DD04
ref=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/IndividualAlignments/References/Streptococcus_sanguinis/Streptococcus_sanguinis_SK405

outdir="$species"_MapDamage_Aug2024

echo -e "Sample \t Referece \t Analysis-ready_Reads \t Mapped_reads \t Q30_reads \t Unique_reads \t Mean_Length \t Breadth of Coverage \t Depth of Coverage \t Theta_mean \t DeltaD_mean \t DeltaS_mean \t Lamba_mean \t Rho_mean \t LogLik_mean \t Theta_std \t DeltaD_std \t DeltaS_std \t Lambda_std \t Rho_std \t LogLik_std" > $MAPDAMAGEPATHS/../"$species"_Results.txt

for sample in $PATHS/Sample-*;
 do
  cd $sample

  ID=$(basename $sample)
  NAME=$(echo $ID |cut -d "-" -f2)
  echo "Running pipeline on sample "${NAME}"..."
  
  mkdir $MAPDAMAGEPATHS/Sample-"$NAME"/$outdir

mapDamage -i $PATHS/Sample-"${NAME}"/"$species"/"${NAME}"_mapped_q30_sorted_rmdup_uniq.bam -r $ref -d $MAPDAMAGEPATHS/Sample-"$NAME"/$outdir --merge-libraries

AnalysisReady=`grep -A 1 "Analysis-ready reads" $PATHS/Sample-"${NAME}"/$species/"${NAME}"_MappingReport.txt | grep -v "Analysis-ready reads"`
Mapped=`grep -A 1 "Total mapped reads" $PATHS/Sample-"${NAME}"/$species/"${NAME}"_MappingReport.txt | grep -v "Total mapped reads"`
Q30Mapped=`grep -A 1 "Total Q30 mapped reads" $PATHS/Sample-"${NAME}"/$species/"${NAME}"_MappingReport.txt | grep -v "Total Q30 mapped reads"`
Unique=`grep -A 1 "Total unique reads" $PATHS/Sample-"${NAME}"/$species/"${NAME}"_MappingReport.txt | grep -v "Total unique reads"`
MeanLength=`grep -A 1 "Average length of reads" $PATHS/Sample-"${NAME}"/$species/"${NAME}"_MappingReport.txt | grep -v "Average length of reads"`
BreadthofCoverage=`grep -A 1 "Breadth of coverage" $PATHS/Sample-"${NAME}"/$species/"${NAME}"_MappingReport.txt | grep -v "Breadth of coverage"`
DepthofCoverage=`grep -A 1 "Depth of coverage" $PATHS/Sample-"${NAME}"/$species/"${NAME}"_MappingReport.txt | grep -v "Depth of coverage"`

ThetaMean=`grep "Mean" $MAPDAMAGEPATHS/Sample-"${NAME}"/$outdir/Stats_out_MCMC_iter_summ_stat.csv | cut -d "," -f2`
DeltaDMean=`grep "Mean" $MAPDAMAGEPATHS/Sample-"${NAME}"/$outdir/Stats_out_MCMC_iter_summ_stat.csv | cut -d "," -f3`
DeltaSMean=`grep "Mean" $MAPDAMAGEPATHS/Sample-"${NAME}"/$outdir/Stats_out_MCMC_iter_summ_stat.csv | cut -d "," -f4`
LambaMean=`grep "Mean" $MAPDAMAGEPATHS/Sample-"${NAME}"/$outdir/Stats_out_MCMC_iter_summ_stat.csv | cut -d "," -f5`
RhoMean=`grep "Mean" $MAPDAMAGEPATHS/Sample-"${NAME}"/$outdir/Stats_out_MCMC_iter_summ_stat.csv | cut -d "," -f6`
LogLikMean=`grep "Mean" $MAPDAMAGEPATHS/Sample-"${NAME}"/$outdir/Stats_out_MCMC_iter_summ_stat.csv | cut -d "," -f7`

ThetaStd=`grep "Std." $MAPDAMAGEPATHS/Sample-"${NAME}"/$outdir/Stats_out_MCMC_iter_summ_stat.csv | cut -d "," -f2`
DeltaDStd=`grep "Std." $MAPDAMAGEPATHS/Sample-"${NAME}"/$outdir/Stats_out_MCMC_iter_summ_stat.csv | cut -d "," -f3`
DeltaSStd=`grep "Std." $MAPDAMAGEPATHS/Sample-"${NAME}"/$outdir/Stats_out_MCMC_iter_summ_stat.csv | cut -d "," -f4`
LambdaStd=`grep "Std." $MAPDAMAGEPATHS/Sample-"${NAME}"/$outdir/Stats_out_MCMC_iter_summ_stat.csv | cut -d "," -f5`
RhoStd=`grep "Std." $MAPDAMAGEPATHS/Sample-"${NAME}"/$outdir/Stats_out_MCMC_iter_summ_stat.csv | cut -d "," -f6`
LogLikStd=`grep "Std." $MAPDAMAGEPATHS/Sample-"${NAME}"/$outdir/Stats_out_MCMC_iter_summ_stat.csv | cut -d "," -f7`

cd ..

# Printing to file
echo -e "$NAME \t $species \t $AnalysisReady \t $Mapped \t $Q30Mapped \t $Unique \t $MeanLength \t $BreadthofCoverage \t $DepthofCoverage \t $ThetaMean \t $DeltaDMean \t $DeltaSMean \t $LambaMean \t $RhoMean \t $LogLikMean \t $ThetaStd \t $DeltaDStd \t $DeltaSStd \t $LambdaStd \t $RhoStd \t $LogLikStd" >> $MAPDAMAGEPATHS/../"$species"_Results.txt


done
