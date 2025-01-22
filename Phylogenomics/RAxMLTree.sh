## Usage: run a phylogenetic tree in RAxML with 500 bootstraps using an existing alignment

# Use alignment cleaned in Gblocks for the tree input
NAME=LargeTreeDD04PlacementCleanedHalfGaps

PATHS=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/Phylogenetics/FinalizingTrees/LargeTreeDD04Placement

raxmlHPC=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/Phylogenetics/standard-RAxML/raxmlHPC

OUTGROUP=Streptococcus_suis_BM407

cd $PATHS

$raxmlHPC -f a -o $OUTGROUP -m GTRGAMMA -p 05149 -x 44297 -# 500 -s "${NAME}".fa -n "${NAME}"
