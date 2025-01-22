## Usage: Create a whole-genome alignment in Mugsy and convert the alignment to .fa using HarvestTools. Resulting FASTA alignment can be cleaned and then used in RAxML.


MUGSY=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/Phylogenetics/mugsy_x86-64-v1r2.3/mugsyenv.sh 

HARVESTTOOLS=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/Phylogenetics/harvesttools-Linux64-v1.2/harvesttools

NAME=LargeTreeDD04Placement

PATHS=/storage/group/LiberalArts/default/lsw132_collab/04_LabMembersFolders/AvaGabrys/Phylogenetics/FinalizingTrees/LargeTreeDD04Placement

cd $PATHS

##################################

SEQUENCES=`cat sequencestoalign"${NAME}".txt`

source $MUGSY

mugsy --directory $PATHS --prefix $NAME $SEQUENCES

##################################

$HARVESTTOOLS -a "${NAME}".maf -M "${NAME}".fa
