# View S. sinensis phylogenetic trees
library(tidyverse)
library(ggtree)
library(dplyr)
library(treeio)
# https://4va.github.io/biodatasci/r-ggtree.html

dates <- read.delim("~/Desktop/MicroARCH/CompetitiveMappingBarPlots/BritishFormattedDates_3-7-2024.txt", header = TRUE, check.names = FALSE)
tree <- read.tree("~/Desktop/Phylogenetics/RAxML_bipartitions.DD04CMTreeCleanedHalfGaps") # Competitive mapping tree
tree <- read.tree("~/Desktop/Phylogenetics/RAxML_bipartitions.SanguinisOutgroupAdditionalReferencesUnstrictCleanedHalfGaps") # Single reference tree n = 0.01
#tree <- read.tree("~/Desktop/Phylogenetics/RAxML_bipartitions.SanguinisOutgroupAdditionalReferencesCleanedHalfGaps08032024") # Single reference tree n = 0.1

#tree$tip.label <- gsub("_DD04_competitivemappingconsensus", "", tree$tip.label)
#tree$tip.label <- gsub("_Streptococcus_spDD04_Strict_consensus", "", tree$tip.label)
tree$tip.label <- gsub("_Streptococcus_spDD04_consensus", "", tree$tip.label)
tree$tip.label <- gsub("Streptococcus_sp", "Streptococus sinensis DD04", tree$tip.label)
tree$tip.label <- gsub("16830_Medieval_SpitalSquare", "16830_Medieval_SpitalSquare_SSQ88", tree$tip.label)

# Add in dates
tipnames <- as.data.frame(tree$tip.label)
colnames(tipnames) <- "#SampleID"
merged <- merge(tipnames, dates, by = "#SampleID", all.x = TRUE)

# Classify modern references
merged <- merged %>%
  mutate(BlockPeriodSimplified = ifelse(is.na(BlockPeriodSimplified) & grepl("sanguinis", `#SampleID`), "Streptococcus sanguinis reference", BlockPeriodSimplified)) %>%
  mutate(BlockPeriodSimplified = ifelse(is.na(BlockPeriodSimplified) & grepl("cristatus", `#SampleID`), "Streptococcus cristatus reference", BlockPeriodSimplified)) %>%
  mutate(BlockPeriodSimplified = ifelse(is.na(BlockPeriodSimplified) & grepl("sinensis", `#SampleID`), "Streptococcus sinensis reference", BlockPeriodSimplified))

# Add in London classification
LondonClassification <- data.frame(
  SampleID = c("16813_Medieval_SpitalSquare_NRF88", "16818_Medieval_SpitalSquare_NRF88", "16830_Medieval_SpitalSquare_NRF88", "16839_Medieval_MertonPriory_MPY86", "16887_Medieval_MertonPriory_MPY86", "16955_Medieval_GuildhallYard_GYE92"),
  Location = c("Yes", "Yes", "Yes", "Yes", "Yes", "Yes"))
colnames(LondonClassification) <- c("#SampleID", "London")
merged <- merge(merged, LondonClassification, by = "#SampleID", all.x = TRUE)
merged$London <- replace_na(merged$London, "No")

DD04tree <- ggtree(tree) %<+% merged + geom_treescale() +
  geom_tiplab(hjust = -0.05) +
  geom_nodelab(hjust = 1.5, vjust = -1) +
  geom_tree(size = 1.5) +
  geom_tippoint(aes(color = BlockPeriodSimplified, shape = London), size = 5) +
  theme(legend.position = "left") +
  hexpand(1, direction = 1) +
  #geom_text(aes(label=node), hjust = -0.3) +
  scale_color_manual(values = c(
    "Pre-Roman Britian" = "#bfd3e6",
    "Roman Britain" = "#9ebcda",
    "Anglo-Saxon or Early Medieval Britain" = "#8c96c6",
    "Norman Britain and the Middle Ages" = "#8c6bb1",
    "Reformation" = "#88419d",
    "Streptococcus sanguinis reference" = "#35978f",
    "Streptococcus sinensis reference" = "#543005",
    "Streptococcus cristatus reference" = "#bf812d")
  )

ggsave ("~/Desktop/Phylogenetics/DD04CMSanguinisOutgroup.pdf", width = 12, height = 8, DD04tree) # Competitive mapping tree
ggsave ("~/Desktop/Phylogenetics/DD04-SanguinisOutgroup-Unstrict.pdf", width = 12, height = 8, DD04tree) # Single reference tree n = 0.01
#ggsave ("~/Desktop/Phylogenetics/DD04-SanguinisOutgroup-08032024.pdf", width = 12, height = 8, DD04tree) # Single reference tree n = 0.1


equallength <- ggtree(tree, branch.length = "none") %<+% merged + geom_treescale() +
  geom_tiplab(hjust = -0.05) +
  geom_nodelab(hjust = 1.5, vjust = -1) +
  geom_tree(size = 1.1) +
  geom_tippoint(aes(color = BlockPeriodSimplified, shape = London), size = 5) +
  theme(legend.position = "left") +
  hexpand(1, direction = 1) +
  #geom_text(aes(label=node), hjust = -0.3) +
  scale_color_manual(values = c(
    "Pre-Roman Britian" = "#bfd3e6",
    "Roman Britain" = "#9ebcda",
    "Anglo-Saxon or Early Medieval Britain" = "#8c96c6",
    "Norman Britain and the Middle Ages" = "#8c6bb1",
    "Reformation" = "#88419d",
    "Streptococcus sanguinis reference" = "#80cdc1",
    "Streptococcus sinensis reference" = "#35978f",
    "Streptococcus cristatus reference" = "#8c510a")
  )

#ggsave ("~/Desktop/Phylogenetics/DD04CM-SanguinisOutgroup-EqualLength.pdf", width = 12, height = 8, equallength) # Competitive mapping tree
ggsave ("~/Desktop/Phylogenetics/DD04-SanguinisOutgroup-EqualLength-Unstrict.pdf", width = 12, height = 8, equallength) # Single reference tree n = 0.1
#ggsave ("~/Desktop/Phylogenetics/DD04-SanguinisOutgroup-EqualLength-08032024.pdf", width = 12, height = 8, equallength) # Single reference tree n = 0.01

##############################################################################################################
# View large Streptococcus phylogenetic trees (see broader placement of sinensis and ancient sequences)

#placementtree <- read.tree("~/Desktop/Phylogenetics/RAxML_bipartitions.LargeTreeDD04PlacementCleanedHalfGaps") # Single reference tree
placementtree <- read.tree("/Users/avagabrys/Desktop/Phylogenetics/RAxML_bipartitions.LargeTreeCMCleanedHalfGaps") # Competitive mapping

placementtree$tip.label <- gsub("Streptococcus_sp_DD04", "Streptococus_sinensis_DD04", placementtree$tip.label)
#placementtree$tip.label <- gsub("_Streptococcus_spDD04_Strict_consensus", "", placementtree$tip.label)
#placementtree$tip.label <- gsub("_Streptococcus_sanguinis_Strict_consensus", "", placementtree$tip.label)
placementtree$tip.label <- gsub("_DD04_competitivemappingconsensus", "", placementtree$tip.label)
placementtree$tip.label <- gsub("_Sanguinis_competitivemappingconsensus", "", placementtree$tip.label)

# Add in classification of ancient sequences
ancientclassification <- data.frame(
  SampleID = c("StHelensonWalls8875", "Hinxton16489", "BreedonOnTheHill12857", "16955_Medieval_GuildhallYard_GYE92", "Yorkshire8895", "Yorkshire8893", "16857_Medieval_MertonPriory_MPY86", "16847_Medieval_MertonPriory_MPY86", "16918_Medieval_StBenetSherehog_ONE94"),
  Type = c("Ancient Streptococcus sinensis", "Ancient Streptococcus sinensis", "Ancient Streptococcus sinensis", "Ancient Streptococcus sinensis", "Ancient Streptococcus sinensis", "Ancient Streptococcus sanguinis", "Ancient Streptococcus sanguinis", "Ancient Streptococcus sanguinis", "Ancient Streptococcus sanguinis"))
colnames(ancientclassification) <- c("#SampleID", "Classification")

tipnames <- as.data.frame(placementtree$tip.label)
colnames(tipnames) <- "#SampleID"
mergedforlargetree <- merge(tipnames, ancientclassification, by = "#SampleID", all.x = TRUE)

largetree <- ggtree(placementtree) %<+% mergedforlargetree +
  geom_treescale(y = -1) +
  geom_tiplab(hjust = -0.1) +
  geom_nodelab(hjust = 1.3, vjust = -0.5) +
  geom_tree(size = 1.1) +
  geom_tippoint(aes(color = Classification), size = 4) +
  scale_color_manual(values = c("Ancient Streptococcus sinensis" = "#543005",
                                 "Ancient Streptococcus sanguinis" = "#35978f")) +
    #theme(legend.position = "left") +
    #hexpand(1, direction = 1) +
    #geom_text(aes(label=node), hjust = -0.3) +
    #scale_color_manual(values = c( "Ancient Streptococcus sinensis" = "#bfd3e6",
    #"Ancient Streptococcus sanguinis" = "#9ebcda"))
  hexpand(0.2, direction = 1)

#ggsave ("~/Desktop/Phylogenetics/LargeStreptococcusTree.pdf", width = 20, height = 8, largetree)
ggsave ("~/Desktop/Phylogenetics/LargeStreptococcusTreeCM.pdf", width = 20, height = 8, largetree)
