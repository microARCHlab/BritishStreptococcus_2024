library(ggplot2)
library(tidyr)
library(svglite)

# Read in competitive mapping read counts and dates
ancientCMrawdata <- read.delim("~/Desktop/MicroARCH/CompetitiveMappingBarPlots/Ancient-AllSpeciesCollectedCompetitiveMappingResults.txt", header = TRUE, check.names = FALSE)
dates <- read.delim("~/Desktop/MicroARCH/CompetitiveMappingBarPlots/BritishFormattedDates_3-7-2024.txt", header = TRUE, check.names = FALSE)

ancientCMrawdata <- ancientCMrawdata[, -c(2:6)]
ancientCMrawdata$Sample <- gsub(" ", "", ancientCMrawdata$Sample)

# Only include samples with >1000 reads for at least one strep species
ancientCM_countfiltered <- ancientCMrawdata[rowSums(ancientCMrawdata[, -(1)] > 1000) >= 1, ]

# Include only the 12 most abundant species on average
speciesaverages <- colMeans(ancientCM_countfiltered[, -(1)])
speciesaveragessorted <- names(ancientCM_countfiltered[, -(1)][order(speciesaverages, decreasing = TRUE)])
ancientCM_countfiltered[, -(1)] <-  ancientCM_countfiltered[, speciesaveragessorted]
colnames(ancientCM_countfiltered) <- c("Sample", speciesaveragessorted)
ancientCM_countspeciesfiltered <- ancientCM_countfiltered[, 1:11]

# Sort samples by early date
dates$AdjustedEarlyDate <- as.numeric(dates$AdjustedEarlyDate)
dates <- dates[order(dates$`AdjustedEarlyDate`, na.last = TRUE, decreasing = FALSE), ]

# Remove those which have no early date (where NA)
dates <- dates[complete.cases(dates$AdjustedEarlyDate), ]
orderedsamplelist <- dates[, 1]
orderedsamplelist <- orderedsamplelist[orderedsamplelist %in% ancientCM_countspeciesfiltered[, 1]]

rownames(ancientCM_countspeciesfiltered) <- ancientCM_countspeciesfiltered[, 1]
ancientCM_sortedfiltered <- ancientCM_countspeciesfiltered[orderedsamplelist, ] 

# Arrange for plotting and plot

ancientCM_long <- ancientCM_sortedfiltered %>% pivot_longer(cols = -Sample, names_to = "Species", values_to = "ReadCount")
ancientCMplot <- ggplot(ancientCM_long, aes(x = factor(Sample, orderedsamplelist), y = ReadCount, fill = Species)) + 
  geom_bar(position="stack", stat="identity", width = 1) +
  #geom_bar(position="fill", stat="identity", width = 1) +
  labs(title = "Ancient Competitive Mapping", x = "Sample", y = "Relative Abundance") +
  guides(fill = guide_legend(reverse = TRUE)) +
  xlab("") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, color = "black", size = 8), 
        axis.text.y = element_text(color = "black", size = 10)) +
  scale_fill_manual(values = c("#000000", "#543005","#8c510a","#bf812d","#dfc27d","#f6e8c3","#f5f5f5","#c7eae5","#80cdc1","#35978f","#01665e","#003c30"))

ggsave("ancientCM-stack.svg", "~/Desktop/MicroARCH/CompetitiveMappingBarPlots", plot = ancientCMplot, device = "svg", width = 13, height = 8)
#ggsave("ancientCM.svg", "~/Desktop/MicroARCH/CompetitiveMappingBarPlots", plot = ancientCMplot, device = "svg", width = 13, height = 8)

## Plot Mutans group (was not included in top 12 most abundant species) ##
# Plot the same samples (ones with >1000 reads for at least one strep species and specific date)
ancientmutansgrouponly <- ancientCM_countfiltered[, c("Sample", "S.mutans", "S.downei", "S.troglodytae", "S.sobrinus")]
# Sort samples by early date
orderedsamplelistancientmutansgroup <- orderedsamplelist[orderedsamplelist %in% ancientmutansgrouponly[, 1]]
rownames(ancientmutansgrouponly) <- ancientmutansgrouponly[, 1]
ancientmutansgrouponly <- ancientmutansgrouponly[orderedsamplelistancientmutansgroup, ] 

# Arrange for plotting and plot
mutansancientCM_long <- ancientmutansgrouponly %>% pivot_longer(cols = -Sample, names_to = "Species", values_to = "ReadCount")
mutansancientCMplot <- ggplot(mutansancientCM_long, aes(x = factor(Sample, orderedsamplelistancientmutansgroup), y = ReadCount, fill = Species)) + 
  geom_bar(position="stack", stat="identity", width = 1) +
  #geom_bar(position="fill", stat="identity", width = 1) +
  labs(title = "Ancient Competitive Mapping", x = "Sample", y = "Relative Abundance") +
  guides(fill = guide_legend(reverse = TRUE)) +
  xlab("") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, color = "black", size = 8), 
        axis.text.y = element_text(color = "black", size = 10)) +
  scale_fill_manual(values = c("#000000", "#543005","#8c510a","#bf812d","#dfc27d","#f6e8c3","#f5f5f5","#c7eae5","#80cdc1","#35978f","#01665e","#003c30"))
ggsave("~/Desktop/MicroARCH/CompetitiveMappingBarPlots/mutansancient.pdf", plot = mutansancientCMplot, width = 13, height = 8)

###########################################################################
healthymodernCMrawdata <- read.delim("~/Desktop/MicroARCH/CompetitiveMappingBarPlots/HealthyModern-AllSpeciesCollectedCompetitiveMappingResults.txt", header = TRUE, check.names = FALSE)
periomodernCMrawdata <- read.delim("~/Desktop/MicroARCH/CompetitiveMappingBarPlots/PerioModern-AllSpeciesCollectedCompetitiveMappingResults.txt", header = TRUE, check.names = FALSE)

combinedmodernCMrawdata <- rbind(healthymodernCMrawdata, periomodernCMrawdata)

combinedmodernCMrawdata <- combinedmodernCMrawdata[, -c(2:6)]
combinedmodernCMrawdata$Sample <- gsub(" ", "", combinedmodernCMrawdata$Sample)

# Only include samples with >1000 reads for at least one strep species
combinedmodernCM_countfiltered <- combinedmodernCMrawdata[rowSums(combinedmodernCMrawdata[, -(1)] > 1000) >= 1, ]

# Include only the 12 most abundant species on average
speciesaverages <- colMeans(combinedmodernCM_countfiltered[, -(1)])
speciesaveragessorted <- names(combinedmodernCM_countfiltered[, -(1)][order(speciesaverages, decreasing = TRUE)])
combinedmodernCM_countfiltered[, -(1)] <-  combinedmodernCM_countfiltered[, speciesaveragessorted]
colnames(combinedmodernCM_countfiltered) <- c("Sample", speciesaveragessorted)
combinedmodernCM_countspeciesfiltered <- combinedmodernCM_countfiltered[, 1:11]
combinedmodernCM_long <- combinedmodernCM_countspeciesfiltered %>% pivot_longer(cols = -Sample, names_to = "Species", values_to = "ReadCount")

combinedmodernCMplot <- ggplot(combinedmodernCM_long, aes(x = Sample, y = ReadCount, fill = Species)) + 
  geom_bar(position="stack", stat="identity", width = 1) +
  #geom_bar(position="fill", stat="identity", width = 1) +
  labs(title = "Modern Competitive Mapping", x = "Sample", y = "Relative Abundance") +
  guides(fill = guide_legend(reverse = TRUE)) +
  xlab("") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, color = "black", size = 8), 
        axis.text.y = element_text(color = "black", size = 10)) +
  scale_fill_manual(values = c("#000000", "#543005","#8c510a","#bf812d","#dfc27d","#f6e8c3","#f5f5f5","#c7eae5","#80cdc1","#35978f","#01665e","#003c30"))

#ggsave("combinedmodernCM.svg", "~/Desktop/MicroARCH/CompetitiveMappingBarPlots", plot = combinedmodernCMplot, device = "svg", width = 7.5, height = 8)
ggsave("combinedmodernCM-stackedTEST.svg", "~/Desktop/MicroARCH/CompetitiveMappingBarPlots", plot = combinedmodernCMplot, device = "svg", width = 7.5, height = 8)

## Plot Mutans group (was not included in top 12 most abundant species) ##
# Plot the same samples (ones with >1000 reads for at least one strep species and specific date)
modernmutansgrouponly <- combinedmodernCM_countfiltered[, c("Sample", "S.mutans", "S.downei", "S.troglodytae", "S.sobrinus")]

# Arrange for plotting and plot
mutansmodernCM_long <- modernmutansgrouponly %>% pivot_longer(cols = -Sample, names_to = "Species", values_to = "ReadCount")
muttansmodernCMplot <- ggplot(mutansmodernCM_long, aes(x = Sample, y = ReadCount, fill = Species)) + 
  geom_bar(position="stack", stat="identity", width = 1) +
  #geom_bar(position="fill", stat="identity", width = 1) +
  labs(title = "Modern Competitive Mapping", x = "Sample", y = "Relative Abundance") +
  guides(fill = guide_legend(reverse = TRUE)) +
  xlab("") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank()) +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, color = "black", size = 8), 
        axis.text.y = element_text(color = "black", size = 10)) +
  scale_fill_manual(values = c("#000000", "#543005","#8c510a","#bf812d","#dfc27d","#f6e8c3","#f5f5f5","#c7eae5","#80cdc1","#35978f","#01665e","#003c30"))

ggsave("~/Desktop/MicroARCH/CompetitiveMappingBarPlots/mutansmodern.pdf", plot = muttansmodernCMplot, width = 13, height = 8)


