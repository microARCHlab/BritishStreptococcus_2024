library(dplyr)
library(tidyverse)

# Read in data
# Coverage summary
rawdata <- read.delim("~/Desktop/MicroARCH/AncientStrepPaper/GenomicIslands/Streptococcus_spDD04_GenomicIslandCoverageSummary.txt", header = TRUE, check.names = FALSE)
#rawdata <- read.delim("~/Desktop/MicroARCH/AncientStrepPaper/GenomicIslands/Streptococcus_sanguinis_GenomicIslandCoverageSummary_Combined.txt", header = TRUE, check.names = FALSE)
# Island attributes
rawislanddata <- read.delim("~/Desktop/MicroARCH/AncientStrepPaper/GenomicIslands/Streptococcus_spDD04_GenomicIslandAttributes.txt", header = TRUE, check.names = FALSE)
#rawislanddata <- read.delim("~/Desktop/MicroARCH/AncientStrepPaper/GenomicIslands/Streptococcus_sanguinis_GenomicIslandAttributes.txt", header = TRUE, check.names = FALSE)
# Sample attributes
rawsampledata <- read.delim("~/Desktop/MicroARCH/AncientStrepPaper/GenomicIslands/Streptococcus_spDD04_SampleAttributes.txt", header = TRUE, check.names = FALSE)
#rawsampledata <- read.delim("~/Desktop/MicroARCH/AncientStrepPaper/GenomicIslands/Streptococcus_sanguinis_SampleAttributes_Combined.txt", header = TRUE, check.names = FALSE)
# Island characterization (made afer filtered islands)
islandcharacterization <- read.delim("~/Desktop/MicroARCH/AncientStrepPaper/GenomicIslands/DD04_Characterization.txt", header = TRUE, check.names = FALSE)
#islandcharacterization <- read.delim("~/Desktop/MicroARCH/AncientStrepPaper/GenomicIslands/Sanguinis_Characterization.txt", header = TRUE, check.names = FALSE)

rownames <- rawdata[, 1]
islandcovs <- rawdata[, -(1)]
rownames(islandcovs) <- rownames
islandcovs <- islandcovs %>%
  mutate_all(~ . * 100)

################################
## Islands ##
islandrownames <- gsub(" ", "", (rawislanddata[, 1]))
islandmetadata <- rawislanddata[, -(1)]
rownames(islandmetadata) <- islandrownames

# Filter to only include islands of specified length
minislandlength <- 2000
minbreadthofcoverage <- 0.60
islands <- islandmetadata %>% filter(`Island Length` > minislandlength)
selected_data <- islandcovs[rownames(islandcovs) %in% rownames(islands), ]
selected_data <- selected_data %>% arrange(rownames(selected_data))

################################
## Samples ##
samplerownames <- gsub(" ", "", (rawsampledata[, 1]))
samplemetadata <- rawsampledata[, -(1)]
rownames(samplemetadata) <- samplerownames

# Filter samples by overall genome coverage
filteredbybreadth <- samplemetadata %>%
  filter(`Whole Genome Coverage` > minbreadthofcoverage)
filteredbybreadth$`Early Date` <- as.numeric(filteredbybreadth$`Early Date`)
filteredbybreadth$`Whole Genome Coverage` <- as.numeric(filteredbybreadth$`Whole Genome Coverage`)
filteredbybreadth$`Relative Abundance` <- as.numeric(filteredbybreadth$`Relative Abundance`)
filteredbybreadth$`Whole Genome Coverage` <- filteredbybreadth$`Whole Genome Coverage` * 100
filteredbybreadth$`Relative Abundance` <- filteredbybreadth$`Relative Abundance` * 100
filteredbybreadthsorted <- filteredbybreadth[order(filteredbybreadth$`Early Date`),]

samplelist <- rownames(filteredbybreadthsorted)
selected_data <- selected_data %>% select(any_of(samplelist))
selected_data <- as.data.frame(t(selected_data))
percent <- selected_data

# Make a version of genomic island data that is presence or absence of the island, not percent covered
preabs <- as.data.frame(apply(selected_data, 2, function(x) ifelse(x > 15, 1, 0)))

# Join genomic island data with metadata
filteredbybreadthsorted$Sample <- rownames(filteredbybreadthsorted)
percent$Sample <- rownames(percent)
preabs$Sample <- rownames(preabs)
percent <- left_join(percent, filteredbybreadthsorted)
preabs <- left_join(preabs, filteredbybreadthsorted)

# Make date a different format
preabs$DateCatNum <- preabs$Date
preabs$DateCatNum[preabs$DateCatNum == " Pre-Roman Britian"] <- 1
preabs$DateCatNum[preabs$DateCatNum == " Roman Britain"] <- 2
preabs$DateCatNum[preabs$DateCatNum == " Anglo-Saxon or Early Medieval Britain"] <- 3
preabs$DateCatNum[preabs$DateCatNum == " Norman Britain and the Middle Ages"] <- 4
preabs$DateCatNum[preabs$DateCatNum == " Reformation"] <- 5

################################
## Statistics ##

# Use glm model to test over time
models <- list() 
row_names <- colnames(preabs[,-(37:46)])
col_names <- c("Intercept Estimate", "Intercept Std. Error", "Intercept z value", "Intercept Pr(>|z|)",
               "Early Date Estimate", "Early Date Std. Error", "Early Date z value", "Early Date Pr(>|z|)", "Whole Genome Coverage Estimate", "Whole Genome Coverage Std. Error", "Whole Genome Coverage z value", "Whole Genome Coverage Pr(>|z|)")
results <- data.frame(matrix(ncol = length(col_names), nrow = length(row_names)))
rownames(results) <- row_names
colnames(results) <- col_names

for (Island in colnames(preabs[,-(37:46)])){
  model <- glm(preabs[[Island]] ~ `Early Date` + `Whole Genome Coverage`, data = preabs, family = binomial)
  model_summary <- summary(model)
  coef_matrix <- model_summary$coefficients
  results[Island, 1:4] <- model_summary$coefficients[1, ]
  results[Island, 5:8] <- model_summary$coefficients[2, ]
  results[Island, 9:12] <- model_summary$coefficients[3, ] }

# Save
write.table(results, file = "~/Desktop/MicroARCH/AncientStrepPaper/GenomicIslands/GIPresenceOverTime.txt",        # Specify your desired file path
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
################################
# Figure out how many islands were completely absent across ancient strains
colsums <- sort(colSums(preabs[,-(37:46)]))
colsums / length(rownames(preabs))
