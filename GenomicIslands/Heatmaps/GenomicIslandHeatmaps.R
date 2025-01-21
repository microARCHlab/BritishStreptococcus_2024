library(colorRamp2)
library(dplyr)
library(tidyverse)
library(ComplexHeatmap)

minislandlength <- 2000
minbreadthofcoverage <- 0.60

# Read in data
# Coverage summary
#rawdata <- read.delim("~/Desktop/MicroARCH/AncientStrepPaper/GenomicIslands/Streptococcus_spDD04_GenomicIslandCoverageSummary.txt", header = TRUE, check.names = FALSE)
rawdata <- read.delim("~/Desktop/MicroARCH/AncientStrepPaper/GenomicIslands/Streptococcus_sanguinis_GenomicIslandCoverageSummary_Combined.txt", header = TRUE, check.names = FALSE)
# Island attributes
#rawislanddata <- read.delim("~/Desktop/MicroARCH/AncientStrepPaper/GenomicIslands/Streptococcus_spDD04_GenomicIslandAttributes.txt", header = TRUE, check.names = FALSE)
rawislanddata <- read.delim("~/Desktop/MicroARCH/AncientStrepPaper/GenomicIslands/Streptococcus_sanguinis_GenomicIslandAttributes.txt", header = TRUE, check.names = FALSE)
# Sample attributes
#rawsampledata <- read.delim("~/Desktop/MicroARCH/AncientStrepPaper/GenomicIslands/Streptococcus_spDD04_SampleAttributes.txt", header = TRUE, check.names = FALSE)
rawsampledata <- read.delim("~/Desktop/MicroARCH/AncientStrepPaper/GenomicIslands/Streptococcus_sanguinis_SampleAttributes_Combined.txt", header = TRUE, check.names = FALSE)
# Island characterization (made afer filtered islands)
#islandcharacterization <- read.delim("~/Desktop/MicroARCH/AncientStrepPaper/GenomicIslands/DD04_Characterization.txt", header = TRUE, check.names = FALSE)
islandcharacterization <- read.delim("~/Desktop/MicroARCH/AncientStrepPaper/GenomicIslands/Sanguinis_Characterization.txt", header = TRUE, check.names = FALSE)

rownames <- rawdata[, 1]
islandcovs <- rawdata[, -(1)]
rownames(islandcovs) <- rownames
islandcovs <- islandcovs %>%
  mutate_all(~ . * 100)

################################
## Islands ##
islandrownames <- gsub(" ", "", (rawislanddata[, 1]))
islandmetadata <- rawislanddata[, -(1)]
##rownames(islandmetadata) <- gsub(" ", "", gsub("GI-1-", "GI-1", gsub("GI-2-", "GI-2", gsub("GI-3-", "GI-3", gsub("GI-4-", "GI-4", islandrownames)))))
rownames(islandmetadata) <- islandrownames

# Filter to only include islands of specified length
islands <- islandmetadata %>% filter(`Island Length` > minislandlength)
islandlist <- rownames(islands)
selected_data <- islandcovs[rownames(islandcovs) %in% islandlist, ]

selected_data <- selected_data %>% arrange(rownames(selected_data))

# Get annotation data formatted
selectedislandmetadata <- islandmetadata[rownames(islandmetadata) %in% islandlist, ] 
selectedislandmetadata <- selectedislandmetadata %>% arrange(rownames(selectedislandmetadata))
selectedislandmetadata$`GC content` <- selectedislandmetadata$`GC content` * 100
annotationrow <- selectedislandmetadata[, c(1,3)]

# Add in island characterizations and make these the new rownames for plotting (all data frames are sorted the same)
islandcharacterization$`Genomic Island` <- gsub(" ", "", islandcharacterization$`Genomic Island`)
selectedislandcharacterization <- islandcharacterization[islandcharacterization$`Genomic Island` %in% islandlist, ]
selectedislandcharacterization <- selectedislandcharacterization %>% arrange(`Genomic Island`)
characterizationnames <- selectedislandcharacterization$`Summarized for Chart`
rownames(selected_data) <- characterizationnames
rownames(annotationrow) <- characterizationnames

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

matrix <- as.matrix(selected_data)
annotationcol <- filteredbybreadthsorted[ , c(1,2,6)]

################################
## Plot Heatmap ##
habottom = HeatmapAnnotation(
  df = annotationcol[, c("Whole Genome Coverage", "Relative Abundance")],
  simple_anno_size = unit(2.5, "mm"),
  col = list(
    #`Whole Genome Coverage` = colorRamp2(c(0, 50, 100), c("#deebf7", "#6baed6", "#08519c")),
    `Whole Genome Coverage` = colorRamp2(c(0, 10, 100), c("#08519c", "white", "#b30000")),
    `Relative Abundance` = colorRamp2(c(15,100), c("#c994c7", "#980043"))),
  annotation_label = list(
    `Whole Genome Coverage` = "Genome Coverage (%)",
    `Relative Abundance` = "Relative Abundance (%)"),
  annotation_legend_param = list(
    `Whole Genome Coverage` = list(at = c(0, 20, 40, 60, 80, 100))#,
    #`Relative Abundance` = list(at = c(0, 20, 40, 60, 80, 100))
  )
)

haleft = rowAnnotation(
  df = annotationrow,
  col = list(
    `GC content` = colorRamp2(c(30, 50), c("#cccccc","#252525")),
    `Island Length` = colorRamp2(c(2000, 45000), c("#d1ac90", "#4a2b14"))
  ),
  annotation_label = list(
    `GC content` = "Average GC (%)",
    `Island Length` = "Island Length (bp)"
  ),
  simple_anno_size = unit(2.5, "mm")
)

hatop = HeatmapAnnotation(
  df = annotationcol$Date,
  annotation_label = "Period",
  #col = list(
  # " Pre-Roman Britian" = "blue",
  # " Roman Britain" = "red",
  # " Anglo-Saxon or Early Medieval Britain" = "yellow",
  # " Norman Britain and the Middle Ages" = "green",
  # " Reformation" = "pink",
  # " U" = "orange"),
  annotation_legend_param = list(
    labels = c(" Pre-Roman Britian (-43 CE)",
               " Roman Britain (43-410 CE)",
               " Anglo-Saxon (410-1066 CE)",
               " Middle Ages (1066-1547)",
               " Reformation (1547-1750 CE)",
               " Healthy Modern",
               " Perio Modern",
               " U"
    ),
    at = c(" Pre-Roman Britian",
           " Roman Britain",
           " Anglo-Saxon or Early Medieval Britain",
           " Norman Britain and the Middle Ages",
           " Reformation",
           " Healthy Modern",
           " Perio Modern",
           " U")
  ))

Heatmap <- ComplexHeatmap::Heatmap(
  matrix, 
  name = "Coverage (%)", 
  #col = colorRamp2(c(0, 50, 100), c("#deebf7", "#6baed6", "#08519c")),
  #col = colorRamp2(c(0, 50, 100), c("#08519c", "white", "#b30000")),
  col = colorRamp2(c(0, 10, 100), c("#08519c", "white", "#b30000")),
  column_split = rep(c("Ancient", "HealthyModern", "PerioModern"), c(16, 11, 7)),
  column_title = "Streptococcus sanguinis Genomic Islands",
  column_order = order(as.numeric(colnames(matrix))),
  bottom_annotation = (habottom),
  left_annotation = (haleft),
  top_annotation = (hatop),
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  heatmap_legend_param = list(
    at = c(0, 20, 40, 60, 80, 100),
    legend_gp = gpar(fontsize = 12)),
  width = 10,
  height = 10
)

#####################################################################
# Additional steps to add in values for the heatmap for easier referencing with the supplemental tables
# Get row order, and sort by this
ht <- draw(Heatmap)
roworder <- row_order(ht)
selectedislandcharacterization_rowlabels <- selectedislandcharacterization[roworder, ]

# Add in row number values in the sorted frame; export for later reference. Also get rid of "-1" that was temporarily added to avoid duplicate names
selectedislandcharacterization_rowlabels$`Summarized for Chart` <- paste0(seq_len(nrow(selectedislandcharacterization_rowlabels)), ". ", selectedislandcharacterization_rowlabels$`Summarized for Chart`)
selectedislandcharacterization_rowlabels$`Summarized for Chart` <- gsub("-1", "", selectedislandcharacterization_rowlabels$`Summarized for Chart`)
write_tsv(selectedislandcharacterization_rowlabels, "~/Desktop/MicroARCH/AncientStrepPaper/GenomicIslands/Sanguinis-sortedcharacterizations.tsv")

# Also add in fraction of genes with GO terms
selectedislandcharacterization_rowlabels$`Summarized for Chart` <- paste0(selectedislandcharacterization_rowlabels$`Summarized for Chart`, selectedislandcharacterization_rowlabels$`Fraction of CDS Genes with GO Terms`)
#view(selectedislandcharacterization_rowlabels)

################################################################
# Make the heatmap again!
# Make sure input order matches again
selectedislandcharacterization_rowlabels <- selectedislandcharacterization_rowlabels %>% arrange(`Genomic Island`)
characterizationnames <- selectedislandcharacterization_rowlabels$`Summarized for Chart`
rownames(selected_data) <- characterizationnames
rownames(annotationrow) <- characterizationnames

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

matrix <- as.matrix(selected_data)
annotationcol <- filteredbybreadthsorted[ , c(1,2,6)]

################################
## Plot Heatmap ##
habottom = HeatmapAnnotation(
  df = annotationcol[, c("Whole Genome Coverage", "Relative Abundance")],
  simple_anno_size = unit(2.5, "mm"),
  col = list(
    #`Whole Genome Coverage` = colorRamp2(c(0, 50, 100), c("#deebf7", "#6baed6", "#08519c")),
    `Whole Genome Coverage` = colorRamp2(c(0, 10, 100), c("#08519c", "white", "#b30000")),
    `Relative Abundance` = colorRamp2(c(15,100), c("#c994c7", "#980043"))),
  annotation_label = list(
    `Whole Genome Coverage` = "Genome Coverage (%)",
    `Relative Abundance` = "Relative Abundance (%)"),
  annotation_legend_param = list(
    `Whole Genome Coverage` = list(at = c(0, 20, 40, 60, 80, 100))#,
    #`Relative Abundance` = list(at = c(0, 20, 40, 60, 80, 100))
  )
)

haleft = rowAnnotation(
  df = annotationrow,
  col = list(
    `GC content` = colorRamp2(c(30, 50), c("#cccccc","#252525")),
    `Island Length` = colorRamp2(c(2000, 45000), c("#d1ac90", "#4a2b14"))
  ),
  annotation_label = list(
    `GC content` = "Average GC (%)",
    `Island Length` = "Island Length (bp)"
  ),
  simple_anno_size = unit(2.5, "mm")
)

hatop = HeatmapAnnotation(
  df = annotationcol$Date,
  annotation_label = "Period",
  #col = list(
  # " Pre-Roman Britian" = "blue",
  # " Roman Britain" = "red",
  # " Anglo-Saxon or Early Medieval Britain" = "yellow",
  # " Norman Britain and the Middle Ages" = "green",
  # " Reformation" = "pink",
  # " U" = "orange"),
  annotation_legend_param = list(
    labels = c(" Pre-Roman Britian (-43 CE)",
               " Roman Britain (43-410 CE)",
               " Anglo-Saxon (410-1066 CE)",
               " Middle Ages (1066-1547)",
               " Reformation (1547-1750 CE)",
               " Healthy Modern",
               " Perio Modern",
               " U"
    ),
    at = c(" Pre-Roman Britian",
           " Roman Britain",
           " Anglo-Saxon or Early Medieval Britain",
           " Norman Britain and the Middle Ages",
           " Reformation",
           " Healthy Modern",
           " Perio Modern",
           " U")
  ))

Heatmap <- ComplexHeatmap::Heatmap(
  matrix, 
  name = "Coverage (%)", 
  #col = colorRamp2(c(0, 50, 100), c("#deebf7", "#6baed6", "#08519c")),
  #col = colorRamp2(c(0, 50, 100), c("#08519c", "white", "#b30000")),
  col = colorRamp2(c(0, 10, 100), c("#08519c", "white", "#b30000")),
  column_split = rep(c("Ancient", "HealthyModern", "PerioModern"), c(16, 11, 7)),
  column_title = "Streptococcus sanguinis Genomic Islands",
  column_order = order(as.numeric(colnames(matrix))),
  bottom_annotation = (habottom),
  left_annotation = (haleft),
  top_annotation = (hatop),
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8),
  heatmap_legend_param = list(
    at = c(0, 20, 40, 60, 80, 100),
    legend_gp = gpar(fontsize = 12)),
  width = 10,
  height = 10
)

pdf("~/Desktop/MicroARCH/AncientStrepPaper/GenomicIslands/SanguinisHM_01212024.pdf", width = 11, height = 8)
#pdf("~/Desktop/MicroARCH/AncientStrepPaper/GenomicIslands/SinensisHM_01212024.pdf", width = 11, height = 8)
draw(Heatmap)
dev.off()

