# Strep Paper | Competitive Mapping Associations
#
# Abigail S. Gancz, 02/06/2024
# Edited 03/28/2024 with Ava's age estimations 
# Edited 08/05/2024 with additional variables
# Edited 09/02/2024 with updated filtering and only deltaD for DD04
# Edited 10/02/2024 to filter for samples that have at least >1000 reads for a single species 
#
# Usage: Run MAASLIN2 on Streptoccocus competitive mapping results to determine if associations exist
#
############################################################################################################################

# Package Installation

############################################################################################################################

#function to help with package calling/installation | DO NOT CHANGE
using<-function(...) {
  libs<-unlist(list(...))
  req<-unlist(lapply(libs,require,character.only=TRUE))
  need<-libs[req==FALSE]
  if(length(need)>0){ 
    install.packages(need)
    lapply(need,require,character.only=TRUE)
  }
}

#if(!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("Maaslin2")

#install.packages("proxy")

#call desired packages, and install if not available 
using("ggplot2", "dplyr", "vegan", "RColorBrewer", "data.table", "Maaslin2", "reshape2", "tidyverse", "magrittr", "readxl")


##############################################################################################################################################

# Set working directory

############################################################################################################################

setwd("C:/Users/aganc/Desktop/Professional/Projects/Strep_Paper/analyses/r_coding")

############################################################################################################################

# Input data

############################################################################################################################

#insert data
data <- read.delim("C:/Users/aganc/Desktop/Professional/Projects/PSU_GraduateSchool_Projects/Strep_Paper/data/AllSpeciesCollectedCompetitiveMappingResults_08052024.txt")

#filter data
data <- data[, c(1,7:30)] #remove mapping info columns, keeping only different species

#insert metadata 
metadata <- read_excel("C:/Users/aganc/Desktop/Professional/Projects/PSU_GraduateSchool_Projects/Strep_Paper/data/Metadata_AllSpeciesCollectedCompetitiveMappingResults_09022024.xlsx")
metadata_filtered <- filter(metadata, metadata$StrepFiltered == "Y")
metadata_filtered <- filter(metadata, metadata$reads_1000 == "Y")

#filter metadata
metadata_filtered <- merge(x=metadata_filtered, y=data, by="Sample") #ensure metadata matches data contents
MoL_metadata_filtered <- filter(metadata_filtered, metadata_filtered$Museum =="MoL") #create metadata subset featuring only MoL samples

#add rownmaes to metadata
rownames(metadata_filtered ) <- metadata_filtered [, 1]
rownames(MoL_metadata_filtered ) <- MoL_metadata_filtered [, 1]

#add rownames to data
rownames(data) <- data[, 1]
data <- data[, -1]

############################################################################################################################

# Run Maaslin2

############################################################################################################################

###########
#Oral Health Variables
###########

#oral_Pathology #NOT ENOUGTH SAMPLES
oral_Pathology_metadata_filtered <- filter(MoL_metadata_filtered, MoL_metadata_filtered$Oral_Pathology == "1" | MoL_metadata_filtered$Oral_Pathology == "0")
table(oral_Pathology_metadata_filtered$Oral_Pathology) 


#Abscess_YN
Abscess_YN_metadata_filtered <- filter(MoL_metadata_filtered, MoL_metadata_filtered$Abscess_YN == "1" | MoL_metadata_filtered$Abscess_YN == "0")
table(Abscess_YN_metadata_filtered$Abscess_YN) 

Abscess_YN_metadata_filtered$Tooth_modified = factor(Abscess_YN_metadata_filtered$Tooth,
                                                     levels = c("Molar", "Bicuspid", "Canine", "Incisor"))

fit_data <- Maaslin2(
  input_data = data,
  input_metadata = Abscess_YN_metadata_filtered,
  output = "C:/Users/aganc/Desktop/Professional/Projects/PSU_GraduateSchool_Projects/Strep_Paper/analyses/r_coding/MAASLIN2_Results_10072024/Abscess_YN",
  fixed_effects = c('Abscess_YN'),
  random_effects = c('Tooth_modified'),
  normalization = "CLR", 
  transform = "LOG",
  max_significance = .05,
  min_prevalence = 0.1,
  min_abundance = 0.0001)


#Caries
Caries_metadata_filtered <- filter(MoL_metadata_filtered, MoL_metadata_filtered$Caries == "1" | MoL_metadata_filtered$Caries == "0")
table(Caries_metadata_filtered$Caries) 

Caries_metadata_filtered$Tooth_modified = factor(Caries_metadata_filtered$Tooth,
                                                 levels = c("Molar", "Bicuspid", "Canine", "Incisor"))

fit_data <- Maaslin2(
  input_data = data,
  input_metadata = Caries_metadata_filtered,
  output = "C:/Users/aganc/Desktop/Professional/Projects/PSU_GraduateSchool_Projects/Strep_Paper/analyses/r_coding/MAASLIN2_Results_10072024/Caries",
  fixed_effects = c('Caries'),
  random_effects = c('Tooth_modified'),
  normalization = "CLR", 
  transform = "LOG",
  max_significance = .05,
  min_prevalence = 0.1,
  min_abundance = 0.0001)


#Periodontitis_YN
Periodontitis_YN_metadata_filtered <- filter(MoL_metadata_filtered, MoL_metadata_filtered$Periodontitis_YN == "1" | MoL_metadata_filtered$Periodontitis_YN == "0")
table(Periodontitis_YN_metadata_filtered$Periodontitis_YN) 

Periodontitis_YN_metadata_filtered$Tooth_modified = factor(Periodontitis_YN_metadata_filtered$Tooth,
                                                           levels = c("Molar", "Bicuspid", "Canine", "Incisor"))

fit_data <- Maaslin2(
  input_data = data,
  input_metadata = Periodontitis_YN_metadata_filtered,
  output = "C:/Users/aganc/Desktop/Professional/Projects/PSU_GraduateSchool_Projects/Strep_Paper/analyses/r_coding/MAASLIN2_Results_10072024/Periodontitis_YN",
  fixed_effects = c('Periodontitis_YN'),
  random_effects = c('Tooth_modified'),
  normalization = "CLR", 
  transform = "LOG",
  max_significance = .05,
  min_prevalence = 0.1,
  min_abundance = 0.0001)


###########
#Systemic Health Variables
###########

#Overall_Pathology
Overall_Pathology_metadata_filtered <- filter(MoL_metadata_filtered, MoL_metadata_filtered$Overall_Pathology== "1" | MoL_metadata_filtered$Overall_Pathology == "0")

table(Overall_Pathology_metadata_filtered$Overall_Pathology) 

Overall_Pathology_metadata_filtered$Overall_Pathology <- as.factor(Overall_Pathology_metadata_filtered$Overall_Pathology)

Overall_Pathology_metadata_filtered$Tooth_modified = factor(Overall_Pathology_metadata_filtered$Tooth,
                                                            levels = c("Molar", "Bicuspid", "Canine", "Incisor"))

fit_data <- Maaslin2(
  input_data = data,
  input_metadata = Overall_Pathology_metadata_filtered,
  output = "C:/Users/aganc/Desktop/Professional/Projects/PSU_GraduateSchool_Projects/Strep_Paper/analyses/r_coding/MAASLIN2_Results_10072024/Overall_Pathology",
  fixed_effects = c('Overall_Pathology'),
  random_effects = c('Tooth_modified'),
  normalization = "CLR", 
  transform = "LOG",
  max_significance = .05,
  min_prevalence = 0.1,
  min_abundance = 0.0001)


#Blood_Disorder
Blood_Disorder_metadata_filtered <- filter(MoL_metadata_filtered, MoL_metadata_filtered$Blood_Disorder== "1" | MoL_metadata_filtered$Blood_Disorder == "0")

table(Blood_Disorder_metadata_filtered$Blood_Disorder) 

Blood_Disorder_metadata_filtered$Tooth_modified = factor(Blood_Disorder_metadata_filtered$Tooth,
                                                         levels = c("Molar", "Bicuspid", "Canine", "Incisor"))

fit_data <- Maaslin2(
  input_data = data,
  input_metadata = Blood_Disorder_metadata_filtered,
  output = "C:/Users/aganc/Desktop/Professional/Projects/PSU_GraduateSchool_Projects/Strep_Paper/analyses/r_coding/MAASLIN2_Results_10072024/Blood_Disorder",
  fixed_effects = c('Blood_Disorder'),
  random_effects = c('Tooth_modified'),
  normalization = "CLR", 
  transform = "LOG",
  max_significance = .05,
  min_prevalence = 0.1,
  min_abundance = 0.0001)


#Non_specific_periostitis
Non_specific_periostitis_metadata_filtered <- filter(MoL_metadata_filtered, MoL_metadata_filtered$Non_specific_periostitis== "1" | MoL_metadata_filtered$Non_specific_periostitis == "0")

table(Non_specific_periostitis_metadata_filtered$Non_specific_periostitis) 

Non_specific_periostitis_metadata_filtered$Non_specific_periostitis <- as.factor(Non_specific_periostitis_metadata_filtered$Non_specific_periostitis)

Non_specific_periostitis_metadata_filtered$Tooth_modified = factor(Non_specific_periostitis_metadata_filtered$Tooth,
                                                                   levels = c("Molar", "Bicuspid", "Canine", "Incisor"))

fit_data <- Maaslin2(
  input_data = data,
  input_metadata = Non_specific_periostitis_metadata_filtered,
  output = "C:/Users/aganc/Desktop/Professional/Projects/PSU_GraduateSchool_Projects/Strep_Paper/analyses/r_coding/MAASLIN2_Results_10072024/Non_specific_periostitis",
  fixed_effects = c('Non_specific_periostitis'),
  random_effects = c('Tooth_modified'),
  normalization = "CLR", 
  transform = "LOG",
  max_significance = .05,
  min_prevalence = 0.1,
  min_abundance = 0.0001)

#Overall_VertebralPathology
Overall_VertebralPathology_metadata_filtered <- filter(MoL_metadata_filtered, MoL_metadata_filtered$Overall_VertebralPathology== "1" | MoL_metadata_filtered$Overall_VertebralPathology == "0")

table(Overall_VertebralPathology_metadata_filtered$Overall_VertebralPathology) 

Overall_VertebralPathology_metadata_filtered$Overall_VertebralPathology <- as.factor(Overall_VertebralPathology_metadata_filtered$Overall_VertebralPathology)

Overall_VertebralPathology_metadata_filtered$Tooth_modified = factor(Overall_VertebralPathology_metadata_filtered$Tooth,
                                                                     levels = c("Molar", "Bicuspid", "Canine", "Incisor"))

fit_data <- Maaslin2(
  input_data = data,
  input_metadata = Overall_VertebralPathology_metadata_filtered,
  output = "C:/Users/aganc/Desktop/Professional/Projects/PSU_GraduateSchool_Projects/Strep_Paper/analyses/r_coding/MAASLIN2_Results_10072024/Overall_VertebralPathology",
  fixed_effects = c('Overall_VertebralPathology'),
  random_effects = c('Tooth_modified'),
  normalization = "CLR", 
  transform = "LOG",
  max_significance = .05,
  min_prevalence = 0.1,
  min_abundance = 0.0001)


#Overall_Joint
Overall_Joint_metadata_filtered <- filter(MoL_metadata_filtered, MoL_metadata_filtered$Overall_Joint== "1" | MoL_metadata_filtered$Overall_Joint == "0")

table(Overall_Joint_metadata_filtered$Overall_Joint) 

Overall_Joint_metadata_filtered$Overall_Joint <- as.factor(Overall_Joint_metadata_filtered$Overall_Joint)

Overall_Joint_metadata_filtered$Tooth_modified = factor(Overall_Joint_metadata_filtered$Tooth,
                                                        levels = c("Molar", "Bicuspid", "Canine", "Incisor"))

fit_data <- Maaslin2(
  input_data = data,
  input_metadata = Overall_Joint_metadata_filtered,
  output = "C:/Users/aganc/Desktop/Professional/Projects/PSU_GraduateSchool_Projects/Strep_Paper/analyses/r_coding/MAASLIN2_Results_10072024/Overall_Joint",
  fixed_effects = c('Overall_Joint'),
  random_effects = c('Tooth_modified'),
  normalization = "CLR", 
  transform = "LOG",
  max_significance = .05,
  min_prevalence = 0.1,
  min_abundance = 0.0001)


###########
#Other Variables
###########

#Tooth
Tooth_metadata_filtered <- filter(metadata_filtered, metadata_filtered$Tooth== "Molar" | metadata_filtered$Tooth == "Bicuspid" | metadata_filtered$Tooth == "Incisor" | metadata_filtered$Tooth == "Canine")
table(Tooth_metadata_filtered$Tooth) 

Tooth_metadata_filtered$Tooth= factor(Tooth_metadata_filtered$Tooth,
                                      levels = c("Molar", "Bicuspid", "Incisor", "Canine"))

fit_data <- Maaslin2(
  input_data = data,
  input_metadata = Tooth_metadata_filtered,
  output = "C:/Users/aganc/Desktop/Professional/Projects/PSU_GraduateSchool_Projects/Strep_Paper/analyses/r_coding/MAASLIN2_Results_10072024/Tooth",
  normalization = "CLR", 
  transform = "LOG",
  fixed_effects = c('Tooth'),
  max_significance = .05,
  min_prevalence = 0.1,
  min_abundance = 0.0001)


#Dominant_Taxa
Dominant_Taxa_metadata_filtered <- filter(metadata_filtered, metadata_filtered$Dominant_Taxa== "Actino" | metadata_filtered$Dominant_Taxa == "Methanobrevibacter" | metadata_filtered$Dominant_Taxa == "Streptococcus")

table(Dominant_Taxa_metadata_filtered$Dominant_Taxa) 

Dominant_Taxa_metadata_filtered$Dominant_Taxa= factor(Dominant_Taxa_metadata_filtered$Dominant_Taxa,
                                                      levels = c("Actino", "Methanobrevibacter", "Streptococcus"))

fit_data <- Maaslin2(
  input_data = data,
  input_metadata = Dominant_Taxa_metadata_filtered,
  output = "C:/Users/aganc/Desktop/Professional/Projects/PSU_GraduateSchool_Projects/Strep_Paper/analyses/r_coding/MAASLIN2_Results_10072024/Dominant_Taxa",
  fixed_effects = c('Dominant_Taxa'),
  random_effects = c('Tooth_modified'),
  normalization = "CLR", 
  transform = "LOG",
  max_significance = .05,
  min_prevalence = 0.1,
  min_abundance = 0.0001)

# BlockPeriodSimplified
BlockPeriodSimplified_metadata_filtered <- filter(metadata_filtered, metadata_filtered$BlockPeriodSimplified != 'U')

table(BlockPeriodSimplified_metadata_filtered$BlockPeriodSimplified)

BlockPeriodSimplified_metadata_filtered$BlockPeriodSimplified= factor(BlockPeriodSimplified_metadata_filtered$BlockPeriodSimplified,
                                                                      levels = c("Pre_Roman_Britian", "Roman_Britain", "Anglo_Saxon_or_Early_Medieval_Britain", "Norman_Britain_and_the_Middle_Ages", "Reformation", "Industrial", "Modern"))

fit_data <- Maaslin2(
  input_data = data,
  input_metadata = BlockPeriodSimplified_metadata_filtered,
  output = "C:/Users/aganc/Desktop/Professional/Projects/PSU_GraduateSchool_Projects/Strep_Paper/analyses/r_coding/MAASLIN2_Results_10072024/BlockPeriodSimplified",
  fixed_effects = c('BlockPeriodSimplified'),
  random_effects = c('Tooth_modified'),
  normalization = "CLR", 
  transform = "LOG",
  max_significance = .05,
  min_prevalence = 0.1,
  min_abundance = 0.0001)

