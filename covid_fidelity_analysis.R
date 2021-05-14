## Relevant Packages
# Load packages
library(easypackages)
library(tidyverse)
library(ggpubr)
library(RSKC)
library(simEd)
library(httr)
library(Rtsne)

## Import and Prepare Dataset
# Load full dataset
all_fidelity <- read.csv("~/Documents/R/Murphy_Lab/COVID-19/Data/ALL_Fidelity.csv")

# Define a vector containing the genes of interest.
genes_subset <- c("ACE2",
                  "DPP4",
                  "TMPRSS2",
                  "NRP1",
                  "ITGB3")

genes_full_list <- c("ACE2", "DPP4", "ANPEP", "CD209", "ENPEP", "CLEC4G", 
                     "CLEC4M", "CEACAM-1", "TMPRSS2", "TMPRSS11d", 
                     "CTSL", "ADAM17", "FURIN", "ITGA5", "ITGB3", "ITGA4",
                     "ITGB6", "ITGA7")

# Filter full dataset to include only the genes of interest and drop the 
# columns that correspond to "Entrez" and "Alias".
fidelity_subset <- all_fidelity %>%
  filter(Gene %in% genes_subset) %>% 
  select(-c(2,3)) 

# Tidy the data frame by separating the previous columns into
# two new columns for the brain region and cell subtype. Assign each of the
# values to a new column for the fidelity scores. The tidy data has each
# column correspond to a varaible and each observation is a new row.
fidelity_subset.long <-  pivot_longer(fidelity_subset, cols = -Gene,
                                      names_to = c("Brain.Region", "Cell.Subtype"),
                                      names_sep = "_",
                                      values_to = "Fidelity")

# Pivot the data frame to a wide view by assigning new columns corresponding to 
# the gene_associated cell subtype for each brain region. Assign the fidelity 
# scores as the values.
fidelity_subset.wide <- pivot_wider(fidelity_subset.long,
                                    names_from = c("Gene", "Cell.Subtype"),
                                    values_from = "Fidelity",
                                    names_sep = "_")

# Define a vector of brain regions of interest.
regions <- c("FCX", "PCX", "OCX", "TCX", "STR", "DI", "HIP", "AMY", "CB")

# Choose only the desired brain regions from the data frame.
fidelity_subset.wide2 <- fidelity_subset.wide %>% 
  filter(Brain.Region %in% regions)

# Drop columns containing missing values.
fidelity <- fidelity_subset.wide2[, colSums(is.na(fidelity_subset.wide2)) < nrow(fidelity_subset.wide2)]
# Use the brain regions to name the rows and remove the brain region column.
myFidelity <- fidelity %>% 
  column_to_rownames("Brain.Region")

## Prepare a dataset containing all genes
# Drop the columns corresponding to Entrez and Alias.
all_fidelity2 <- all_fidelity %>% 
  select(-c(2,3))

all_fidelity2.long <- pivot_longer(all_fidelity2, cols = -Gene,
                                   names_to = c("Brain.Region", "Cell.Subtype"),
                                   names_sep = "_",
                                   values_to = "Fidelity")

all_fidelity2.wide <- pivot_wider(all_fidelity2.long, 
                                  names_from = c("Gene", "Cell.Subtype"),
                                  values_from = "Fidelity",
                                  names_sep = "_")

# Remove the last row corresponding to mean expression.
all_fidelity2.wide2 <- all_fidelity2.wide[-22, ]

# Drop columns containing missing values.
fidelity_full <- all_fidelity2.wide2[, colSums(is.na(all_fidelity2.wide2)) < nrow(all_fidelity2.wide2)]
# Use the brain regions to name the rows and remove the brain region column.
myFidelity_full <- fidelity_full %>% 
  data.frame(row.names = fidelity_full$Brain.Region) %>% 
  select(-1)

#### Exploratory Analysis
# Use data only related to the FCX region and ACE2, TMPRSS2 genes.
plot_df <- fidelity_subset.long %>% 
  filter(Gene == "ACE2" | Gene == "TMPRSS2") %>% 
  filter(Brain.Region == "FCX")

ggplot(plot_df, aes(x = Cell.Subtype, y = Fidelity, fill = Gene)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal()