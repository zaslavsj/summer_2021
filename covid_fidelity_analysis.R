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
  filter(Brain.Region %in% regions) %>% 
  slice(match(regions, Brain.Region)) # Reorder rows to match entries in "regions" vector

# Drop columns containing missing values.
fidelity <- fidelity_subset.wide2[, colSums(is.na(fidelity_subset.wide2)) < nrow(fidelity_subset.wide2)]

# Use the brain regions to name the rows and remove the brain region column.
myFidelity <- fidelity %>% 
  column_to_rownames("Brain.Region")

## Prepare a dataset containing all COVID-related genes
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
# Remove row corresponding to 'ALL'.
fidelity_full <- fidelity_full %>% 
  filter(!(Brain.Region %in% "ALL"))

# Use the brain regions to name the rows and remove the brain region column.
myFidelity_full <- fidelity_full %>% 
  data.frame(row.names = fidelity_full$Brain.Region) %>% 
  select(-1)

#### Exploratory Analysis
# Use data only related to the FCX region and ACE2, TMPRSS2 genes.
plot_df <- fidelity_subset.long %>% 
  filter(Gene == "ACE2" | Gene == "TMPRSS2") %>% 
  filter(Brain.Region == "FCX")

# Rough plot to check that we can reproduce previously made plots.
(ggplot(plot_df, aes(x = Cell.Subtype, y = Fidelity, fill = Gene)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal())

## tSNE
# Use fidelity for tSNE; it has brain regions as a column, not row names
tsne <- Rtsne(as.matrix(fidelity[,2:ncol(fidelity)]), perplexity = 1)

# tsne_out: the two dimensions and corresponding regions
tsne_out <- tsne$Y %>%
  data.frame(regions) %>%
  rename(Brain.Region = regions, V1 = X1, V2 = X2) # Rename columns


## RSKC

while (T) {
  
  # Assign the values of 2, 4, 6 and 8 to 'clust_vect'.
  clust_vect <- c(2,3,4)
  
  # Assign an empty list to 'rskc_list'.
  rskc_list <- list()
  
  # Assign 8 colours to 'col_vect'.
  col_vect <- c("#FF0000",
                "#0000FF",
                "#00FF00",
                "#A020F0",
                "#FFA500",
                "#FFFF00",
                "#A65628",
                "#F781BF")
  
  # Assign a value of 0 to 'counter'.
  counter <- 0
  
  # Assign an empty list to 'boxplot_list'.
  boxplot_list <- list()
  
  # Assign an empty list to 'tsne_list'.
  tsne_list <- list()
  
  # Assign an empty list to 'weight_list'.
  weight_list <- list()
  
  # For 'i' -- the current number of clusters -- in 'clust_vect'...
  for (i in clust_vect) {
    # i = 2
    # Increment 'counter' with a value of 1.
    counter = counter + 1
    
    # Perform RSKC for whatever-the-value-of-'i'-is many clusters using 'myFidelity'.
    # Assign RSKC's output as an entry in 'rskc_list'. Use myFidelity for 
    # RSKC; it has brain region as rows and gene_celltype as columns.
    rskc_list[[counter]] <- RSKC(myFidelity, 
                                 ncl = i,
                                 alpha = 0.1,
                                 L1 = sqrt(ncol(myFidelity)))
    
    
    # Convert the row names of 'myFidelity' to a column
    # called 'Brain.Region' and store it in 'gene_and_regions'.
    gene_and_region <- myFidelity %>% 
      rownames_to_column("Brain.Region")
    
    # For the current object in 'rskc_list' convert the cluster labels
    # into characters, and assign them to a new column called 'cluster_labels'
    # in 'gene_and_region'.  
    gene_and_region$cluster_labels <- rskc_list[[counter]]$labels %>% 
      as.character()
    
    # Merge the first 4 columns of 'tsne_out' with 'gene_and_region' according
    # to their shared 'Brain.Region' column, and assign to 
    # 'tsne_gene_region_clusts'. 
    tsne_gene_region_clusts <- merge(tsne_out,
                                     gene_and_region,
                                     by = "Brain.Region")
    
    
    # Create a tSNE scatter plot where each point is colour-coded according to
    # its designated RSKC cluster and assign this figure to 'tsne_scatter'.
    tsne_scatter <- ggplot(tsne_gene_region_clusts,
                           aes(V1,
                               V2,
                               fill = cluster_labels)) +
      geom_point(shape = 21, size = 3)+ 
      scale_fill_manual(values = col_vect[1:i],
                        labels = 1:i,
                        name = "Cluster") +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.background = element_rect(fill = NA, 
                                            colour = "white"),
            panel.border = element_blank(),
            axis.line = element_line(),
            legend.position = 'bottom',
            legend.background = element_rect(fill = NA,
                                             colour = NA), 
            legend.title.align=0.5) +
      guides(fill=guide_legend(nrow = 2,
                               ncol = 4, 
                               byrow = TRUE)) +
      labs(x="V1",
           y="V2") 
    
    # Assign 'tsne_scatter' as an entry in 'tsne_list'.
    tsne_list[[counter]] <- tsne_scatter
    
    # Order the weights for the current item in 'rskc_list' from largest
    # to smallest, extract the names of the proteins in this order,
    # convert this object into a data frame, and store this info in
    # an object 'weight_df' in a column called 'gene_celltype'.
    weight_df <- sort(rskc_list[[counter]]$weights,
                      decreasing = T) %>% 
      names() %>% 
      
      as.data.frame() %>% 
      
      rename('gene_celltype' = ".")
    
    # Assign the ordered weights for the current item in 'rskc_list' into
    # 'weight_df', in a column called 'weight'.
    weight_df$weight <- sort(rskc_list[[counter]]$weights,
                             decreasing = T) %>% 
      
      unname() 
    
    # Impose a factor order on the contents of 'weight_df$gene_celltype' in 
    # the current order. 
    weight_df$gene_celltype <- factor(weight_df$gene_celltype,
                                      weight_df$gene_celltype)
    
    # Create a bar graph of the RSKC weights for each gene and cell type 
    # ordered from largest to smallest. Assign this graph to an object, 
    # 'weight_bars'.
    weight_bars <- weight_df %>% 
      
      ggplot(aes(x = gene_celltype,
                 y = weight)) + 
      
      theme_classic() +
      
      geom_bar(stat = 'identity') +
      
      theme(axis.text.x = element_text(angle = 45,
                                       vjust = 0.6)) +
      
      scale_y_continuous(expand = c(0,
                                    0)) +
      
      ggtitle(paste0("RSKC Weights for K = ",i)) +
      
      xlab("") +
      
      ylab("Weights\n")
    
    # Assign 'weight_bars' as the current entry into 'weight_list'. 
    weight_list[[counter]] <- weight_bars
    
    # Print the current object in 'tsne_scatter' as a tiff in your working
    # directory.
    png(paste0('RSKC_scatter_k=',i,'.png'), 
        units = "in",
        width = 9, 
        height = 8,
        res = 300#,compression = 'lzw'
    )
    
    print(tsne_scatter)
    
    dev.off()
    
    # Print the current object in 'weight_bars' as a tiff in your working
    # directory.
    png(paste0('RSKC_weights_for5_genes_k=',i,'.png'),
        units = "in",
        width = 5+2,
        height = 4+1.5,
        res = 300 #,compression = 'lzw'
    )  
    
    print(weight_bars)
    
    dev.off()
    
  }
  
  # break out of the while-loop when for-loop is done. 
  break
  
}
