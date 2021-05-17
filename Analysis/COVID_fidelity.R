# Load packages
library(easypackages)
library(tidyverse)
library(ggpubr)
library(RSKC)
library(simEd)
library(httr)
library(Rtsne)
library(factoextra)

### Import and Prepare Datasets ###
# Load full dataset
all_fidelity <- read.csv("~/Documents/Lab/COVID/ALL_Fidelity.csv")

# Make a new dataframe containing only genes of interest
#   Filter full dataset to include only the genes of interest 
#   Drop the columns that correspond to "Entrez" and "Alias"
genes_of_interest <- c("ACE2",
                       "DPP4",
                       "TMPRSS2",
                       "NRP1",
                       "ITGB3")

fidelity_subset <- all_fidelity %>%
  filter(Gene %in% genes_of_interest) %>%
  select(-Entrez, -Alias) # Drop other columns

# mgimond.github.io - advanced pivot_longer options
df2.long <- pivot_longer(fidelity_subset, cols = -Gene,
                         names_to = c("Brain.Region", "Cell.Subtype"),
                         names_sep = "_",
                         values_to = "Fidelity")

df2.wide <- pivot_wider(df2.long,
                        names_from = c("Gene", "Cell.Subtype"),
                        values_from = "Fidelity",
                        names_sep = "_")

# Define a vector of brain regions of interest
#   Then use it to choose only the desired brain regions from the data frame
regions_of_interest <- c("AMY", "BF", "CB", "CLA", "DI", 
                         "FCX", "GP", "HIP", "IN", "LIM", 
                         "MED", "MID", "OCX", "PCX", "PON", 
                         "STR", "TCX", "WM")

df_fidelity <- df2.wide %>% 
  filter(Brain.Region %in% regions_of_interest) %>%
  slice(match(regions_of_interest, Brain.Region)) %>% # Reorder rows to match entries in "regions" vector
  select(-contains("Percentile")) # Drop columns with NA values


### Exploratory Analysis ###
plot_df <- df2.long %>% 
  filter(Gene == "ACE2" | Gene == "TMPRSS2") %>% 
  filter(Brain.Region == "FCX")

ggplot(plot_df, aes(x = Cell.Subtype, y = Fidelity, fill = Gene)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Cell Type") + 
  theme_minimal()


#####################################
## RSKC
#####################################

# Use df_fidelity2 for RSKC; has brain region as rows and gene_celltype as columns
df_fidelity2 <- column_to_rownames(df_fidelity, "Brain.Region")

while (T) {
  
  # Assign the values of 2, 4, 6 and 8 to 'clust_vect'.
  clust_vect <- c(2,3,4,5,6,7,8)
  
  # Assign an empty list to 'rskc_list'.
  rskc_list <- list()
  
  # Assign a value of 0 to 'counter'.
  counter <- 0
  
  # Assign an empty list to 'weight_list'.
  weight_list <- list()
  
  # For 'i' -- the current number of clusters -- in 'clust_vect'...
  for (i in clust_vect) {
    # i = 2
    # Increment 'counter' with a value of 1.
    counter = counter + 1
    
    # Perform RSKC for whatever-the-value-of-'i'-is many clusters using 'myProt'.
    #    Assign RSKC's output as an entry in 'rskc_list'. 
    rskc_list[[counter]] <- RSKC(df_fidelity2, 
                                ncl = i,
                                alpha = 0.1,
                                L1 = NULL)
    
    
    # Convert the row names of 'myProt' to a column
    #   called 'Identifier' and store it in 'proteins_and_ids'.
    gene_and_region <- df_fidelity2 %>% 
      rownames_to_column("Brain.Region")
    
    # For the current object in 'rskc_list' convert the cluster labels
    #    into characters, and assign them to a new column called 'cluster_labels'
    #    in 'protein_and_ids'.  
    gene_and_region$cluster_labels <- rskc_list[[counter]]$labels %>% 
      as.character()
    
    # Order the weights for the current item in 'rskc_list' from largest
    #   to smallest, extract the names of the proteins in this order,
    #   convert this object into a data frame, and store this info in
    #   an object 'weight_df' in a column called 'protein'.
    weight_df <- sort(rskc_list[[counter]]$weights,
                      decreasing = T) %>% 
      names() %>% 
      
      as.data.frame() %>% 
      
      rename('protein' = ".")
    
    # Assign the ordered weights for the current item in 'rskc_list' into
    #  'weight_df', in a column called 'weight', 
    weight_df$weight <- sort(rskc_list[[counter]]$weights,
                             decreasing = T) %>% 
      
      unname() 
    
    # Impose a factor order on the contents of 'weight_df$protein' in 
    #   the current order. 
    weight_df$protein <- factor(weight_df$protein,
                                weight_df$protein)
    
    # Create a bar graph of the RSKC weights for each protein ordered from
    #   biggest to largest. Assign this graph to an object, 'weight_bars'
    weight_bars <- weight_df %>% 
      
      ggplot(aes(x = protein,
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
    
    # Print the current object in 'weight_bars' as a tiff in your working
    #   directory.
    png(paste0('RSKC_weights_for7_synprot_k=',i,'.png'),
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


#####################################
## Elbow Plot
#####################################



#####################################
## tSNE
#####################################

# Create vector of the weights obtained from RSKC and assign to 'weights'
#   and make empty matrix 'weighted_fidelity' for new weighted fidelity scores 
weights <- as.matrix(rskc_list[[1]]$weights)
weighted_fidelity <- matrix(nrow = 18, ncol = 20)

# Multiply df_fidelity by corresponding weights obtained from RSKC
for (i in 1:20){
  weighted_fidelity[,i] <- df_fidelity2[,i]*weights[i]
}

# Run tsne on weighted fidelity scores, and assign to 'tsne'
tsne <- Rtsne(weighted_fidelity, perplexity = 5)

# tsne_out: the two dimensions and corresponding regions
tsne_out <- tsne$Y %>%
  data.frame(regions_of_interest) %>%
  rename(Brain.Region = regions_of_interest, V1 = X1, V2 = X2) #rename columns

while (T) {
  
  # Have to manually change the number of clusters to be identified
  #   and has to match the number previously assigned to 'clust_vect' in above while loop for RSKC
  clust_vect <- c(8)
  
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
  
  # Assign an empty list to 'tsne_list'.
  tsne_list <- list()
  
  # For 'i' -- the current number of clusters -- in 'clust_vect'...
  for (i in clust_vect) {
    # i = 2
    # Increment 'counter' with a value of 1.
    counter = counter + 1
    
    # Merge 'tsne_out' with 'gene_and_region' according
    #   to their shared 'Brain.Region' column, and assign to 'tsne_genes_regions_clusts'. 
    tsne_genes_regions_clusts <- merge(tsne_out,
                                       gene_and_region,
                                       by = "Brain.Region")
    
    
    # Create a tSNE scatter plot where each point is colour-coded according to
    #   its designated RSKC cluster and assign this figure to 'tsne_scatter'.
    tsne_scatter <- ggplot(tsne_genes_regions_clusts,
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
            legend.title.align=0.5)+
      guides(fill=guide_legend(nrow = 2,
                               ncol = 4, 
                               byrow = TRUE)) +
      labs(x="V1", y="V2") 
    
    # Assign 'tsne_scatter' as an entry in 'tsne_list'.
    tsne_list[[counter]] <- tsne_scatter

    # Print the current object in 'tsne_scatter' as a tiff in your working
    #   directory.
    png(paste0('RSKC_scatter_k=',i,'.png'), 
        units = "in",
        width = 9, 
        height = 8,
        res = 300#,compression = 'lzw'
    )
    
    print(tsne_scatter)
    
    dev.off()
    
  }
  
  # break out of the while-loop when for-loop is done. 
  break
  
}
