# This script was used for a specific comparison image: STRUCTURE K2 versus LEA 
#   K3 for the A1A2 filtered dataset. It is provided for the sake of documentation
#   but no effort has been made to make this script generic like LEA.R and 
#   STRUCTURE.R.
# The code in this script should be considered more of an example on how the
#   combined image was generated rather than a script that can be run with no
#   modification.
# Users of this script will need to change:
#   lines 37-45:   input data: data_file, population_file, LEA_project
#   line 96:       species of interest
#   lines 102-108: labels and colors
#   line 151:      selected K for LEA
#   line 208:      output file name and dimensions
# Users will need to verify
#   line 175:      fuzzy match between truncated STRUCTURE names and LEA names
# Users will need to decide
#   lines 198-199: whether to include LEA y-axis labels

library(stringr)
library(dplyr)
library(tibble)
library(reshape2)
library(gtools)
library(ggplot2)
library(forcats)
library(ggstance)
library(ggpubr)
library(RColorBrewer)
library(dendextend)
library(gridGraphics)
library(ggplotify)
library(ggtext)
library(glue)
library(LEA)

# STRUCTURE A1A2 filtered K = 2, rep = 1
data_file = "structure.2.1.strct_out_f"

# Population data about the filtered VCF
population_file = "../populations.txt"   

# LEA run for A1A2 filtered + plink thinning
# make sure "projDir" in the .snmfProject file 
# is correct relative to the .snmf directory
LEA_project = "../filtered.snmfProject"

## FUNCTIONS

# read STRUCTURE output and remove everything before and after
# the membership proportion data, then read the membership
# proportion data as a table and drop the row number and 
# missing % columns
get_data = function(structure_in) {
  text = readChar(structure_in, file.info(structure_in)$size)
  individual_info = gsub(":", "", 
                         sub(" \n\n\n.*", "", 
                             sub(".*Inferred clusters\n  ", "", text)
                         )
  )
  individual_info = read.table(text = individual_info) %>% mutate(species = V2,
                                                                  population = V4) %>%
    select(-V2, -V4, -V1, -V3)
  return(individual_info)
}

# melt data to get individual rows per cluster
# keep species and population
# convert VN to N (V1 to 1, V2 to 2, ...) as cluster names
# if value < 1%, then don't label it
# else label it as rounded percent
melt_data = function(population_info, analysis_type) {
  melt(population_info, id.vars = c("species", "population")) %>% 
    mutate(variable = as.numeric(str_replace(variable, "V", "")) - ifelse(analysis_type == "STRUCTURE", 4, 0),
           label = ifelse(value < 0.01, 
                          "", 
                          ifelse (value == 1, 
                                  "",
                                  paste0(round(value*100, 3), "%")
                          )
           )
    ) %>%
    select(species, population, cluster=variable, percent=value, label)
}

# plot a horizontal bar chart
plot_data_hbarchart = function(data, order, analysis_type) {
  # remove 0 from x-label
  perc_lab <- function(x) {
    if_else(x != 0, 
            scales::label_percent()(abs(x)), 
            ""
    )
  }
  data = merge(data, order, by.x="species", by.y="order") %>%
    arrange(rank, cluster) %>% 
    mutate(species = ifelse(species %in% c("A1037", "A1123", "A1148", "A2038"),
                            glue("<span style = 'color:red'>{species}</span>"),
                            glue("<span style = 'color:black'>{species}</span>"))
    ) %>%
    mutate(species = fct_reorder(as.factor(species), rank))
  
  if(analysis_type == "STRUCTURE") {
    labels = c("*G. arboreum*", "*G. herbaceum*")
    colors = c("1" = "blue4", "2" = "green4")
  } else {
    labels = c("*G. arboreum*", "*G. herbaceum*", "*G. herbaceum*")
    colors = c("1" = "royalblue2", "2" = "blue4", "3" = "green4")
  }
  
  
  ## plot the data with membership proportion as the x-axis and species as the y-axis
  ggplot(data, aes(x = percent, y = species, rank)) +
    
    # geom_bar_h is a horizontal bar plot without coord_flip(), which doesn't work with free
    # scales and space in the facets
    geom_barh(aes(fill = fct_rev(as.factor(cluster))), stat = "identity") +
    scale_fill_manual("Legend", labels = labels, values = colors) +
    
    # add labels to to bars for percents > 99 and != 100
    # geom_text(aes(x = percent, label = label), color = "white", size = 2, position = position_stack(vjust = 0.5)) +
    
    # set title and axes labels and legend title
    labs(title = paste0("Membership in ", length(data$cluster %>% unique()), " Inferred Clusters"), 
         x = "Membership", 
         y = "Species", 
         fill = "Inferred Cluster") +
    
    # set x-tick labels
    scale_x_continuous(labels=perc_lab, expand=c(0,0)) +
    
    # set theme and elements
    theme_minimal() + 
    theme(axis.text.y = element_markdown(),
          legend.text = element_markdown()) +
    guides(fill = guide_legend(reverse = TRUE)) +
    
    # facet by population
    facet_grid(population ~., scales = "free_y", space = "free_y")
}

# get raw data for membership proportion from STRUCTURE file
STRUCTURE_data = get_data(data_file)

# melt and reformat the data
STRUCTURE_data = melt_data(STRUCTURE_data, "STRUCTURE")

# order data by dendrogram of LEA equivalent
project = load.snmfProject(LEA_project)

# Change this K if running something different.
best = which.min(cross.entropy(project, K = 3))

# add population info to individuals in VCF order
species_info = read.table(population_file, header = T) %>% 
  rename(species = individual)

# add cluster proportions to species_info
LEA_data = as.data.frame(cbind(species_info, Q(project, 3, best)))

# get order for LEA data
cluster = dist(LEA_data %>% 
                 select(-population) %>% 
                 column_to_rownames("species")
)
dend <- hclust(cluster, "ave") %>% as.dendrogram()
order = LEA_data[order.dendrogram(dend),]$species
order = as.data.frame(order) %>% rownames_to_column("rank") %>%
  mutate(rank = as.numeric(rank))

# melt cluster info
LEA_data = melt_data(LEA_data, "LEA")

# rebuild STRUCTURE name truncation
# verify this is a correct mapping before proceeding
non_truncated = stringdist_left_join(STRUCTURE_data %>% filter(str_detect(species, "^S")), 
                                     order %>% filter(str_detect(order, "^S")), 
                                     by = c("species" = "order"), max_dist = 4) %>%
  mutate(species = NULL) %>%
  rename(species = order)

# add LEA order to STRUCTURE data
STRUCTURE_data = rbind(left_join(STRUCTURE_data, order, by = c("species" = "order")) %>% filter(is.na(rank) == F), 
                       non_truncated) %>% 
  select(-rank)

# add this to the barcharts to move the legend in the combined image
move_legend = theme(legend.position = "bottom", legend.direction = "horizontal")

# add this to the LEA plot to remove the y-axis labels in the combined image
remove_y = rremove("y.text") + rremove("y.title")

# build barcharts
LEA_barchart = plot_data_hbarchart(LEA_data, order, "LEA")
STRUCTURE_barchart = plot_data_hbarchart(STRUCTURE_data, order, "STRUCTURE")

# make the comparison, use the more detailed legend from the LEA plot
comparison_barchart = ggarrange(STRUCTURE_barchart + move_legend,
                                LEA_barchart + move_legend + remove_y, # No y-axis labels for LEA
                                # LEA_barchart + move_legend, # Keep y-axis labels for LEA
                                labels = c("STRUCTURE", "LEA"),
                                vjust = 0,
                                nrow = 1, 
                                legend = "bottom", 
                                common.legend = T,
                                legend.grob = get_legend(LEA_barchart + move_legend)) + 
  theme(plot.margin = margin(1, 0, 0, 0, "cm")) 

ggsave("comparison.filtered.STRUCTURE_K2.LEA_K3.png", comparison_barchart, width = 20, height = 15)