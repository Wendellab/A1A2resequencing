#!/usr/bin/env R
### SCRIPT DESCRIPTION ###
# Plot a horizontal barchart of membership proportions from STRUCTURE results.

## NOTES
# You may need to change `scale_fill_manual` in `plot_data_hbarchart` to match
#   the number of clusters in your analysis.
# You may need to change the interest variable to match your individuals of
#   interest.
# You may need to change the size of the final image.

## SECTIONS
# SETUP
# ARGUMENTS
# FUNCTIONS
# ANALYSIS
# DISCARDED CODE

## SETUP
library(stringr)
library(dplyr)
library(reshape2)
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
library(tibble)

## ARGUMENTS

# data_file is the STRUCTURE input data
# population_file is a tab-delimited file containing accession-population mappings
args = commandArgs(trailingOnly=TRUE)
data_file = args[1]
population_file = args[2]

# If you aren't running this script from the command-line, set these values to
#   your input data and comment the args lines above.
# data_file = "structure.2.1.strct_out_f"
# population_file = "populations.txt"

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
melt_data = function(population_info) {
  melt(population_info, id.vars = c("species", "population")) %>%
    mutate(variable = as.numeric(str_replace(variable, "V", "")) - 4,
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


# Plot a horizontal barchart of the accessions faceted by population.
# The y-axis is sorted by the order variable, which is assumed to be
#   an ordering based on how the accessions cluster. See below for more
#   information.
# The labels and colors are two variables that describe how the clusters
#   should be colored and labeled in the legend. If they aren't provided, then
#   ggplot with color the clusters automatically and label them with 1, 2, 3, ...
plot_data_hbarchart = function(data, order, interest, labels = NULL, colors = NULL) {

  if (length(unique(data$cluster)) != length(colors)) {
    stop("You do not have the same number of colors and clusters.")
  }

  if (length(unique(data$cluster)) != length(labels)) {
    stop("You do not have the same number of labels and clusters.")
  }

  if (length(colors) != length(labels)) {
    stop("You do not have the same number of labels and colors")
  }

  # Remove 0 from x-axis labels.
  perc_lab <- function(x) {
    if_else(x != 0,
            scales::label_percent()(abs(x)),
            ""
    )
  }

  # Reorganize the data based on the order variable.
  # The y-axis label of species in the interest variable are colored red.
  data = merge(data, order, by.x="species", by.y="order") %>%
    arrange(rank, cluster) %>%
    mutate(species = ifelse(species %in% interest,
                            glue("<span style = 'color:red'>{species}</span>"),
                            glue("<span style = 'color:black'>{species}</span>")
    )
    ) %>%
    mutate(species = fct_reorder(as.factor(species), rank))


  # Plot the data with membership proportion as the x-axis and species on the
  #   y-axis.
  barchart = ggplot(data, aes(x = percent, y = species)) +

    # geom_bar_h is a horizontal bar plot from ggstance. Using coord_flip() to
    # flip geom_bar has issues with the facet axes.
    geom_barh(aes(fill = fct_rev(as.factor(cluster))), stat = "identity") +

    # The line below adds labels to to bars for percents > 99 and != 100.
    # geom_text(aes(x = percent, label = label), size = 2, position = position_stack(vjust = 0.5)) +

    # set title and axes labels and legend title
    labs(title = paste0("Membership in ", length(data$cluster %>% unique()), " Inferred Clusters"),
         x = "Membership",
         y = "Species",
         fill = "Inferred Cluster") +

    # set x-tick labels
    scale_x_continuous(labels=perc_lab, expand=c(0,0)) +

    # set theme and elements
    theme_minimal() +
    theme(axis.text.y = element_markdown(),      # Color species of interest on y-axis
          legend.text = element_markdown()) +    # Render markdown in legend.
    guides(fill = guide_legend(reverse = TRUE)) +

    # Facet the barchart by population.
    # Scales and space are free so that there aren't gaps in the barcharts.
    facet_grid(population ~., scales = "free_y", space = "free_y")

  if (is.null(labels) & is.null(colors)) {
    return(barchart)
  } else if (is.null(labels) | is.null(colors)) {
    stop("You must set both labels and colors. You cannot set only one.")
  } else {
    return(barchart + scale_fill_manual("Legend", labels = labels, values = colors))

  }
}

## ANALYIS

# Get raw data for membership proportion from STRUCTURE file.
raw_data = get_data(data_file)

# Melt the raw data.
data = melt_data(raw_data)

# Cluster the data based on how similarity of membership proportion.
cluster = dist(raw_data %>%
                 select(-population) %>%
                 column_to_rownames("species")
)
dend <- hclust(cluster, "ave") %>% as.dendrogram()
order = raw_data[order.dendrogram(dend),]$species
order = as.data.frame(order) %>% rownames_to_column("rank") %>%
  mutate(rank = as.numeric(rank))

# Set accessions of interest. Set to c() if there are none.
# interest = c("A1037", "A1123", "A1148", "A2038")
interest = c()

# Set colors and labels.
# If you don't set these, ggplot2 will set the colors automatically.
# The labels will be cluster number.
# The colors were used when there were three clusters. When there are two
#   clusters, blue4 represented G. herbaceum.
#   G. arboreum            = green4
#   G. herbaceum           = blue4
#   G. herbaceum (Chinese) = royalblue2

# Notice the Markdown italics formatting for the labels.
labels = c("*G. arboreum*", "*G. herbaceum*")
colors = c("1" = "blue4", "2" = "green4")

# Plot barchart and save to same location as data_file but with ".barchart.png"
#   extension.
ggsave(sub(pattern = "(.*)\\..*$", replacement = "\\1.barchart.png", data_file),
       plot_data_hbarchart(data, order, interest, labels, colors),
       width = 15,
       height = 15)

## DISCARDED CODE
# These functions were originally used to assign colors to each population
#   and produce a dendrogram of the clustering.
# The information in the cluster dendrogram was eventually added to the barchart
#   by sorting the y-axis labels by the clustering order, and this code was
#   moved to this section for documentation purposes.
# There are no guarantees that this code works with the code above.

# get colors from ColorBrewer based on number of clusters
get_colors = function(population_in) {
  populations = read.table(population_in, sep = "\t", header = T)

  # ColorBrewer complains if you choose a set that has a minimum of three colors
  # when you only have two clusters
  if (populations %>% select(population) %>% unique() %>% nrow() > 2){
    data.frame(population = populations %>% select(population) %>% unique() %>% pull(),
               color = I(brewer.pal(populations %>% select(population) %>% unique() %>% nrow(),
                                    name = 'Set1')))
  } else {
    data.frame(population = populations %>% select(population) %>% unique() %>% pull(),
               color = c("#E41A1C", "#377EB8"))
  }
}

# plot a dendrogram showing how species cluster
plot_data_dendro = function(raw_data, population_info) {

  # get colors for color bar
  bar_colors = population_info %>%
    select(species, color) %>%
    arrange(species) %>%
    unique() %>%
    select(color) %>%
    pull()

  # calculate distance of data
  cluster = dist(raw_data %>%
                   select(-population) %>%
                   column_to_rownames("species")
  )

  # make dendrogram
  dend <- hclust(cluster, "ave") %>% as.dendrogram()
  summary(dend)
  raw_data[order.dendrogram(dend),]$species
  # plot dendrogram
  dendro_plot = function() {
    par(mar=c(15,5,1,1))
    plot.new()
    dend %>%
      set("leaves_pch", 19)  %>%
      set("nodes_cex", 0.7) %>%
      plot(axes=FALSE)
    colored_bars(colors = bar_colors, dend = dend, rowLabels = "Population")
    legend("topright", legend = unique(population_info$population),
           fill=unique(population_info$color), cex=0.8)
  }

  # store plot in variable and return
  grid.echo(dendro_plot)
  p = grid.grab()
  dev.off()
  return(p)
}
