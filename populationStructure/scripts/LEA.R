#!/usr/bin/env R
### SCRIPT DESCRIPTION ###
# Run LEA or load a saved analysis and plot a horizontal barchart of membership
# proportions from "best" results.

## NOTES
# This script is likely best run interactively if you want to set the colors
#   to something specific. If you set `scale_fill_manual` to n cluster/color
#   pairings and don't have that many clusters, then the script will fail.
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
library(tibble)
library(dplyr)
library(LEA)
library(reshape2)
library(stringr)
library(ggplot2)
library(forcats)
library(ggstance)
library(RColorBrewer)
library(dendextend)
library(gridGraphics)
library(ggtext)
library(glue)

## ARGUMENTS

# data_file is the VCF input data
# population_file is a tab-delimited file containing accession-population mappings
args = commandArgs(trailingOnly=TRUE)
data_file = args[1]
population_file = args[2]
num_cpus = args[3]

# If you aren't running this script from the command-line, set these values to
#   your input data and comment the args lines above.
# data_file = "../all/A1A2.all.LD_filtered.vcf"
# population_file = "populations.all.txt"
# num_cpus = 16

## FUNCTIONS

# melt data to get individual rows per cluster
# keep species and population
# convert VN to N (V1 to 1, V2 to 2, ...) as cluster names
# if value < 1%, then don't label it
# else label it as rounded percent
melt_data = function(population_info) {
  melt(population_info, id.vars = c("species", "population")) %>%
    mutate(variable = as.numeric(str_replace(variable, "V", "")),
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

## ANALYSIS

# Strip the extension from data_file for later use.
file_base = file.path(sub(pattern = "(.*)\\..*$", replacement = "\\1", data_file))

# Write .geno file in working directory using data_file name for LEA analysis.
# vcf2geno(data_file, paste0(file_base, ".geno"))

# Run analysis on the .geno file for K = 1 - K = 10 with 10 reps per K and 16 CPUs.
project = snmf(geno_file,
               K = 1:10,
               entropy = T,
               repetitions = 10,
               project = "new",
               CPU = num_cpus)

# Alternatively, if the analysis has already been done, load the .snmfProject file.
# Loading the data from a saved project requires that the .snmfProject file and
#   the .snmf folder be in the same location.
# This script saves them to file_base.snmf and file_base.snmfProject.
# project = load.snmfProject("../all/A1A2.all.LD_filtered.snmfProject")

# Get best Ks by finding minimum min, mean, and max cross entropy
summary_info = summary(project)$crossEntropy

# Plot the mean cross entropy.
# Get the mean from summary_info and convert the rownames to column "K".
# Remove "K = " from the K column.
# Rename the "summary_info[2, ]" column to "mean".
# Pass this information to ggplot.
# Add a horizontal red line at the minimum cross entropy.
# Plot the cross entropy by K as a scatterplot.
# Set the x-axis labels from 1-10.
# Change the x-axis and y-axis labels and the title.
# Remove the legend.
summary_plot = as.data.frame(summary_info[2,]) %>%
  rownames_to_column("K") %>%
  mutate(K = as.numeric(str_replace(K, "K = ", ""))) %>%
  rename(mean = `summary_info[2, ]`) %>%
  ggplot() +
  geom_hline(aes(yintercept = min(summary_info[2,]), color = "red")) +
  geom_point(aes(x = K, y = mean)) +
  scale_x_continuous(breaks = c(1:10)) +
  labs(title = "Cross-entropy versus K", x = "K", y = "mean cross-entropy") +
  theme(legend.position = "none")

# Save the summary image to file_base.summary.png.
ggsave(paste0(file_base, ".summary.png"),
       summary_plot,
       width = 5,
       height = 5)

# Find the smallest cross-entropy for mininum, mean, and maximum cross-entropy.
# Remove duplicates.
# Find the smallest K among these values.
# You could also set this manually.
K = min(
  c(which.min(summary_info[1,]),
    which.min(summary_info[2,]),
    which.min(summary_info[3,])) %>%
    unique()
)

# Find the replicate with the smallest cross-entropy for the selected value of K.
best = which.min(cross.entropy(project, K = K))

# Get population information.
species_info = read.table(population_file, header = T) %>% rename(species = individual)

# Add cluster proportions to species_info.
raw_data = as.data.frame(cbind(species_info, Q(project, K, best)))

# Melt cluster info.
population_info = melt_data(raw_data)

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
# You'll probably need to save the file without any special labels or colors,
#   identify the clusters, and then set labels and colors.
# labels = c("*G. arboreum*", "*G. herbaceum*", "*G. herbaceum*")
# colors = c("1" = "blue4", "2" = "green4", "3" = "royalblue2")

# Plot barchart and save to same location as data_file but with
# ".K(K).barchart.png" extension, where (K) is the smallest K.
ggsave(paste0(file_base, ".K" , K, ".test.barchart.png"),
       # plot_data_hbarchart(population_info, order, interest, labels, colors),
       plot_data_hbarchart(population_info, order, interest),
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
