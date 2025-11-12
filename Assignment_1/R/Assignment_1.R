# Assignment_1 Bartek Puzio 1146247 ####

# Load all packages for this script
# I took the time to add proper headers for each section of this assignment, please open and collapse each section as needed for a structured experience.
library(tidyverse)
library(tidyr)
library(vegan)
library(readr)
library(dplyr)
library(ggplot2)
library(maps)
library(styler)
library(rlang)

## Housekeeping ####

# Load data file into Rstudio, remove the # from the next line and run styler if you want because I already did:)
# setwd(Assignment_1.zip/Assignment_1/R)
# hopefully the working directory is appropriately set. The pathways used to call the data file are sound but needs to be set in the correct directory to work. I use the environment to find the file with the R script and using the blue gear "more" setting, set working directory.
# styler::style_file("Assignment_1.R")
data.b <- read_tsv(file = "../data/water_bear.tsv")

#Check for parsing issues (TA suggestion / quality check)
parse_issues <- problems(data.b)  # Lists rows with parsing issues
cat("Total rows with parsing issues:", nrow(parse_issues), "\n")
head(parse_issues, 20)

#Summary statistics for exploratory analysis to gain an understanding of data structure
summary_stats <- data.b %>%
  summarise(
    n_total = n(),
    n_bins = n_distinct(bin_uri),
    n_species = n_distinct(species, na.rm = TRUE),
    ratio_bins_species = n_bins / n_species
  )
print(summary_stats)

# Data manipulation + cleaning ####
# The Coordinates in the dataset columns had square brackets, a comma, and were not separated into 2 columns by latitude and longitude. This code removes the special characters and separates the coordinates into 2 columns and recognizes the values as numeric. The original cood column remains in the dataset for reproducabiltiy

data.b <- data.b %>%
  separate(coord, into = c("lat", "long"), sep = ",", remove = FALSE)

data.b$lat <- gsub("\\[", "", data.b$lat)
data.b$long <- gsub("\\]", "", data.b$long)

data.b$lat <- as.numeric(data.b$lat)
data.b$long <- as.numeric(data.b$long)

## Graph - Specimen collection site on Globe ####
# Map plotting with cleaned data to avoid NA errors. The map visualizes specimen collection sites across the globe, overlaying points on country polygons to show geographic distribution.
# Filter first, then pass to ggplot
data_clean <- data.b %>% 
  filter(!is.na(lat) & !is.na(long))

ggplot() +
  geom_polygon(
    data = world_map, aes(x = long, y = lat, group = group),
    fill = "lightblue", color = "gray70", alpha = 0.5
  ) +
  geom_point(
    data = data_clean,  # <- now guaranteed to have lat/long
    aes(x = long, y = lat),
    color = "red", alpha = 0.5
  ) +
  labs(title = "Specimen Locations on Globe", x = "Longitude", y = "Latitude") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

## Data objects for Graph 2 ####
# Statistics

data.b.ratio <- data.b %>% filter(!is.na(lat), lat > 0)

n_bins <- n_distinct(data.b.ratio$bin_uri, na.rm = TRUE)
n_species <- n_distinct(data.b.ratio$species, na.rm = TRUE)
ratio <- n_bins / n_species

c(n_bins = n_bins, n_species = n_species, ratio = ratio)

# Result: 1.81. This suggests that there are more distinctive clusters, or mitochondrial lineages, than named species in this dataset. This can be an indicator of cryptic species diversity or populations that have been geographically isolated for some time.



data.b.clean <- data.b %>%
  filter(
    !is.na(bin_uri),
    !is.na(`country/ocean`),
    !is.na(species),
    `country/ocean` != "Unrecoverable",
    !is.na(lat),
    lat >= 0
  )

bin_country <- data.b.clean %>%
  group_by(`country/ocean`) %>%
  summarise(
    n_specimens = n(),
    n_bins = n_distinct(bin_uri),
    n_species = n_distinct(species, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_bins))


## Graph - BIN richness vs Named species from Countries North of the Equator####
bin_country %>%
  slice_max(n_bins, n = 15) %>%
  select(`country/ocean`, n_bins, n_species) %>%
  pivot_longer(
    cols = c(n_bins, n_species),
    names_to = "metric", values_to = "value"
  ) %>%
  ggplot(aes(x = reorder(`country/ocean`, value), y = value, fill = metric)) + # https://r4ds.hadley.nz/data-visualize for help and Chatgpt Ai for debugging
  geom_col(position = "dodge") +
  coord_flip() +
  labs(
    x = "Country",
    y = "Count",
    fill = "",
    title = "BIN Richness vs. Named Species (Top 15 Countries, Equator & North)"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

  ## Define reusable histogram function ####

plot_lat_hist <- function(df, bin_col = "bin_uri", lat_filter = 0, title_suffix = "", unique_only = TRUE) {
  bin_sym <- rlang::sym(bin_col)
  
  df_plot <- df
  if (unique_only) {
    df_plot <- df_plot %>% distinct(!!bin_sym, .keep_all = TRUE)
  }
  
  df_plot %>%
    filter(lat > lat_filter) %>%
    ggplot(aes(x = lat)) +
    geom_histogram(bins = 100, fill = "blue", color = "white") +
    labs(
      title = paste("Distribution of", title_suffix),
      x = "Latitude (°)", y = "Frequency"
    ) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
## Graph - Distribution of BIN's by Latitude North of the Equator ####
  
plot_lat_hist(data.b, title_suffix = "BINs by Latitude (North of Equator)", unique_only = FALSE)

##Correlation test : provides a quantitative assessment of the relationship between geographic position and biodiversity
df_lat <- richness.bands %>%
  mutate(lat_mid = (as.numeric(sub("-.*", "", lat_band)) + as.numeric(sub(".*-", "", lat_band)))/2)

cor_test <- cor.test(df_lat$lat_mid, df_lat$BIN_richness)
cor_test  # ==> Gives correlation coefficient and p-value

#Results: correlation coefficient (r = 0.07) and a p-value of 0.88.This means the relationship between latitude and BIN richness is very weak and statistically non-significant.The small sample size and broad confidence interval indicate limited statistical power, meaning that even if a trend exists, the available data are insufficient to confirm it.

## Graph - distribution of unique bins north of the equator by latitude  ####
plot_lat_hist(data.b, title_suffix = "Unique BINs North of the Equator Latitude", unique_only = TRUE)

## Graph - BIN richness across latitude bands north of the equator ####

data.b.bands <- data.b %>%
  filter(!is.na(lat), lat >= 0) %>%
  mutate(
    lat_band = cut(
      lat,
      breaks = seq(0, 90, by = 10),
      include.lowest = TRUE, right = FALSE,
      labels = paste0(seq(0, 80, 10), "-", seq(10, 90, 10))
    )
  )

richness.bands <- data.b.bands %>%
  filter(!is.na(bin_uri)) %>%
  group_by(lat_band) %>%
  summarise(BIN_richness = n_distinct(bin_uri), .groups = "drop")

ggplot(richness.bands, aes(x = lat_band, y = BIN_richness)) +
  geom_col(fill = "blue") +
  theme_minimal() +
  labs(
    title = "Tardigrade BIN richness across latitude bands \n (Northern Hemisphere)",
    x = "Latitude band (°N)",
    y = "Number of unique BINs"
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 200))



