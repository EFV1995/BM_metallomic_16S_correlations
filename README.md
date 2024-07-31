# BM_metallomic_16S_heatmap correlations
# This repository contains the script utilized to perform the heatmap with metal and metalloids and microbiota abundances including diversity indices.

# Load required libraries
library(dplyr)
library(tibble)
library(forcats)
library(ComplexHeatmap)
library(ggplot2)

# Data frames are split_data_frames$File_3_micro$early and split_data_frames$File_2_mineral$early
microbiota <- x1
dim(x1)
Diversities <- Infant_apha_diversities
metals <- File_2_mineral
head(File_3_micro)
head(Infant_apha_diversities)
head(File_2_mineral)

# Cleaning data
microbiota$Others <- NULL
metals$cluster_time_interaction <- NULL
Diversities$time <- NULL
Diversities <- column_to_rownames(Diversities, var = "Sample")

# Select only numeric variables from microbiota
numeric_microbiota <- x1 %>% select_if(is.numeric)
dim(numeric_microbiota)

# Set up environmental variables (metals)
sel_env <- c("Al", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "As", "Mo", "Cd", "Sb", "Tl", "Pb")
sel_env_label <- list(
  'Al' = "Aluminum",
  'V' = "Vanadium",
  'Cr' = "Chromium",
  'Mn' = "Manganese",
  'Fe' = "Iron",
  'Co' = "Cobalt",
  'Ni' = "Nickel",
  'Cu' = "Copper",
  'Zn' = "Zinc",
  'As' = "Arsenic",
  'Mo' = "Molybdenum",
  'Cd' = "Cadmium",
  'Sb' = "Antimony",
  'Tl' = "Thallium",
  'Pb' = "Lead"
)

# Transpose microbiota data (optional)
microbiota <- t(microbiota)
dim(microbiota)

# Filter out samples with fewer counts (optional)
x <- microbiota[rowSums(microbiota) > 0, ]
dim(x)

# Apply normalization (log-relative transformation)
x <- log((x + 1) / (rowSums(x) + ncol(x)))
dim(x)

# Order columns based on colSums (optional)
x <- x[, order(colSums(x), decreasing = TRUE)]
dim(x)

# Extract top N taxa (optional)
N <- 20  # You can change this to your desired number
taxa_list <- colnames(x)[1:N]
taxa_list <- taxa_list[!grepl("Unknown", taxa_list)]  # Remove "__Unknown__"
taxa_list <- taxa_list[!grepl("Ralstonia", taxa_list)]  # Remove "__Unknown__"
taxa_list <- taxa_list[!grepl("Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", taxa_list)]  # Remove "__Unknown__"
N <- length(taxa_list)

x <- x[, colnames(x) %in% taxa_list]
dim(x)

# Prepare data for correlation analysis
x <- x
y <- metals

# Set rownames as a column name for the microbiota
x$Sample <- rownames_to_column(x, var="Sample")
head(x)

# Merge dataframes
XY <- merge(x, y, by = "row.names", all = TRUE)
XY <- cbind(x, y)
dim(XY)

# Remove NAs (optional)
XY <- na.omit(XY)
dim(XY)

# Split data into x and y
x <- XY[, c(1, 2:9)]  # Select columns corresponding to the 'x' data frame
y <- XY[, c(10, 11:24)]  # Select columns corresponding to the 'y' data frame

# Grouping information
groups <- XY$time  # Use "Time" variable for grouping

# Both must be numeric data frames
x <- x %>% select_if(is.numeric)
y <- y %>% select_if(is.numeric)

x <- as.data.frame(x)
y <- as.data.frame(y)

# Choose correlation method
method <- "spearman"

# Initialize an empty list to store the results
results_list <- list()

for (i in colnames(x)) {
  for (j in colnames(y)) {
    for (k in unique(groups)) {
      a <- as.numeric(x[groups == k, i])
      b <- as.numeric(y[groups == k, j])
      
      if (!any(is.na(a)) && !any(is.na(b)) && length(a) > 0 && length(b) > 0) {
        tmp <- data.frame(
          Taxa = i,
          Metals = j,
          Correlation = cor(a, b, use = "everything", method = method),
          Pvalue = cor.test(a, b, method = method)$p.value,
          Time = k,
          stringsAsFactors = FALSE
        )
        
        # Append to the results list
        results_list <- append(results_list, list(tmp))
      }
    }
  }
}

# Combine all results into a single data frame
df <- do.call(rbind, results_list)

# Convert necessary columns to numeric
df$Pvalue <- as.numeric(df$Pvalue)
df$Correlation <- as.numeric(df$Correlation)

# Adjust p-values for multiple comparisons
df$AdjPvalue <- p.adjust(df$Pvalue, method = "BH")

# Generate labels for significant values
df$Significance <- cut(df$Pvalue, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("***", "**", "*", ""))

# Display the final data frame
print(df)

# Reorganize Metals and bacteria data based on concentrations and abundances
average_concentrations <- colMeans(metals)
ordered_metals <- names(average_concentrations)[order(-average_concentrations)]
df$Metals <- factor(df$Metals, levels = ordered_metals)

average_abundances <- colMeans(microbiota)
ordered_bacteria <- names(average_abundances)[order(average_abundances)]
df$Taxa <- factor(df$Taxa, levels = ordered_bacteria)

# Create heatmap
Metal_labeller <- function(variable, value) {
  return(label_value(value))  # Using label_value from forcats
}

p <- ggplot(aes(x = Time, y = Taxa, fill = Correlation), data = df) +
  geom_tile() +
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  geom_text(aes(label = Significance), color = "black", size = 3) +
  labs(y = NULL, x = NULL, fill = method) +
  facet_grid(. ~ Metals, scale = "free", space = "free_x", labeller = Metal_labeller)

# Display heatmap
print(p)

# Save plot as a PDF
pdf("Correlation_Metals_diver21_03_2024.pdf", height = 8, width = 22)
print(p)
dev.off()


