# BM_metallomic_16S_heatmap correlations
This repository contains the script utilized to perform the heatmap with metal and metalloids and microbiota abundances including diversity indices.
library(dplyr)
library(tibble)

# Cleaning
more_data$Código_Familia <- NULL
infant_genus$Group <- NULL
infant_genus$Origen <- NULL
infant_genus$Muestra <- NULL
microbiota
#changing names
microbiota <- infant_genus
metals <- BM_minerals
z <- more_data
microbiota <- column_to_rownames(microbiota, "Sample")
microbiota$time <- NULL
# Select only numeric variables from microbiota
numeric_microbiota <- microbiota %>% select_if(is.numeric)

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

# Filter out samples with fewer counts (optional)
microbiota <- microbiota[rowSums(numeric_microbiota) > 0, ]

# Apply normalization (log-relative transformation)
# Apply the log transformation (choose normalization)
x <- log((numeric_microbiota + 1) / (rowSums(numeric_microbiota) + ncol(numeric_microbiota)))

install.packages("compositions") #package for CLR normalization
library(compositions)

x <- clr(microbiota)
# Order columns based on colSums (optional)
x <- x[, order(colSums(x), decreasing = TRUE)]

# Extract top N taxa (optional)
N <- 10  # You can change this to your desired number
taxa_list <- colnames(x)[1:N]
taxa_list <- taxa_list[!grepl("Unknown", taxa_list)]  # Remove "__Unknown__"
N <- length(taxa_list)
x <- data.frame(x[, colnames(x) %in% taxa_list])

# Change names to x and y
x <- Infant_div #microbiota or diversities
dim(x)
x$time <- NULL
y <- metals

dim(x)
dim(y)
dim(z)
#make sure both dataframes have same length (optinal)
#before merging make sure that Sample is in columns
x <- rownames_to_column(x,"Sample")

XY <- merge(x, y, by = "Sample", all = TRUE)
XYZ <- merge(XY, z, by = "Sample", all = TRUE)

dim(XY)
dim(XYZ)
head(XY)
head(XYZ)
table(XYZ$time.x)
#remove Nas (optional)
XYZ <- na.omit(XYZ)

# Print the result
print(XYZ)
view(XYZ)
XY
# Assuming 'Row.names' is the common column after merging
x <- XYZ[, c(1:8)]  # Select columns corresponding to the 'x' data frame
y <- XYZ[, c(1, 9:24)]  # Select columns corresponding to the 'y' data frame
z <- XYZ[, c(1,9,25:28)] #sample data for adjusting p-values for spearman test

dim(x)
str(x)
str(y)
str(z)
view(y)
view(x)
view(z)
# Grouping information
groups <- z$time.y  # Use "Time" variable for grouping
length(groups)

# Both must be numeric data frames
x <- x %>% select_if(is.numeric)
y <- y %>% select_if(is.numeric)

#check dimension of grouping variable, X, and, Y

# Choose correlation method
# Choose correlation method
method <- "spearman"

df <- NULL

for (i in colnames(x)) {
  for (j in colnames(y)) {
    for (k in unique(groups)) {
      a <- as.numeric(x[groups == k, i])
      b <- as.numeric(y[groups == k, j])
      
      if (!any(is.na(a)) && !any(is.na(b)) && length(a) > 0 && length(b) > 0) {
        # Calculate correlation coefficient and p-value
        cor_coef <- cor(a, b, use = "everything", method = method)
        p_value <- cor.test(a, b, method = method)$p.value
        
        # Apply Fisher transformation
        fisher_transform <- 0.5 * log((1 + cor_coef) / (1 - cor_coef))
        
        # Create a data frame with results
        tmp <- data.frame(Variable1 = i,
                          Variable2 = j,
                          Correlation = cor_coef,
                          P_Value = p_value,
                          Group = k,
                          Fisher_Transform = fisher_transform)
        
        if (is.null(df)) {
          df <- tmp
        } else {
          df <- rbind(df, tmp)
        }
      }
    }
  }
}

# Print the resulting data frame
print(df)

warnings()
# Create a data frame
df <- data.frame(row.names = NULL, df)
colnames(df) <- c("Taxa", "Metals", "Correlation", "Pvalue", "Time","Fisher_tranform")
df$Pvalue <- as.numeric(as.character(df$Pvalue))
df$AdjPvalue <- rep(0, dim(df)[1])
df$Correlation <- as.numeric(as.character(df$Correlation))


# Adjust p-values for multiple comparisons
# Use appropriate adjustment method based on your preference
# 1 -> don't adjust
# 2 -> adjust Metals + Type
# 3 -> adjust Taxa + Type
# 4 -> adjust Taxa, Sex, Parto, and Tiempo_de_lactancia
# 5 -> adjust Metals, Sex, Parto, and Tiempo_de_lactancia
adjustment <- 4

if (adjustment %in% c(4, 5)) {
  # Additional adjustments for external variables (Sex, Parto, Lactancia)
  for (k in unique(z$Sex)) {
    for (l in unique(z$Parto))
        sel <- z$Sex == k & z$Parto == l
        df$AdjPvalue[sel] <- p.adjust(df$Pvalue[sel], method = "BH")
      }
    }
df

#adjust only for delivery/parto (optional)
adjustment <- 4

if (adjustment %in% c(4, 5)) {
  # Additional adjustments for external variable (Parto)
  for (l in unique(z$Parto)) {
    sel <- z$Parto == l
    df$AdjPvalue[sel] <- p.adjust(df$Pvalue[sel], method = "BH")
  }
}

# Generate labels for significant values
df$Significance <- cut(df$AdjPvalue, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), label = c("***", "**", "*", ""))

# Remove NAs
df <- df[complete.cases(df), ]
dim(df)
head(df)
# Reorganize Metals data based on appearance (optional)
df$Metals <- factor(df$Metals, as.character(df$Metals))
----
#Reorganize Metals data based on concentrations in milk
# Calculate average concentrations for each metal
num_metal <- metals
num_metal$time <- NULL
num_metal$Sample <- NULL
average_concentrations <- colMeans(num_metal)
average_concentrations
# Order Metals based on average concentrations
ordered_metals <- names(average_concentrations)[order(-average_concentrations)]

# Reorganize Metals data based on the new order
df$Metals <- factor(df$Metals, levels = ordered_metals)

----
  #Reorganize bacteria data based on abundances on infant feaces
  # Calculate average abundances for each bacteria

average_concentrations <- colMeans(x)
average_concentrations
# Order taxa based on average concentrations
ordered_metals <- names(average_concentrations)[order(average_concentrations)]
---

# Reorganize bacteria data based on the new order
df$Taxa <- factor(df$Taxa, levels = ordered_metals)
-----
library(forcats)
# Updated labeller function
Metal_labeller <- function(variable, value) {
  return(label_value(value))  # Using label_value from forcats
}
library(ggplot2)
# Create ggplot heatmap (if there are Nas go back to remove NAs in df)
p <- ggplot(aes(x = Time, y = Taxa, fill = Correlation), data = df)
p <- p + geom_tile() + scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C") 
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p <- p + geom_text(aes(label = Significance), color = "black", size = 3) + labs(y = NULL, x = NULL, fill = method)
p <- p + facet_grid(. ~ Metals, drop = TRUE, scale = "free", space = "free_x", labeller = Metal_labeller)

#change name
Adjusted_heatmap_diversityCLRadj <- p
Adjusted_heatmap_deliverymodeCLR <- p

Heatmap_apha_indices <- p
# Save plot as a PDF
pdf("Correlation_Metals_Microbiota.pdf", height = 8, width = 22)
print(p)
dev.off()

#to keep
print(Adjusted_heatmap_deliverymodeCLR)
print(Adjusted_heatmap_diversityCLRadj)
#testing
print(Adjusted_heatmap_deliverymode)
print(Adjusted_heatmap_sexpartlact)
print(Heatmap_early_mature)
print(Heatmap_apha_indices)

# Assuming df is your dataframe this is a faceted barplot (another way to visualize the correlation values)
faceted_barplot_bacteria <- ggplot(df, aes(x = Metals, y = Correlation, fill = Time)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Taxa, scales = "free", nrow = 4, ncol = 5) +
  labs(x = "Metals", y = "Correlation", fill = "Time") +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),  # Adjust the angle as needed
    plot.title = element_text(size = 12),  # Adjust the title font size
    legend.text = element_text(size = 8),  # Adjust the legend font size
    strip.text = element_text(size = 8)  # Adjust the Taxa names font size
  )

print(faceted_barplot_indices)
print(faceted_barplot_bacteria)


