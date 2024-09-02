# Clear the environment to remove any existing objects
rm(list = ls())

# Load necessary libraries for data manipulation, visualization, and analysis
library(dplyr)
library(ggplot2)
library(tidyverse)
library(purrr)
library(pheatmap)
library(heatmaply)
library(DT)

# Load and preprocess the data
data <- read.csv("/Users/lazizasamov/Desktop/Borrelia_RNA-SEQ/Borrelia_data 2.csv")
data <- data %>% select(-X)  # Remove the unnecessary 'X' column generated during data import
data[, 3:11] <- data[, 3:11] + 1  # Add 1 to counts to avoid issues during log transformation

# Convert data to long format for easier manipulation and analysis
data_long <- data %>%
  gather(stage, counts, Early.1:Stationary.3, factor_key = TRUE) %>%
  mutate(
    Sample = case_when(
      stage == 'Early.1' ~ 'Early Sample 1',
      stage == 'Early.2' ~ 'Early Sample 2',
      stage == 'Early.3' ~ 'Early Sample 3',
      stage == 'Mid.1' ~ 'Mid Sample 1',
      stage == 'Mid.2' ~ 'Mid Sample 2',
      stage == 'Mid.3' ~ 'Mid Sample 3',
      stage == 'Stationary.1' ~ 'Stationary Sample 1',
      stage == 'Stationary.2' ~ 'Stationary Sample 2',
      stage == 'Stationary.3' ~ 'Stationary Sample 3'
    ),
    stage = factor(stage, 
                   levels = c("Early.1", "Early.2", "Early.3", 
                              "Mid.1", "Mid.2", "Mid.3", 
                              "Stationary.1", "Stationary.2", "Stationary.3"),
                   labels = c("Early", "Early", "Early", 
                              "Mid", "Mid", "Mid", 
                              "Stationary", "Stationary", "Stationary"))
  )

# Add a log-transformed counts column for better scaling in visualization and analysis
data_long <- data_long %>%
  mutate(logCounts = log10(counts))

# Calculate mean and standard deviation for a specific gene across stages (Test for the app)
df_mean_std <- data_long %>%
  filter(Gene == "BB_0001_-_hypothetical_protein") %>%
  group_by(stage) %>%
  summarise(mean = mean(logCounts), sd = sd(logCounts)) %>%
  as.data.frame()

# Plot the mean and standard deviation along with individual logCounts for visualization
ggplot(df_mean_std, aes(x = stage, y = mean)) +
  geom_point(size = 4, shape = 24, colour = "lightblue", fill = "darkturquoise") +
  geom_point(data = data_long %>% filter(Gene == 'BB_0001_-_hypothetical_protein'),
             aes(x = stage, y = logCounts), size = 3, shape = 21, fill = 'deepskyblue3', colour = 'deepskyblue3') +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), size = 0.9, width = 0.4, linetype = 2, colour = 'red') +
  theme_bw()

# Add F and P values columns to the dataset for statistical analysis
data_long <- data_long %>%
  mutate(F_value = NA, P_value = NA)

# Calculate F and P values using ANOVA for each gene and collect significant results
p_f_values <- data.frame(Gene = factor(), F_value = numeric(), P_value = numeric())

for (i in unique(data_long$Gene)) {
  model <- lm(counts ~ stage, data = data_long %>% filter(Gene == i))
  anv <- anova(model)
  if (anv[1, 5] < 0.05) {  # Filter for statistically significant results (p < 0.05)
    df <- data.frame(Gene = factor(i), F_value = anv[1, 4], P_value = anv[1, 5])
    p_f_values <- bind_rows(p_f_values, df)
  }
}

# Merge statistical results with the long-format data
data_stat <- merge(p_f_values, data_long, by = 'Gene') %>%
  select(-F_value.y, -P_value.y) %>%
  rename(F_values = F_value.x, P_values = P_value.x)

# Reshape data for heatmap visualization using pivot_wider
heatmap_df <- data_stat %>%
  select(Gene, Sample, logCounts, F_values, P_values) %>%
  pivot_wider(names_from = Sample, values_from = logCounts)

# Convert the resulting tibble to a data frame, set row names to Gene, and drop the Gene column
heatmap_df <- as.data.frame(heatmap_df)
rownames(heatmap_df) <- heatmap_df$Gene
heatmap_df <- heatmap_df %>% select(-Gene)

# Order columns to match desired heatmap structure
col_order <- c('Early Sample 1','Early Sample 2','Early Sample 3', 
               'Mid Sample 1', 'Mid Sample 2', 'Mid Sample 3', 
               'Stationary Sample 1', 'Stationary Sample 2', 'Stationary Sample 3', 
               'P_values', 'F_values')
heatmap_df <- heatmap_df[, col_order]

# Plot the heatmap using heatmaply, with the first 9 columns (samples) as input
heatmaply(heatmap_df[1:9], showticklabels = FALSE)

# Testing ANOVA model for a specific gene (Test for the app))
model <- lm(counts ~ stage, data = data_long %>% filter(Gene == 'BB_0007_-_hypothetical_protein'))
dat <- anova(model)

# Display the ANOVA results using a DataTable
datatable(dat, options = list(dom = 't'))



