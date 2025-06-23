##########################################################################
# Script Name:    negativebin_model.R
# Author:         Maite Benlloch Montón
# Date:           [2025-06-23]
# Description:    Estimates dispersion from real RNA-seq count data, simulates 
#                 biological replicates using a Negative Binomial model, and 
#                 compares the original vs simulated distributions using 
#                 density plots and a Kolmogorov–Smirnov test.
# 'edgeR' approach
#
# Input:          Count matrices (genes x samples) from RNA-seq datasets.
# Output:         Density plots comparing original and simulated data,
#                 summary of KS test statistics (D and p-value).
#
# Dependencies:   edgeR, tidyverse, ggplot2, dplyr, tidyr
#
# Usage Example:  
#   result <- analyze_dataset(counts_matrix, "Title of Dataset")
#   print(result$plot)
#   cat(result$conclusion)
#
# Notes:
# - Assumes gene/transcript names are row names of the input count matrix.
# - Simulation is based on tagwise dispersion estimated with edgeR.
# - Set a seed for reproducibility of simulations.
#
##########################################################################

#=============================================================================================
# FUNCTION TO GENERATE SIMULATED REPLICATES
#=============================================================================================
library(edgeR)
library(tidyverse)

# Function to generate Negative Binomial replicates given a mean and dispersion
generate_replicates_nb <- function(mu, dispersion, numReps = 5) {
  t(sapply(seq_along(mu), function(i) {
    rnbinom(numReps, mu = mu[i], size = 1 / dispersion[i])
  }))
}

# Function to estimate parameters and generate biological replicates
generate_biological_replicates <- function(counts) {
  library(edgeR)
  
  # Step 1: Calculate percentage of zeros per column (sample)
  zero_percentages <- colMeans(counts == 0)
  min_sparsity <- min(zero_percentages)
  max_sparsity <- max(zero_percentages)
  
  # Step 2: Estimate NB parameters from real data
  group <- rep(1, ncol(counts))
  dge <- DGEList(counts = counts, group = group)
  dge <- normLibSizes(dge, method = "TMMwsp")
  dge <- estimateDisp(dge, robust = TRUE)
  
  mu <- rowMeans(dge$counts)
  dispersion <- dge$tagwise.dispersion
  
  num_replicates <- ncol(counts)
  set.seed(123)
  replicates <- generate_replicates_nb(mu, dispersion, numReps = num_replicates)
  colnames(replicates) <- colnames(counts)
  rownames(replicates) <- rownames(counts)
  
  # Step 3: Introduce dropout adjusted to real sparsity
  for (j in seq_len(ncol(replicates))) {
    target_sparsity <- runif(1, min = min_sparsity, max = max_sparsity)
    current_zero_count <- sum(replicates[, j] == 0)
    target_zero_count <- floor(nrow(replicates) * target_sparsity)
    zeros_needed <- target_zero_count - current_zero_count
    
    if (zeros_needed > 0) {
      candidates <- which(replicates[, j] > 0)
      genes_to_zero <- sample(candidates, min(zeros_needed, length(candidates)), replace = FALSE)
      replicates[genes_to_zero, j] <- 0
    }
  }
  
  return(replicates)
}

# Generate simulated replicates for each dataset
rep_k_ont <- generate_biological_replicates(counts_k_ont)
rep_k_isoseq <- generate_biological_replicates(counts_k_isoseq)
rep_b_ont <- generate_biological_replicates(counts_b_ont)
rep_b_isoseq <- generate_biological_replicates(counts_b_isoseq)

#=============================================================================================
# FUNCTION TO ANALYZE REAL VS SIMULATED DISTRIBUTIONS
#=============================================================================================
analyze_dataset <- function(original_counts, simulated_replicates, title, alpha = 0.05) {
  library(tidyverse)
  library(ggplot2)
  
  # Transform original and simulated counts: log10(count + 1)
  df_original <- as.data.frame(log10(original_counts + 1)) %>%
    rownames_to_column("Transcript") %>%
    pivot_longer(-Transcript, names_to = "sample", values_to = "value") %>%
    mutate(type = "Original")
  
  df_replicates <- as.data.frame(log10(simulated_replicates + 1)) %>%
    rownames_to_column("Transcript") %>%
    pivot_longer(-Transcript, names_to = "sample", values_to = "value") %>%
    mutate(type = "Simulated biological replicates")
  
  df_combined <- bind_rows(df_original, df_replicates)
  
  # Density plot
  plot <- ggplot(df_combined, aes(x = value, fill = type, color = type)) +
    geom_density(alpha = 0.3) +
    theme_minimal(base_size = 14) +
    labs(title = title, x = "log10(count + 1)", y = "Density") +
    scale_fill_manual(values = c("red", "skyblue")) +
    scale_color_manual(values = c("red", "skyblue"))
  
  # Kolmogorov–Smirnov test
  ks <- ks.test(df_original$value, df_replicates$value)
  
  result <- list(
    title = title,
    D = ks$statistic,
    p_value = ks$p.value,
    conclusion = if (ks$p.value >= alpha) {
      "Similar distributions (fail to reject H0)"
    } else {
      "Different distributions (reject H0)"
    },
    plot = plot
  )
  
  return(result)
}

# Apply analysis to each dataset
res_k_ont_ana     <- analyze_dataset(counts_k_ont,     rep_k_ont,     "Kidney ONT: Original vs Simulated")
res_k_isoseq_ana  <- analyze_dataset(counts_k_isoseq,  rep_k_isoseq,  "Kidney IsoSeq: Original vs Simulated")
res_b_ont_ana     <- analyze_dataset(counts_b_ont,     rep_b_ont,     "Brain ONT: Original vs Simulated")
res_b_isoseq_ana  <- analyze_dataset(counts_b_isoseq,  rep_b_isoseq,  "Brain IsoSeq: Original vs Simulated")

# Show plots
print(res_k_ont_ana$plot)
print(res_k_isoseq_ana$plot)
print(res_b_ont_ana$plot)
print(res_b_isoseq_ana$plot)

# Save plots to folder
plots_list <- list(
  res_k_ont_plot = res_k_ont_ana$plot,
  res_k_isoseq_plot = res_k_isoseq_ana$plot,
  res_b_ont_plot = res_b_ont_ana$plot,
  res_b_isoseq_plot = res_b_isoseq_ana$plot
)

output_dir <- "C:/Users/mbenm/Downloads/TFG/scripts_def/Visualizations/model"

for (name in names(plots_list)) {
  ggsave(
    filename = file.path(output_dir, paste0(name, "_plot.png")),
    plot = plots_list[[name]],
    width = 8,
    height = 6,
    units = "in",
    dpi = 300
  )
}

# Kolmogorov–Smirnov Test Results Summary
cat("Kolmogorov–Smirnov Test Results:\n\n")
for (res in list(res_k_ont_ana, res_k_isoseq_ana, res_b_ont_ana, res_b_isoseq_ana)) {
  cat("→", res$title, "\n")
  cat("  D =", round(res$D, 4), " | p-value =", signif(res$p_value, 4), "\n")
  cat("  Conclusion:", res$conclusion, "\n\n")
}

# Alternative analyze_dataset function using vectorized log transform and plot
analyze_dataset <- function(original_counts, simulated_replicates, title, alpha = 0.05) {
  library(tidyverse)
  library(ggplot2)
  
  # Transform data to log10(count + 1) vectors
  original_log <- log10(as.vector(as.matrix(original_counts)) + 1)
  simulated_log <- log10(as.vector(as.matrix(simulated_replicates)) + 1)
  
  # Density plot without combining dataframes
  plot <- ggplot() +
    geom_density(aes(x = original_log), color = "red", fill = "red", alpha = 0.3) +
    geom_density(aes(x = simulated_log), color = "skyblue", fill = "skyblue", alpha = 0.3) +
    theme_minimal(base_size = 14) +
    labs(
      title = title,
      x = "log10(count + 1)",
      y = "Density"
    ) +
    annotate("text", x = Inf, y = Inf, label = "Red: Original\nBlue: Simulated",
             hjust = 1.1, vjust = 1.5, size = 4, color = "black")
  
  # Kolmogorov–Smirnov test
  ks <- ks.test(original_log, simulated_log)
  
  result <- list(
    title = title,
    D = ks$statistic,
    p_value = ks$p.value,
    conclusion = if (ks$p.value >= alpha) {
      "Similar distributions (fail to reject H0)"
    } else {
      "Different distributions (reject H0)"
    },
    plot = plot
  )
  
  return(result)
}

#=============================================================================================
# FUNCTION TO PLOT DENSITY FOR EACH REPLICATE PAIR INDIVIDUALLY
#=============================================================================================
plot_density_pairs <- function(original_counts, simulated_counts, dataset_name = "Dataset") {
  # Assign generic names if missing
  if (is.null(colnames(original_counts))) {
    colnames(original_counts) <- paste0("Sample", seq_len(ncol(original_counts)))
  }
  
  # Ensure simulated counts have the same column names if missing
  if (is.null(colnames(simulated_counts))) {
    colnames(simulated_counts) <- colnames(original_counts)
  }
  
  # Convert original counts to long format with log10(count + 1)
  df_original <- as.data.frame(log10(original_counts + 1)) %>%
    mutate(Transcript = rownames(original_counts)) %>%
    pivot_longer(-Transcript, names_to = "Sample", values_to = "Value") %>%
    mutate(Type = "Original")
  
  # Convert simulated counts similarly
  df_simulated <- as.data.frame(log10(simulated_counts + 1)) %>%
    mutate(Transcript = rownames(simulated_counts)) %>%
    pivot_longer(-Transcript, names_to = "Sample", values_to = "Value") %>%
    mutate(Type = "Simulated")
  
  # Combine both datasets
  df_combined <- bind_rows(df_original, df_simulated)
  
  # Plot with facets per sample
  ggplot(df_combined, aes(x = Value, fill = Type, color = Type)) +
    geom_density(alpha = 0.3) +
    facet_wrap(~ Sample, ncol = 3) +
    scale_fill_manual(values = c("Original" = "red", "Simulated" = "skyblue")) +
    scale_color_manual(values = c("Original" = "red", "Simulated" = "skyblue")) +
    theme_minimal(base_size = 13) +
    labs(
      title = paste0(dataset_name, ": Comparison of Original Biological Replicates vs Simulated"),
      x = "log10(count + 1)",
      y = "Density"
    )
}

# Kidney - ONT
colnames(counts_k_ont) <- c("K31", "K32", "K33", "K34", "K35")
colnames(rep_k_ont) <- colnames(counts_k_ont)
plot_density_pairs(counts_k_ont, rep_k_ont, dataset_name = "Kidney - ONT")

# Kidney - IsoSeq
colnames(counts_k_isoseq) <- c("K31", "K32", "K33", "K34", "K35")  # Adjust as needed
colnames(rep_k_isoseq) <- colnames(counts_k_isoseq)
plot_density_pairs(counts_k_isoseq, rep_k_isoseq, dataset_name = "Kidney - IsoSeq")

# Brain - ONT
colnames(counts_b_ont) <- c("B31", "B32", "B33", "B34", "B35")  # Adjust as needed
colnames(rep_b_ont) <- colnames(counts_b_ont)
plot_density_pairs(counts_b_ont, rep_b_ont, dataset_name = "Brain - ONT")

# Brain - IsoSeq
colnames(counts_b_isoseq) <- c("B31", "B32", "B33", "B34", "B35")  # Adjust as needed
colnames(rep_b_isoseq) <- colnames(counts_b_isoseq)
plot_density_pairs(counts_b_isoseq, rep_b_isoseq, dataset_name = "Brain - IsoSeq")

#=============================================================================================
# VIOLIN PLOTS FOR SIMULATED DATA
#=============================================================================================
library(ggplot2)
library(tidyr)
library(dplyr)

plot_violin_simulated <- function(simulated_counts, dataset_name = "Dataset") {
  # Ensure column names exist
  if (is.null(colnames(simulated_counts))) {
    colnames(simulated_counts) <- paste0("Sample", seq_len(ncol(simulated_counts)))
  }
  
  # Convert to long format with log10(count + 1)
  df_simulated <- as.data.frame(log10(simulated_counts + 1)) %>%
    mutate(Transcript = rownames(simulated_counts)) %>%
    pivot_longer(-Transcript, names_to = "Sample", values_to = "Value")
  
  # Violin plot
  ggplot(df_simulated, aes(x = Sample, y = Value, fill = Sample)) +
    geom_violin(trim = FALSE, alpha = 0.6) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5, position = position_dodge(width = 0.9)) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste0(dataset_name, ": Simulated Biological Replicates Distributions"),
      x = "Biological Replicates",
      y = "log10(count + 1)"
    ) +
    theme(legend.position = "none")  # Hide legend since violins are identified by X axis
}

plot_violin_simulated(rep_k_ont, dataset_name = "Kidney - ONT")
plot_violin_simulated(rep_k_isoseq, dataset_name = "Kidney - IsoSeq")
plot_violin_simulated(rep_b_ont, dataset_name = "Brain - ONT")
plot_violin_simulated(rep_b_isoseq, dataset_name = "Brain - IsoSeq")

#=============================================================================================
# VIOLIN PLOTS FOR ORIGINAL DATA
#=============================================================================================
plot_violin_original <- function(original_counts, dataset_name = "Dataset") {
  # Ensure column names exist
  if (is.null(colnames(original_counts))) {
    colnames(original_counts) <- paste0("Sample", seq_len(ncol(original_counts)))
  }
  
  # Convert to long format with log10(count + 1)
  df_original <- as.data.frame(log10(original_counts + 1)) %>%
    mutate(Transcript = rownames(original_counts)) %>%
    pivot_longer(-Transcript, names_to = "Sample", values_to = "
  
  # Violin plot
  ggplot(df_original, aes(x = Sample, y = Value, fill = Sample)) +
    geom_violin(trim = FALSE, alpha = 0.6) +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5, position = position_dodge(width = 0.9)) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste0(dataset_name, ": Original Biological Replicates Distributions"),
      x = "Biological Replicates",
      y = "log10(count + 1)"
    ) +
    theme(legend.position = "none")
}

plot_violin_original(counts_k_ont, dataset_name = "Kidney - ONT")
plot_violin_original(counts_k_isoseq, dataset_name = "Kidney - IsoSeq")
plot_violin_original(counts_b_ont, dataset_name = "Brain - ONT")
plot_violin_original(counts_b_isoseq, dataset_name = "Brain - IsoSeq")
