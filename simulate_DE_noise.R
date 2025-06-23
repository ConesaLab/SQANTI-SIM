##########################################################################
# Script Name:    simulate_DE_noise.R
# Author:         Maite Benlloch Montón
# Date:           [2025-06-23]
# Description:    Simulates gene expression data with specified fold changes 
#                 and adds biological noise. Generates plots to visualize the
#                 effect of noise and expression patterns across conditions.
#
# 'edgeR' approach
# Input:          Raw count matrix (`counts`), assumed to be transcripts x samples.
# Output:         List containing simulated counts, noisy counts, fold changes,
#                 and ggplot2 density plots.
#
# Dependencies:   edgeR, tidyverse, reshape2
#
# Usage Example:  
#   result <- simulate_and_visualize_with_noise(counts_matrix, "Sample Title")
#   print(result$facet_plot)
#   print(result$noise_plot)
#
# Notes:
# - Assumes gene identifiers are row names of the input matrix.
# - Adjust fold change parameters (pct_up, pct_down, minFC, maxFC) as needed.
#
##########################################################################

simulate_and_visualize_with_noise <- function(counts, title,
                                              num_replicates = 2,
                                              pct_up = 0.30, pct_down = 0.10,
                                              minFC = 0.25, maxFC = 4,
                                              noise = 0.5) {
  library(edgeR)
  library(tidyverse)
  library(reshape2)
  
  # --- Estimate NB parameters from real data ---
  group <- rep(1, ncol(counts))
  dge <- DGEList(counts = counts, group = group)
  dge <- normLibSizes(dge, method = "TMMwsp")
  dge <- estimateDisp(dge, robust = TRUE)
  
  mu_base <- rowMeans(dge$counts)
  dispersion <- dge$tagwise.dispersion
  n_genes <- length(mu_base)
  gene_ids <- rownames(counts)
  
  # --- Calculate actual sparsity ---
  zero_percentages <- colMeans(counts == 0)
  min_sparsity <- min(zero_percentages)
  max_sparsity <- max(zero_percentages)
  
  # --- Helper function to simulate replicates with dropout ---
  generate_biological_replicates_DE_noise <- function(mu, dispersion, num_reps, min_sp, max_sp) {
    replicates <- t(sapply(seq_along(mu), function(i) {
      rnbinom(num_reps, mu = mu[i], size = 1 / dispersion[i])
    }))
    colnames(replicates) <- paste0("Sample_", 1:num_reps)
    rownames(replicates) <- gene_ids
    
    for (j in seq_len(ncol(replicates))) {
      target_sparsity <- runif(1, min = min_sp, max = max_sp)
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
  
  # --- Simulate control group ---
  set.seed(123)
  control_replicates <- generate_biological_replicates_DE_noise(mu_base, dispersion, num_replicates, min_sparsity, max_sparsity)
  colnames(control_replicates) <- paste0("Control_", 1:num_replicates)
  
  # --- Simulate tumor group with fold changes ---
  set.seed(42)
  group_assignment <- sample(
    c(rep("up", round(pct_up * n_genes)),
      rep("down", round(pct_down * n_genes)),
      rep("normal", n_genes - round(pct_up * n_genes) - round(pct_down * n_genes)))
  )
  names(group_assignment) <- gene_ids
  
  df_fc <- data.frame(Gene = gene_ids, Group = group_assignment, FC = 1)
  df_fc$FC[df_fc$Group == "down"] <- runif(sum(df_fc$Group == "down"), min = 0.0001, max = minFC)
  df_fc$FC[df_fc$Group == "up"]   <- runif(sum(df_fc$Group == "up"), min = maxFC, max = 100)
  
  tumor_mu <- mu_base * df_fc$FC
  tumor_mu[tumor_mu < 0] <- 0
  
  tumor_replicates <- generate_biological_replicates_DE_noise(tumor_mu, dispersion, num_replicates, min_sparsity, max_sparsity)
  colnames(tumor_replicates) <- paste0("Tumor_", 1:num_replicates)
  
  # --- Final simulated matrix ---
  simulated_counts <- cbind(control_replicates, tumor_replicates)
  
  # --- Add biological noise ---
  set.seed(123)
  noisy_counts <- simulated_counts + rnorm(length(simulated_counts), mean = 0, sd = noise)
  dim(noisy_counts) <- dim(simulated_counts)
  colnames(noisy_counts) <- colnames(simulated_counts)
  rownames(noisy_counts) <- rownames(simulated_counts)
  noisy_counts[noisy_counts < 0] <- 0
  noisy_counts <- round(noisy_counts)
  
  # --- Plot 1: Per-replicate distribution (with noise) ---
  counts_df <- as.data.frame(noisy_counts)
  counts_df$Gene <- rownames(counts_df)
  counts_long <- melt(counts_df, id.vars = "Gene", variable.name = "Sample", value.name = "Count")
  counts_long$Group <- ifelse(grepl("Control", counts_long$Sample), "Control", "Tumor")
  counts_long$logCount <- log10(counts_long$Count + 1)
  
  facet_plot <- ggplot(counts_long, aes(x = logCount, color = Group, fill = Group)) +
    geom_density(alpha = 0.3) +
    facet_wrap(~Sample, ncol = 4) +
    theme_minimal() +
    labs(
      title = paste0("Control vs Tumor distribution: ", title),
      x = "log10(counts + 1)", y = "Density"
    )
  
  # --- Function to visualize noise effect in control replicates ---
  plot_noise_effect_facet <- function(sim_counts, noisy_counts, sample_name, n_reps = 2) {
    sim_controls <- sim_counts[, 1:n_reps, drop = FALSE]
    noisy_controls <- noisy_counts[, 1:n_reps, drop = FALSE]
    
    df_sim <- as.data.frame(sim_controls)
    df_noisy <- as.data.frame(noisy_controls)
    df_sim$Gene <- rownames(sim_controls)
    df_noisy$Gene <- rownames(noisy_controls)
    
    sim_long <- melt(df_sim, id.vars = "Gene", variable.name = "Sample", value.name = "Count")
    noisy_long <- melt(df_noisy, id.vars = "Gene", variable.name = "Sample", value.name = "Count")
    
    sim_long$Sample <- paste0(sim_long$Sample, "_NoNoise")
    noisy_long$Sample <- paste0(noisy_long$Sample, "_Noise")
    
    combined <- rbind(sim_long, noisy_long)
    combined$logCount <- log10(combined$Count + 1)
    
    ggplot(combined, aes(x = logCount, color = Sample, fill = Sample)) +
      geom_density(alpha = 0.3) +
      facet_wrap(~Sample, ncol = 4) +
      theme_minimal() +
      labs(
        title = paste0("Comparison of control replicates without and with noise: ", sample_name),
        x = "log10(count + 1)",
        y = "Density"
      ) +
      theme(legend.position = "none")
  }
  
  # Pass the sample title to the function
  noise_plot <- plot_noise_effect_facet(control_replicates, noisy_counts[, 1:num_replicates], title)
  
  return(list(
    facet_plot = facet_plot,
    noise_plot = noise_plot,
    simulated_counts = simulated_counts,
    noisy_counts = noisy_counts,
    fold_changes = df_fc
  ))
}


# Simulate and visualize for "Kidney ONT"
res_kidney_DE_noise <- simulate_and_visualize_with_noise(counts_k_ont, "Kidney ONT")

# View the two plots:
print(res_kidney_DE_noise$facet_plot)
print(res_kidney_DE_noise$noise_plot)

# Simulate and visualize for "Brain ONT"
res_brain_DE_noise <- simulate_and_visualize_with_noise(counts_b_ont, "Brain ONT")
print(res_brain_DE_noise$facet_plot)
print(res_brain_DE_noise$noise_plot)

# Simulate and visualize for "Kidney IsoSeq"
res_kidney_isoseq_DE_noise <- simulate_and_visualize_with_noise(counts_k_isoseq, "Kidney IsoSeq")
print(res_kidney_isoseq_DE_noise$facet_plot)
print(res_kidney_isoseq_DE_noise$noise_plot)

# Simulate and visualize for "Brain IsoSeq"
res_brain_isoseq_DE_noise <- simulate_and_visualize_with_noise(counts_b_isoseq, "Brain IsoSeq")
print(res_brain_isoseq_DE_noise$facet_plot)
print(res_brain_isoseq_DE_noise$noise_plot)
