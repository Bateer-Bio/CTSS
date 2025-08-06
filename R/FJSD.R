#' Calculate the Value at Two-Thirds Position in Sorted Vector
#' @param Numeric vector of input values
#' @return Numeric value at the index = length(numbers) * 2/3 after descending sort
#' @export

two_thirds<- function(numbers){
  sorted_numbers <- sort(numbers, decreasing = TRUE)
  top_two_thirds <- sorted_numbers[length(sorted_numbers) * 2/3]
  return(top_two_thirds)
}

#' Compute FoldChange Jensen-Shannon Divergence (FJSD) Score
#'
#' Calculates gene-specific FJSD scores by combining expression fold-change with
#' Jensen-Shannon Divergence from an idealized distribution. Handles zero values by
#' substitution with 1e-10 prior to normalization.
#'
#' @param Exp Numeric matrix of expression values (genes in columns, samples in rows)
#' @param n0 Integer count of target cell-type samples
#' @param n1 Integer count of background samples
#' @param celltype Character string naming the target cell type
#' @param FC_set Method for fold-change calculation: "median", "mean",
#'              "two_thirds" (default), or "no_replicate"
#'
#' @return A list containing the following elements:
#' - FJSD_matrix: Matrix of FJSD scores (celltype x genes)
#' - FJSD_df: Dataframe of results including gene names, FJSD scores, 1-JSD values,
#'            and log fold-change (sorted by FJSD)
#' @export

FJSD_score <- function(Exp, n0, n1, celltype, FC_set = "two_thirds"){
  Exp[Exp == 0] <- 1e-10
  Exp_norm <- data.frame(apply(Exp, 2, function(x) x/sum(x)), row.names = rownames(Exp))
  #idealized vetcor
  v <- c(rep(1/n0, n0), rep(1e-10, n1))
  #FJSD score initialization
  JSD_score <- c(rep(0, ncol(Exp)))
  FJSD_score <- c(rep(0, ncol(Exp)))
  logFC <- c(rep(0, ncol(Exp)))
  #FJSD score calculate
  for (i in 1:ncol(Exp_norm)) {
    P <- Exp_norm[, i]
    Q <- v
    M <- 0.5*(P+Q)
    #K-L
    KL_P_M <- sum(P * log2(P / M))
    KL_Q_M <- sum(Q * log2(Q / M))
    JSD_P_Q <- 0.5*KL_P_M+0.5*KL_Q_M
    JSD_score[i] <- 1-JSD_P_Q
    if (FC_set == "median") {
      FC <- max(log2((median(Exp[1:n0, i])+1)/(median(Exp[(n0 + 1):nrow(Exp), i])+1)),0)
    }
    if (FC_set == "mean") {
      FC <- max(log2((mean(Exp[1:n0, i])+1)/(mean(Exp[(n0 + 1):nrow(Exp), i])+1)),0)
    }
    if (FC_set == "two_thirds") {
      FC <- max(log2((two_thirds(Exp[1:n0, i])+1)/(two_thirds(Exp[(n0 + 1):nrow(Exp), i])+1)), 0)
    }
    if (FC_set == "no_replicate") {
      FC <- max(log2((Exp[1:n0, i]+1)/mean(Exp[(n0 + 1):nrow(Exp), i]+1)), 0)
    }
    FJSD_score[i] <- (FC*(1-JSD_P_Q))^(1/2)
    logFC[i] <- FC
  }
  FJSD_matrix <- data.frame(FJSD = t(FJSD_score))
  rownames(FJSD_matrix) <- celltype
  colnames(FJSD_matrix) <- colnames(Exp)
  FJSD_df <- data.frame(gene = colnames(Exp), FJSD = FJSD_score, One_Minus_JSD = JSD_score, logFC = logFC)
  FJSD_df <- FJSD_df[order(FJSD_df$FJSD, decreasing = TRUE), ]
  FJSD_list <- list(FJSD_matrix = FJSD_matrix, FJSD_df = FJSD_df)
  return(FJSD_list = FJSD_list)
}

#' Perform Cell-Type-Specific FJSD Analysis
#'
#' Computes FJSD scores for all cell types against background, incorporates significance testing
#' through random sampling simulations, and adjusts p-values via Benjamini-Hochberg.
#'
#' @param Exp_raw Raw expression matrix (genes in rows, samples in columns)
#' @param cell_type Dataframe with two columns: "sample" (sample IDs) and "celltype" (annotations)
#' @param FC_set Fold-change calculation method (see FJSD_score)
#' @param n_simulations Number of random background simulations (default=50) for p-value calculation
#'
#' @return A list containing the following elements:
#' - FJSD_list: List per cell type with detailed results (FJSD_df, p-values, adjusted p-values)
#' - FJSD_matrix: Aggregated FJSD score matrix (celltypes x genes)
#' @export

FJSD <- function(Exp_raw, cell_type, FC_set = "two_thirds", n_simulations = 50){
  Exp_raw <- t(Exp_raw)
  FJSD_list <- list()
  FJSD_df <- data.frame()
  p_values <- data.frame(matrix(ncol = 2, nrow = ncol(Exp_raw)))
  colnames(p_values) <- c("gene","p_val")
  rownames(p_values) <- colnames(Exp_raw)
  p_values$gene <- colnames(Exp_raw)
  BH_values <- data.frame(matrix(ncol = 2, nrow = ncol(Exp_raw)))
  colnames(BH_values) <- c("gene","p_val_adj")
  rownames(BH_values) <- colnames(Exp_raw)
  BH_values$gene <- colnames(Exp_raw)
  flag <- 1
  celltype_all <- unique(cell_type$celltype)
  celltype_num <- length(celltype_all)
  for (i in 1:celltype_num) {
    celltype_i_sample <- cell_type$sample[cell_type$celltype == celltype_all[i]]
    Exp0 <- Exp_raw[rownames(Exp_raw) %in% celltype_i_sample, ]
    Exp1 <- Exp_raw[!rownames(Exp_raw) %in% celltype_i_sample, ]
    Exp <- rbind(Exp0, Exp1)
    n0 <- nrow(Exp0)
    n1 <- nrow(Exp1)
    celltype <- celltype_all[i]
    FJSD_i <- FJSD_score(Exp, n0, n1, celltype, FC_set)
    FJSD_list[[i]] <- FJSD_i$FJSD_df
    names(FJSD_list)[i] <- celltype
    if (flag == 1) {
      FJSD_matrix <- FJSD_i$FJSD_matrix
    } else {
      FJSD_matrix <- rbind(FJSD_matrix, FJSD_i$FJSD_matrix)
    }
    true_FJSD <- FJSD_i$FJSD_df
    background_FJSD <- data.frame(gene = true_FJSD$gene, true = true_FJSD$FJSD)
    for (j in 1:n_simulations) {
      random_samples <- sample(rownames(Exp), n0)
      random_Exp0 <- Exp[rownames(Exp) %in% random_samples, ]
      random_Exp1 <- Exp[!rownames(Exp) %in% random_samples, ]
      random_Exp <- rbind(random_Exp0, random_Exp1)
      random_FJSD <- FJSD_score(random_Exp, n0, n1, celltype, FC_set)$FJSD_df
      col_name <- paste0(j, "_simulation")
      random_FJSD1 <- setNames(data.frame(gene = random_FJSD$gene, random_FJSD$FJSD), c("gene", col_name))
      background_FJSD <- merge(background_FJSD, random_FJSD1, by = 'gene')
    }
    for (gene in colnames(Exp)) {
      gene_idx <- which(background_FJSD$gene == gene)
      gene_background <- as.numeric(background_FJSD[gene_idx, -c(1, 2)])
      gene_mean <- mean(gene_background)
      gene_sd <- sd(gene_background)
      gene_true_FJSD <- background_FJSD[gene_idx, 2]
      gene_z_score <- (gene_true_FJSD - gene_mean) / (gene_sd + 1e-10)
      gene_p_value <- pnorm(gene_z_score, lower.tail = FALSE)
      p_values[gene, 2] <- gene_p_value
    }
    BH_values[, 2] <- p.adjust(p_values[, 2], method = "BH")
    FJSD_list[[i]] <- suppressMessages(left_join(FJSD_list[[i]], p_values))
    FJSD_list[[i]] <- suppressMessages(left_join(FJSD_list[[i]], BH_values))
    flag <- flag + 1
    print(paste("Cell type:", celltype, "finished!"))
  }
  FJSD <- list(FJSD_list = FJSD_list, FJSD_matrix = FJSD_matrix)
  return(FJSD = FJSD)
}

#' Filter Top FJSD Genes by Statistical Significance
#' @param FJSD_list List of data frames containing FJSD results for each cell type.
#'                  Each data frame must contain columns: 'gene', 'FJSD', 'p_val', 'p_val_adj'.
#' @param Top Maximum number of top genes to return per cell type.
#' @param FJSD_cutoff Minimum FJSD score threshold. Genes with FJSD < cutoff are excluded. Default: 0.
#' @param p_val_cutoff Maximum unadjusted p-value threshold. Default: 0.05.
#' @param p_adj_cutoff Maximum Benjamini-Hochberg adjusted p-value threshold. Default: 0.05.
#'
#' @return A filtered list of data frames mirroring input structure, containing only genes passing all thresholds.
#' @export

Top_FJSD <- function(FJSD_list, Top = 20, FJSD_cutoff = 0, p_val_cutoff = 0.05, p_adj_cutoff = 0.05){
  Top_FJSD_list <- list()
  for (i in 1:length(FJSD_list)) {
    FJSD_df <- FJSD_list[[i]]
    FJSD_df_filt <- FJSD_df %>%
      filter(FJSD >= FJSD_cutoff & p_val < p_val_cutoff & p_val_adj < p_adj_cutoff)
    if (nrow(FJSD_df_filt) < Top) {
      FJSD_df_filt <- FJSD_df_filt
    }else{
      FJSD_df_filt <- FJSD_df_filt[1:Top, ]
    }
    Top_FJSD_list[[i]] <- FJSD_df_filt
    names(Top_FJSD_list)[i] <- names(FJSD_list)[i]
  }
  return(Top_FJSD_list)
}

#' Determine Optimal FJSD Cutoff via Silhouette Coefficient Analysis
#'
#' Computes optimal FJSD score cutoff by maximizing cluster separation (measured by silhouette coefficient)
#' of cell types using genes selected at varying FJSD cutoffs. Evaluates separation quality across cutoff range.
#'
#' @param Exp_raw Raw expression matrix (genes x samples) for silhouette calculation.
#' @param cell_type Dataframe with two columns: "sample" (sample IDs) and "celltype" (annotations)
#' @param FJSD FJSD result object
#' @param start Starting FJSD cutoff value for evaluation range. Default: 0.
#' @param end Ending FJSD cutoff value for evaluation range. Default: 1.8.
#' @param step Increment between consecutive cutoffs. Default: 0.1.
#'
#' @return A list containing the following elements:
#' - Optimal_FJSD: Data frame with optimal cutoff and corresponding silhouette coefficient
#' - Plot: ggplot object showing silhouette coefficient vs. FJSD cutoffs
#' @export

Optimal_FJSD_Cutoff <- function(Exp_raw, cell_type, FJSD, start = 0, end = 1.8, step = 0.1){
  library(ggplot2)
  clusters <- as.character(cell_type$celltype)
  SC_list <- list()
  flag <- 1
  for (i in seq(start, end, by = step)) {
    marker <- Top_FJSD(FJSD$FJSD_list, Top = Inf, FJSD_cutoff = i, p_val_cutoff = 0.05, p_adj_cutoff = 0.05)
    marker_df <- do.call(rbind, marker)
    marker_gene <- unique(marker_df$gene)
    Exp1 <- Exp_raw[marker_gene, ]
    Exp1.matrix <- dist(x = t(Exp1))
    marker_sil <- cluster::silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = Exp1.matrix)
    SC <- as.data.frame(marker_sil[, 3])
    colnames(SC) <- "SC"
    SC_mean <- mean(SC$SC)
    SC_df <- data.frame(FJSD_cutoff = i, SC_mean = SC_mean)
    SC_list[[flag]] <- SC_df
    flag <- flag + 1
  }
  SC_score <- do.call(rbind, SC_list)
  Optimal_FJSD <- SC_score[which.max(SC_score$SC_mean), ]
  Plot <- ggplot(SC_score, aes(x = FJSD_cutoff, y = SC_mean)) +
    theme_classic() +
    geom_line() +
    geom_point() +
    labs(title = "Average Silhouette Coefficient across Different FJSD Cutoff Values", x = "FJSD Cutoff Values", y = "Average Silhouette Coefficient") +
    theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
          plot.title = element_text(hjust = 0.5, size = 20),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 10),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 10)) +
    scale_x_continuous(breaks = seq(0, 3, by = 0.2))
  Optimal_FJSD_Cutoff <- list(Optimal_FJSD = Optimal_FJSD, Plot = Plot)
  return(Optimal_FJSD_Cutoff)
}

#' Volcano Plot Visualization for FJSD Results
#'
#' @param quary_FJSD Data frame of FJSD results for a single cell type with columns:'gene', 'FJSD', 'p_val_adj'
#' @param FJSD_Cutoff Threshold line for FJSD scores. Default: 0.
#' @param title Plot title. Default: "Master TFs Analysis".
#'
#' @return ggplot object of the volcano plot with labeled genes and threshold annotations.
#' @export

FJSD_geom_point <- function(quary_FJSD, FJSD_Cutoff = 0, title = "Master TFs Analysis") {
  library(ggplot2)
  library(ggrepel)
  quary_FJSD$p_val_adj[quary_FJSD$p_val_adj < 1e-30] <- 1e-30
  quary_FJSD$log10_p_val_adj <- -log10(quary_FJSD$p_val_adj)
  significant_genes <- subset(quary_FJSD, log10_p_val_adj > -log10(0.05) & FJSD > FJSD_Cutoff)
  top_log10_genes <- head(significant_genes[order(-significant_genes$log10_p_val_adj), ], 10)
  top_FJSD_genes <- head(quary_FJSD[order(-quary_FJSD$FJSD), ], 10)
  top_genes <- unique(rbind(top_log10_genes, top_FJSD_genes))
  top_genes <- subset(top_genes, gene %in% significant_genes$gene)
  Plot <- ggplot(quary_FJSD, aes(x = FJSD, y = log10_p_val_adj)) +
    geom_point(aes(color = ifelse(log10_p_val_adj >= -log10(0.05) & FJSD > FJSD_Cutoff, "Significant", "Non-significant")), size = 3) +
    theme_classic() +
    scale_color_manual(name = "Group", values = c("Significant" = "red", "Non-significant" = "gray")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    geom_vline(xintercept = FJSD_Cutoff, linetype = "dashed", color = "#F3B169") +
    geom_text_repel(data = top_genes, aes(label = gene),
                    color = "black", size = 5, box.padding = 0.5, point.padding = 0.5, segment.color = 'black') +
    labs(title = title, x = "FJSD Score", y = "-log10(adj P value)") +
    theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
          plot.title = element_text(hjust = 0.5, size = 30),
          axis.title.x = element_text(size = 25),
          axis.title.y = element_text(size = 25),
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          legend.position = "right",
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 18)) +
    annotate("text", x = max(quary_FJSD$FJSD) * 0.95, y = -log10(0.05) - 1, label = "-log10(0.05)", color = "blue", size = 5) +
    annotate("text", x = FJSD_Cutoff - 0.1, y = max(quary_FJSD$log10_p_val_adj) - 1, label = FJSD_Cutoff, color = "#F3B169", size = 5)
  return(Plot)
}

#' Visualize Gene Expression Patterns Across Samples
#'
#' @param Exp_raw Raw expression matrix (genes x samples)
#' @param gene_name Character vector of genes to visualize
#' @param cell_type Dataframe with two columns: "sample" (sample IDs) and "celltype" (annotations)
#' @param title Plot title. Default: "Gene Expression".
#'
#' @return ggplot object showing expression trends across samples
#' @export

FJSD_geom_line <- function(Exp_raw, gene_name, cell_type = FALSE, title = "Gene Expression") {
  library(ggplot2)
  library(tidyr)
  if (is.null(cell_type) || identical(cell_type, FALSE)) {
    gene_df <- data.frame(sample = colnames(Exp_raw), t(Exp_raw[gene_name,]))
    gene_long <- tidyr::gather(gene_df, gene, expression, -sample)
    gene_long$sample <- factor(gene_long$sample, levels = unique(gene_long$sample))
    gene_names <- unique(gene_long$gene)
    Plot <- ggplot(gene_long, aes(x = sample, y = expression, group = gene, color = gene)) +
      theme_classic() +
      geom_line() +
      geom_point() +
      labs(title = title, x = "Samples", y = "Expression") +
      scale_color_discrete(name = "Genes", limits = gene_name) +
      theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
            plot.title = element_text(hjust = 0.5, size = 20),
            axis.title.x = element_text(size = 15),
            axis.title.y = element_text(size = 15),
            axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
            axis.text.y = element_text(size = 10),
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 10))
  } else {
    if (!is.data.frame(cell_type)) {
      stop("cell_type must be a data frame or FALSE")
    }
    cell_type <- cell_type %>%
      group_by(celltype) %>%
      mutate(celltype1 = paste0(celltype, "_", row_number())) %>%
      ungroup()
    colname <- data.frame(sample = colnames(Exp_raw))
    cell_type1 <- left_join(colname, cell_type, by = "sample")
    gene_df <- data.frame(sample = cell_type1$celltype1, t(Exp_raw[gene_name,]))
    gene_long <- tidyr::gather(gene_df, gene, expression, -sample)
    gene_long$sample <- factor(gene_long$sample, levels = unique(gene_long$sample))
    gene_names <- unique(gene_long$gene)
    Plot <- ggplot(gene_long, aes(x = sample, y = expression, group = gene, color = gene)) +
      theme_classic() +
      geom_line() +
      geom_point() +
      labs(title = title, x = "Samples", y = "Expression") +
      scale_color_discrete(name = "Genes", limits = gene_name) +
      theme(panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
            plot.title = element_text(hjust = 0.5, size = 20),
            axis.title.x = element_text(size = 15),
            axis.title.y = element_text(size = 15),
            axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
            axis.text.y = element_text(size = 10),
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 10))
  }
  return(Plot)
}







