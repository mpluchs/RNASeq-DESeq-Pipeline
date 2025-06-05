# DESeq2 Single Comparison Analysis Script
# This script is designed to run a single comparison between two groups

# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if we have enough arguments
if (length(args) < 4) {
  stop("Usage: Rscript deseq2_single_comparison.R <group1_dir> <group2_dir> <output_dir> <comparison_name>")
}

# Parse arguments
group1_dir <- args[1]
group2_dir <- args[2]
output_dir <- args[3]
comparison_name <- args[4]

# Create output directory if it doesn't exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Print parameters for confirmation
cat("Processing comparison:", comparison_name, "\n")
cat("Group 1 directory:", group1_dir, "\n")
cat("Group 2 directory:", group2_dir, "\n")
cat("Output directory:", output_dir, "\n")

# Function to list all files in a directory (tab-delimited count files)
get_files <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    stop(paste("Directory does not exist:", dir_path))
  }
  
  # First, list all files in the directory to see what's there
  all_files <- list.files(dir_path, full.names = FALSE, recursive = TRUE)
  cat("Files in directory:", paste(all_files, collapse=", "), "\n")
  
  # Look specifically for CSV files
  files <- list.files(dir_path, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)
  
  if (length(files) == 0) {
    stop(paste("No CSV files found in", dir_path, 
               "\nPlease check the directory contents and file extensions."))
  }
  
  cat("Found these CSV files that will be used:\n")
  for (f in files) {
    cat(" -", basename(f), "\n")
  }
  
  return(files)
}

# Function to read count files in CSV format
read_count_files <- function(files) {
  # Initialize variables
  count_data_list <- list()
  
  # Read each file
  for (file_path in files) {
    sample_name <- basename(file_path)
    cat("Reading file:", sample_name, "\n")
    
    # Try to read the file
    tryCatch({
      # First try comma-separated
      data <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)
      
      # Print columns to help debugging
      cat("  Columns in file:", paste(colnames(data), collapse=", "), "\n")
      
      # Extract Geneid and Count columns
      if ("Geneid" %in% colnames(data) && "Count" %in% colnames(data)) {
        gene_ids <- data$Geneid
        counts <- data$Count
        
        # Store data
        count_data <- data.frame(
          gene_id = gene_ids,
          count = counts
        )
        count_data_list[[sample_name]] <- count_data
      } else if ("gene_id" %in% colnames(data) && "Count" %in% colnames(data)) {
        gene_ids <- data$gene_id
        counts <- data$Count
        
        # Store data
        count_data <- data.frame(
          gene_id = gene_ids,
          count = counts
        )
        count_data_list[[sample_name]] <- count_data
      } else {
        # Try to guess columns - first column as gene ID, last as count
        if (ncol(data) >= 2) {
          gene_ids <- data[[1]]
          counts <- data[[ncol(data)]]
          
          # Store data
          count_data <- data.frame(
            gene_id = gene_ids,
            count = counts
          )
          count_data_list[[sample_name]] <- count_data
        } else {
          cat("  Error: Cannot identify gene ID and count columns in file\n")
        }
      }
    }, error = function(e) {
      cat("  Error reading file:", e$message, "\n")
    })
  }
  
  # Check if we have any data
  if (length(count_data_list) == 0) {
    stop("Failed to read any count data")
  }
  
  # Find common genes across all samples
  cat("Finding common genes across all samples...\n")
  gene_lists <- lapply(count_data_list, function(x) x$gene_id)
  common_genes <- Reduce(intersect, gene_lists)
  
  cat("Found", length(common_genes), "common genes across all samples\n")
  
  # Create count matrix with only common genes
  cat("Creating count matrix...\n")
  counts_matrix <- matrix(0, nrow = length(common_genes), ncol = length(count_data_list))
  rownames(counts_matrix) <- common_genes
  colnames(counts_matrix) <- names(count_data_list)
  
  for (i in seq_along(count_data_list)) {
    sample_name <- names(count_data_list)[i]
    sample_data <- count_data_list[[i]]
    
    # Find indices of common genes in this sample
    indices <- match(common_genes, sample_data$gene_id)
    
    # Add counts to matrix
    counts_matrix[, i] <- sample_data$count[indices]
  }
  
  return(counts_matrix)
}

# Define rowVars function for PCA plotting
rowVars <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...) / (ncol(x) - 1)
}

# Get files for both groups
group1_files <- get_files(group1_dir)
group2_files <- get_files(group2_dir)
all_files <- c(group1_files, group2_files)

# Read count data
cat("Reading count data from all files...\n")
counts_matrix <- read_count_files(all_files)
cat("Count matrix dimensions:", nrow(counts_matrix), "x", ncol(counts_matrix), "\n")

# Create condition factor
cat("Setting up experimental design...\n")
condition <- factor(c(rep("group1", length(group1_files)), rep("group2", length(group2_files))))
sample_names <- colnames(counts_matrix)

# Verify that the number of samples matches the number of conditions
if (length(condition) != ncol(counts_matrix)) {
  stop("Mismatch between number of samples (", ncol(counts_matrix), 
       ") and number of conditions (", length(condition), ")")
}

# Create column data
coldata <- data.frame(
  row.names = sample_names,
  condition = condition
)

# Double-check that row names of coldata match column names of counts_matrix
cat("Verifying sample names match between count matrix and condition data...\n")
cat("Count matrix column names:", paste(colnames(counts_matrix), collapse=", "), "\n")
cat("Column data row names:", paste(rownames(coldata), collapse=", "), "\n")

# Make sure they match exactly
all_match <- all(rownames(coldata) == colnames(counts_matrix))
cat("Do all names match exactly?", all_match, "\n")

if (!all_match) {
  cat("Warning: Sample names don't match exactly. Fixing...\n")
  # Fix the coldata to match counts_matrix
  coldata <- data.frame(
    row.names = colnames(counts_matrix),
    condition = condition
  )
  cat("Fixed column data row names:", paste(rownames(coldata), collapse=", "), "\n")
}

# Create DESeq dataset
cat("Creating DESeq dataset...\n")
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = coldata,
  design = ~condition
)

# Set reference level
cat("Setting reference level...\n")
dds$condition <- relevel(dds$condition, ref = "group1")

# Filter low count genes
cat("Filtering low count genes...\n")
dds <- dds[rowSums(counts(dds)) >= 10, ]
cat("Genes remaining after filtering:", nrow(dds), "\n")

# Run DESeq
cat("Running DESeq...\n")
dds <- DESeq(dds)

# Get results
cat("Extracting results...\n")
res <- results(dds, alpha = 0.05)

# Summarize results
cat("\nResults summary:\n")
print(summary(res))

# Save results to file
cat("Saving results to CSV...\n")
output_file <- file.path(output_dir, paste0(comparison_name, "_DESeq2.csv"))
resOrdered <- res[order(res$padj), ]
resdata <- data.frame(gene=rownames(resOrdered), resOrdered)

# Add normalized counts
norm_counts <- counts(dds, normalized=TRUE)
for (col in colnames(norm_counts)) {
  resdata[[col]] <- norm_counts[resdata$gene, col]
}

write.csv(resdata, file=output_file, row.names=FALSE)
cat("Results saved to:", output_file, "\n")

# Create p-value histogram
cat("Creating p-value histogram...\n")
pdf_file <- file.path(output_dir, paste0(comparison_name, "_pvalue_hist.pdf"))
pdf(pdf_file)
hist(res$pvalue, breaks=20, main=paste("P-value histogram -", comparison_name), 
     xlab="p-value", col="lightblue")
dev.off()

# MA plot
cat("Creating MA plot...\n")
ma_file <- file.path(output_dir, paste0(comparison_name, "_MA_plot.pdf"))
pdf(ma_file)
DESeq2::plotMA(res, alpha = 0.05, main = comparison_name)
dev.off()

# PCA plot
cat("Creating PCA plot...\n")
pca_file <- file.path(output_dir, paste0(comparison_name, "_PCA_plot.pdf"))
pdf(pca_file)

# Variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# PCA plotting function
plotPCA <- function(ntop=500) {
  # Get most variable genes
  Pvars <- rowVars(assay(vsd))
  select <- order(Pvars, decreasing = TRUE)[seq_len(min(ntop, length(Pvars)))]
  
  # Calculate PCA
  PCA <- prcomp(t(assay(vsd)[select, ]), scale = TRUE, center = TRUE)
  percentVar <- round(100 * PCA$sdev^2 / sum(PCA$sdev^2), 1)
  
  # Create data frame for plotting
  dataGG <- data.frame(
    PC1 = PCA$x[, 1],
    PC2 = PCA$x[, 2],
    condition = colData(vsd)$condition,
    sample = rownames(colData(vsd))
  )
  
  # Create plot
  p <- ggplot(dataGG, aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 4) +
    ggtitle(paste(comparison_name, "- Top", ntop, "Variable Genes")) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    theme_bw() +
    geom_label_repel(aes(label = sample), box.padding = 0.5, max.overlaps = 20)
  
  print(p)
}

# Create PCA plots with different numbers of top genes
plotPCA(500)
plotPCA(1000)
dev.off()

# Create summary file with key statistics
cat("Creating summary file...\n")
summary_file <- file.path(output_dir, paste0(comparison_name, "_summary.txt"))
sink(summary_file)
cat("Comparison:", comparison_name, "\n")
cat("Group 1:", basename(group1_dir), "\n")
cat("Group 2:", basename(group2_dir), "\n")
cat("Total genes analyzed:", nrow(dds), "\n")
cat("Significant genes (FDR < 0.05):", sum(res$padj < 0.05, na.rm = TRUE), "\n")
cat("Up-regulated genes:", sum(res$log2FoldChange > 0 & res$padj < 0.05, na.rm = TRUE), "\n")
cat("Down-regulated genes:", sum(res$log2FoldChange < 0 & res$padj < 0.05, na.rm = TRUE), "\n")
cat("\nTop 20 differentially expressed genes:\n")
top_genes <- head(resOrdered[!is.na(resOrdered$padj), ], 20)
print(data.frame(
  Gene = rownames(top_genes),
  log2FoldChange = round(top_genes$log2FoldChange, 2),
  padj = format(top_genes$padj, scientific = TRUE, digits = 3)
))
sink()

cat("\nAnalysis complete! All output saved to:", output_dir, "\n")
cat("Summary statistics saved to:", summary_file, "\n")

# Write a simple status file to indicate completion
status_file <- file.path(output_dir, paste0(comparison_name, "_status.txt"))
write(paste("Completed successfully at", Sys.time()), file = status_file)