# Tint Added for annotation

#!/usr/bin/env Rscript

# Load libraries
suppressMessages(library(Seurat))
suppressMessages(library(SingleR))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(celldex))

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: SeuratAnnotate.R <Path_to_STARsoloOutputs> <Path_to_Output_CSV>")
}

sample_base_dir <- args[1]      # STARsolo output folder (one subfolder per sample)
output_csv <- args[2]           # Output file (final_cell_info.csv)

# Load default reference
ref <- celldex::HumanPrimaryCellAtlasData()

# Initialize list
seurat_objects <- list()

# Get sample directories
sample_dirs <- list.dirs(sample_base_dir, full.names = TRUE, recursive = FALSE)

# Read each sample
for (sample_dir in sample_dirs) {
  # Read the filtered STARsolo output: matrix.mtx + features.tsv + barcodes.tsv
  counts <- Read10X(data.dir = file.path(sample_dir, "Gene", "filtered"))
  
  # Create Seurat object
  seu_sample <- CreateSeuratObject(counts = counts)
  
  # OPTIONAL: Calculate mitochondrial percentage if MT genes available
  if (any(grepl("^MT-", rownames(seu_sample)))) {
    seu_sample[["percent.mt"]] <- PercentageFeatureSet(seu_sample, pattern = "^MT-")
    seu_sample <- subset(seu_sample, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  } else {
    seu_sample <- subset(seu_sample, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
  }
  
  # SCTransform normalization
  seu_sample <- SCTransform(seu_sample, verbose = FALSE, conserve.memory = TRUE)
  
  # Save object
  seurat_objects[[basename(sample_dir)]] <- seu_sample
}

# Integration
features <- SelectIntegrationFeatures(object.list = seurat_objects)
seurat_integrated <- IntegrateData(object.list = seurat_objects, features.to.integrate = features)

# PCA
seurat_integrated <- RunPCA(seurat_integrated, npcs = 50, verbose = FALSE)

# Find significant PCs
stdv <- seurat_integrated[["pca"]]@stdev
percent_stdv <- (stdv / sum(stdv)) * 100
cumulative <- cumsum(percent_stdv)
co1 <- which(cumulative > 90 & percent_stdv < 5)[1]
co2 <- sort(which((percent_stdv[1:(length(percent_stdv) - 1)] - percent_stdv[2:length(percent_stdv)]) > 0.1), decreasing = TRUE)[1] + 1
min_pc <- min(co1, co2)

# Clustering
seurat_integrated <- FindNeighbors(seurat_integrated, dims = 1:min_pc)
seurat_integrated <- FindClusters(seurat_integrated, resolution = 1)

# Convert to SingleCellExperiment
sce_integrated <- as.SingleCellExperiment(seurat_integrated)

# SingleR annotation
singleR_results <- SingleR(test = sce_integrated, ref = ref, labels = ref$label.fine)

# Add SingleR annotations to Seurat
seurat_integrated$SingleR_Annotation <- singleR_results$pruned.labels

# Save final barcode and cell type information
cell_info <- data.frame(
  Barcode = colnames(seurat_integrated),
  CellType = seurat_integrated$SingleR_Annotation
)

# Write to CSV
write.csv(cell_info, file = output_csv, row.names = FALSE)