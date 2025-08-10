#--------------------------------------------------------------
# Load required libraries
#--------------------------------------------------------------
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(rtracklayer)
library(Matrix)
BiocManager::install("Signac", force = TRUE)
set.seed(1234)

#--------------------------------------------------------------
# Define file paths and genome build
#--------------------------------------------------------------
h5_file <- "GSM5342787_scA_filtered_peak_bc_matrix.h5"
fragments_file <- "fragments.tsv.gz"
blacklist_file <- "hg38-blacklist.v2.bed"
genome_build <- "hg38"

#--------------------------------------------------------------
# Load scATAC-seq count matrix from 10x Genomics (.h5)
#--------------------------------------------------------------
atac_counts <- Read10X_h5(h5_file)

#--------------------------------------------------------------
# Create fragment object
#--------------------------------------------------------------
fragments <- CreateFragmentObject(
  path = fragments_file,
  cells = colnames(atac_counts)
)

#--------------------------------------------------------------
# Create ChromatinAssay and Seurat object
#--------------------------------------------------------------
atac_assay <- CreateChromatinAssay(
  counts = atac_counts,
  fragments = fragments,
  sep = c(":", "-"),
  genome = genome_build,
  min.cells = 10,
  min.features = 200
)

atac_obj <- CreateSeuratObject(
  counts = atac_assay,
  assay = "peaks",
  project = "GSE175621_ATAC"
)


library(AnnotationHub)
ah <- AnnotationHub()

# Search for the Ensembl 98 EnsDb for Homo sapiens on AnnotationHub
query(ah, "EnsDb.Hsapiens.v98")
ensdb_v98 <- ah[["AH75011"]]
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = ensdb_v98)

# change to UCSC style since the data was mapped to hg38
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
# add the gene information to the object
Annotation(atac_obj) <- annotations

atac_obj[['peaks']]
# compute nucleosome signal score per cell
atac_obj <- NucleosomeSignal(object = atac_obj)

# compute TSS enrichment score per cell
atac_obj <- TSSEnrichment(object = atac_obj)


#
# scATAC_qc_metrics.R

# Read tables
pf <- read.table("passed_filters.txt", col.names = c("barcode", "passed_filters"))
prf <- read.table("peak_region_fragments.txt", col.names = c("barcode", "peak_region_fragments"))
brf <- read.table("blacklist_region_fragments.txt", col.names = c("barcode", "blacklist_region_fragments"))

# Merge all
qc <- Reduce(function(x, y) merge(x, y, by = "barcode", all = TRUE), list(pf, prf, brf))
qc[is.na(qc)] <- 0  # Replace missing with 0

# Compute percentages
qc$pct_reads_in_peaks <- qc$peak_region_fragments / qc$passed_filters * 100
qc$blacklist_ratio <- qc$blacklist_region_fragments / qc$peak_region_fragments

# Save
write.table(qc, "scATAC_qc_metrics.txt", sep = "\t", quote = FALSE, row.names = FALSE)
rownames(qc) <- qc$barcode
qc$barcode <- NULL
# Add columns to meta.data
atac_obj <- AddMetaData(atac_obj, metadata = qc)

#
options(repr.plot.width=14, repr.plot.height=4)
VlnPlot(
  object = atac_obj,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)


DensityScatter(atac_obj, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
options(repr.plot.width=14, repr.plot.height=4)

atac_obj <- subset(
  x = atac_obj,
  subset = peak_region_fragments > 2000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)
atac_obj


devtools::install_github("cellgeni/sceasy")
# Tell reticulate where to find the conda executable
reticulate::conda_binary("/Users/abdelazizawad/miniforge3/condabin/conda")
options(reticulate.conda.binary = "/Users/abdelazizawad/miniforge3/bin/conda")
library(reticulate)
# Now, use_condaenv() should be able to find the environment
use_condaenv("sc_env", conda = "/Users/abdelazizawad/miniforge3/bin/conda", required = TRUE)
library(sceasy)
#use_python("/Users/abdelazizawad/miniforge3/bin/python", required = TRUE)
sc <- import("scanpy", convert = FALSE)
scvi <- import("scvi", convert = FALSE)
adata <- convertFormat(pbmc, from="seurat", to="anndata", main_layer="counts", assay="peaks", drop_single_values=FALSE)
print(adata) # Note generally in Python, dataset conventions are obs x var
