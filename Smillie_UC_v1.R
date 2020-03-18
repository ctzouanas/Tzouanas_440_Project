require(Seurat, quietly = T)
library(data.table)
library(ggplot2)
library(gridExtra)
library(gtable)
library(plotly)
library(scales)
library(spam)
library(dplyr)
library(patchwork)

project_name = "20_440_Smillie_UC"
figdir <-  "Plot_Folder"
objdir <- "Object_Folder"
print('hello darkness my old friend')

spam.pathname <- "RawData_Folder/gene_sorted-Imm.matrix.mtx"
gene.names.sheet <- "RawData_Folder/Imm.genes.tsv"
barcodes.sheet <- "RawData_Folder/Imm.barcodes2.tsv"

colon.spam <- read.MM(spam.pathname)
colon.gene.names <- read.table(gene.names.sheet, header = F)
colon.barcodes <- read.table(barcodes.sheet, header = F)

print(head(colon.spam))
print(head(colon.gene.names))
print(head(colon.barcodes))
print('I\'ve come to talk with you again')

n.genes <- colon.spam@dimension[1]
n.cells <- colon.spam@dimension[2]
print(n.cells)
downsample.factor <- 0.05
n.cells.to.sample <- floor(n.cells*downsample.factor)
cell.indices <- sample.int(n.cells, n.cells.to.sample)
downsampled.spam <- colon.spam[, cell.indices]
downsampled.barcodes <- colon.barcodes[cell.indices,1]

test.matrix <- as.matrix(downsampled.spam)
test.df <- as.data.frame(test.matrix, row.names = colon.gene.names, col.names = downsampled.barcodes)

rownames(test.df) <- colon.gene.names[[colnames(colon.gene.names)[1]]]
colnames(test.df) <- as.character(downsampled.barcodes)

smillie.downsampled.imm <- CreateSeuratObject(counts = test.df, project = "CT_JH_20_440_UC_Smillie")
smillie.downsampled.imm[["percent.mt"]] <- PercentageFeatureSet(smillie.downsampled.imm, pattern = "^MT-")

smillie.downsampled.imm <- SCTransform(smillie.downsampled.imm, verbose = FALSE)

smillie.downsampled.imm <- RunPCA(smillie.downsampled.imm, verbose = FALSE)
smillie.downsampled.imm <- RunUMAP(smillie.downsampled.imm, dims = 1:10, verbose = FALSE)

smillie.downsampled.imm <- FindNeighbors(smillie.downsampled.imm, dims = 1:10, verbose = FALSE)
smillie.downsampled.imm <- FindClusters(smillie.downsampled.imm, verbose = FALSE)

metadata.df <- NULL
downsampled.IDs.char <- as.character(downsampled.barcodes)
split.IDs.list <- lapply(downsampled.IDs.char, strsplit, "\\.")
metadata.df <- as.data.frame(split.IDs.list)
colnames(metadata.df) <- downsampled.IDs.char
rownames(metadata.df) <- c("Patient", "Location", "Cell Barcode")
metadata.df <- as.data.frame(t(metadata.df))

smillie.downsampled.imm$patient.id <- as.character(metadata.df$Patient)

colitis.patients <- c("N19", "N26", "N23", "N14", "N12", "N9", "N7", "N106", "N44", "N49", "N539", "N52", "N111", "N58", "N661", "N24", "N50", "N897")
healthy.patients <- c("N8", "N10", "N11", "N13", "N15", "N16", "N17", "N18", "N20", "N21", "N51", "N46")

logic.health <- smillie.downsampled.imm$patient.id %in% healthy.patients
health.status.vec <- vector(mode = "character", length = length(logic.health))
health.status.vec[logic.health] = "healthy"
health.status.vec[!logic.health] = "colitis"
smillie.downsampled.imm$health.status <-  health.status.vec

saveRDS(smillie.downsampled.imm, paste0(objdir, '/', project_name, '_initial_clustered_withmetadata.rds'), compress= T)

jpeg(paste(figdir,"/",project_name,"_umap_split_by_health.jpg", sep = ""))
DimPlot(smillie.downsampled.imm, split.by = "health.status") + NoLegend()
dev.off()