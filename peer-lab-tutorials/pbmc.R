library(ArchR)
set.seed(1)

# Configure
addArchRThreads(threads = 31) 
addArchRGenome('hg19')


# ################################################################################################
# Arrow files and project 

# Input files
data.dir <- '~/projects/kat/endoderm-atac/data/10x_pbmcs/'
inputFiles <- c(sprintf("%s/atac_pbmc_10k_nextgem_fragments.tsv.gz", data.dir)
               )
names(inputFiles) <- c('10X_PBMCs10k'
                       )
# Create Arrow files
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 1, #Dont set this too high because you can always increase later
  filterFrags =3000, 
  addTileMat = TRUE,
  addGeneScoreMat = FALSE,
  excludeChr = c('chrM'),
  removeFilteredCells = TRUE
)


# Create project
proj_name <- "PBMCs_10k"
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = proj_name,
)




# # Doublet scores
# doubScores <- addDoubletScores(
#   input = proj,
#   k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
#   knnMethod = "LSI",
#   LSIMethod = 1
# )




# # Visualizations to check for doublets 
# proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSIPreDoublet")
# proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSIPreDoublet", name='UMAPPreDoublet')

# # Plot doublets
# p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "DoubletScore", embedding = "UMAPPreDoublet")
# plotPDF(p1, name = "Plot-Doublets.pdf",
#         ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)






# ################################################################################################
# PReprocesing

# SVD, Clustering, UMAP
res <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI", scaleDims=FALSE, force=TRUE)
proj <- res[[1]]
var_features <- res[[2]]

proj <- addClusters(input = proj, reducedDims = "IterativeLSI")
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI")

# GEne scores with selected features
# Artificial black list to exclude all non variable features
chrs <- getChromSizes(proj)
var_features_gr <- GRanges(var_features$seqnames, IRanges(var_features$start, var_features$start + 500))
blacklist <- setdiff(chrs, var_features_gr)
proj <- addGeneScoreMatrix(proj, matrixName='GeneScoreMatrix', force=TRUE, blacklist=blacklist)



# Peaks 
proj <- addGroupCoverages(proj, maxFragmentLength=147)
proj <- addReproduciblePeakSet(proj)
# Counts
proj <- addPeakMatrix(proj, maxFragmentLength=147, ceiling=10^9)

# Save 
proj <- saveArchRProject(ArchRProj = proj)




# ################################################################################################
# Prelim annotations

# Plots 
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "nFrags", embedding = "UMAP")
plotPDF(p1,p2, p3, name = "Plot-UMAP-Sample-Clusters.pdf",
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)




# ################################################################################################

# Export
write.csv(getReducedDims(proj), sprintf('%s/export/svd.csv', proj_name), quote=FALSE)
write.csv(getEmbedding(proj), sprintf('%s/export/umap.csv', proj_name), quote=FALSE)
write.csv(getCellColData(proj), sprintf('%s/export/cell_metadata.csv', proj_name), quote=FALSE)


# Gene scores
gene.scores <- getMatrixFromProject(proj)
scores <- assays(gene.scores)[['GeneScoreMatrix']]
scores <- as.matrix(scores)
rownames(scores) <- rowData(gene.scores)$name
write.csv(scores, sprintf('%s/export/all_cells_gene_scores.csv', proj_name), quote=FALSE)

# Peak counts
peaks <- getPeakSet(proj)
peak.counts <- getMatrixFromProject(proj, 'PeakMatrix')
# Reorder peaks 

# Chromosome order
chr_order <- sort(seqlevels(peaks))
reordered_features <- list()
for(chr in chr_order)
    reordered_features[[chr]] = peaks[seqnames(peaks) == chr]
reordered_features <- Reduce("c", reordered_features)    



# Export counts
counts <- assays(peak.counts)[['PeakMatrix']]
writeMM(counts, sprintf('%s/export/peak_counts/counts.mtx', proj_name))
write.csv(colnames(peak.counts), sprintf('%s/export/peak_counts/cells.csv', proj_name), quote=FALSE)
names(reordered_features) <- sprintf("Peak%d", 1:length(reordered_features))
write.csv(as.data.frame(reordered_features), sprintf('%s/export/peak_counts/peaks.csv', proj_name), quote=FALSE)



