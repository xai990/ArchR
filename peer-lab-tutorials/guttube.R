library(ArchR)
set.seed(1)

# Configure
addArchRThreads(threads = 31) 
addArchRGenome('mm10')


# ################################################################################################
# Arrow files and project 

# Input files
data.dir <- '~/projects/kat/endoderm-atac/data/cellranger/'
inputFiles <- c(sprintf("%s/Lib1_Ant-1/fragments.tsv.gz", data.dir),
                sprintf("%s/Lib2_Ant-2/fragments.tsv.gz", data.dir),
                sprintf("%s/Lib3_Post-1/fragments.tsv.gz", data.dir),
                sprintf("%s/Lib4_Post-2/fragments.tsv.gz", data.dir),

                sprintf("%s/Lib5_guttube_GFPpos_1/fragments.tsv.gz", data.dir),
                sprintf("%s/Lib6_guttube_GFPneg_1/fragments.tsv.gz", data.dir),
                sprintf("%s/LIb7_guttube_GFPpos_2/fragments.tsv.gz", data.dir),
                sprintf("%s/Lib8_guttube_GFPneg_2/fragments.tsv.gz", data.dir)
               )
names(inputFiles) <- c('Anterior_Rep1',
                       'Anterior_Rep2', 
                       'Posterior_Rep1', 
                       'Posterior_Rep2', 

                       'GFPPos_Rep1', 
                       'GFPNeg_Rep1', 
                       'GFPPos_Rep2', 
                       'GFPNeg_Rep2'
                       )
# Create Arrow files
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 5000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  maxFrags = 1e+6,
  excludeChr = c('chrM'),
  removeFilteredCells = TRUE
)

# Doublet scores
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "LSI",
  LSIMethod = 1
)

# Create project
proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "Guttube_merged",
  copyArrows = FALSE
)
proj_dir <- "Guttube_merged/export/"

# Visualizations to check for doublets 
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSIPreDoublet")
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSIPreDoublet", name='UMAPPreDoublet')

# Plot doublets
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "DoubletScore", embedding = "UMAPPreDoublet")
plotPDF(p1, name = "Plot-Doublets.pdf",
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)




# ################################################################################################
# PReprocesing

# Fitlering doublets
proj <- filterDoublets(ArchRProj = proj)

# SVD, Clustering, UMAP
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI", scaleDims=FALSE)
proj <- addClusters(input = proj, reducedDims = "IterativeLSI")
proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI")
# Phenography
proj <- addCellColData(ArchRProj = proj, 
    data = sprintf("C%d", read.csv('Guttube_merged/export/all_cells_phenograph.csv')[, 2]), name = 'Phenograph', 
    cells = getCellNames(proj))

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

# Compartment annotations 
cells <- getCellNames(proj)
compartments <- rep("Other", length(cells))
names(compartments) <- cells
clusters <- getCellColData(proj)[, 'Phenograph']
compartments[clusters %in% c('C10')] = 'YsE'
compartments[clusters %in% c('C0', 'C3', 'C4', 'C12', 'C13', 'C17', 'C19')] = 'Mesoderm'
compartments[clusters %in% c('C1', 'C2', 'C5', 'C6', 'C7', 'C8', 'C9', 'C16')] = 'Gut'
proj <- addCellColData(proj, compartments, 'Compartment', names(compartments))

p <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Compartment", embedding = "UMAP")
plotPDF(p, name = "Compartments.pdf",
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

# Export 
# Gene scores
gene.scores <- getMatrixFromProject(proj)
scores <- assays(gene.scores)[['GeneScoreMatrix']]
rownames(scores) <- rowData(gene.scores)$name
write.csv(t(scores), 'Guttube_merged/export/all_cells_gene_scores.csv', quote=FALSE)


# Save 
proj <- saveArchRProject(ArchRProj = proj)






























































# ################################################################################################
# Gut tube subset 
gut_cells = names(compartments)[compartments == 'Gut']
gt_proj <- subsetArchRProject(proj, gut_cells, 'Guttube_endoderm')


# ################################################################################################
# PReprocesing

# SVD, Clustering, UMAP
res <- addIterativeLSI(ArchRProj = gt_proj, useMatrix = "TileMatrix", name = "IterativeLSI", scaleDims=FALSE, force=TRUE)
gt_proj <- res[[1]]
var_features <- res[[2]]

# GEne scores with selected features
# Artificial black list to exclude all non variable features
chrs <- getChromSizes(gt_proj)
blacklist <- setdiff(chrs, var_features)
gt_proj <- addGeneScoreMatrix(gt_proj, matrixName='GeneScoreMatrix', force=TRUE, blacklist=blacklist)


# Clustering and UMAP
gt_proj <- addClusters(input = gt_proj, reducedDims = "IterativeLSI", force=TRUE)
gt_proj <- addUMAP(ArchRProj = gt_proj, reducedDims = "IterativeLSI", force=TRUE)
gt_proj <- saveArchRProject(ArchRProj = gt_proj)



# ################################################################################################
# Plots



p1 <- plotEmbedding(ArchRProj = gt_proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = gt_proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = gt_proj, colorBy = "cellColData", name = "nFrags", embedding = "UMAP")
plotPDF(p1,p2, p3, name = "Plot-UMAP-Sample-Clusters.pdf",
        ArchRProj = gt_proj, addDOC = FALSE, width = 5, height = 5)


# Plot marker genes
gt_proj <- addImputeWeights(gt_proj)
markerGenes  <- c(
    "Nkx2-1",  
    "Nkx2-5",
    "Foxg1", 
    "Hoxb1",
    "Hoxc9" 
  )
p <- plotEmbedding(
    ArchRProj = gt_proj, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(gt_proj)
)
plotPDF(plotList = p, 
    name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", 
    ArchRProj = gt_proj, 
    addDOC = FALSE, width = 5, height = 5)


p <- plotEmbedding(
    ArchRProj = gt_proj, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
)
plotPDF(plotList = p, 
    name = "Plot-UMAP-Marker-Genes-W-UnImputed.pdf", 
    ArchRProj = gt_proj, 
    addDOC = FALSE, width = 5, height = 5)




# ################################################################################################
# Export

# Export
write.csv(getReducedDims(gt_proj), 'Guttube_endoderm/export/svd.csv', quote=FALSE)
write.csv(getEmbedding(gt_proj), 'Guttube_endoderm/export/umap.csv', quote=FALSE)
write.csv(getCellColData(gt_proj), 'Guttube_endoderm/export/cell_metadata.csv', quote=FALSE)


# Gene scores
gene.scores <- getMatrixFromProject(gt_proj)
scores <- assays(gene.scores)[['GeneScoreMatrix']]
scores <- as.matrix(scores)
rownames(scores) <- rowData(gene.scores)$name
write.csv(scores, 'Guttube_endoderm/export/gene_scores.csv')





# Peaks 
gt_proj <- addGroupCoverages(gt_proj, maxFragmentLength=147)
gt_proj <- addReproduciblePeakSet(gt_proj, maxPeaks=250000)
peaks <- getPeakSet(gt_proj)
# Counts
gt_proj <- addPeakMatrix(gt_proj, maxFragmentLength=147)
peak.counts <- getMatrixFromProject(gt_proj, 'PeakMatrix')

# Export counts
counts <- assays(peak.counts)[['PeakMatrix']]
writeMM(counts, 'Guttube_endoderm/export/peak_counts/counts.mtx')
write.csv(colnames(peak.counts), 'Guttube_endoderm/export/peak_counts/cells.csv', quote=FALSE)
names(peaks) <- sprintf("Peak%d", 1:length(peaks))
write.csv(as.data.frame(peaks), 'Guttube_endoderm/export/peak_counts/peaks', quote=FALSE)





# RNA integration
# Read in data 
rna_data_dir <- '~/projects/kat/endoderm-atac/data/rna/'
rna_counts <- readMM(sprintf("%s/cell_counts.mtx", rna_data_dir))
rna_cell_metadata <- read.csv(sprintf("%s/cell_metadata.csv", rna_data_dir))
rownames(rna_cell_metadata) <- rna_cell_metadata[, 1]
rna_gene_metadata <- read.csv(sprintf("%s/gene_metadata.csv", rna_data_dir))
rownames(rna_gene_metadata) <- rna_gene_metadata[, 1]
seRNA <- SummarizedExperiment(list(counts=t(rna_counts)), 
                              colData=rna_cell_metadata, rowData=rna_gene_metadata)

# outDir <- 'Guttube_endoderm/'

gt_proj <- addGeneIntegrationMatrix(
    ArchRProj = gt_proj, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    sampleCellsATAC = 30000,
    sampleCellsRNA = 30000,
    groupRNA = "CellType",
    nameCell = "RNACell",
    nameGroup = "RNACellType",
    nameScore = "RNACellTypeScore",
    force = TRUE
)

# Plot 
pal <- paletteDiscrete(values = colData(seRNA)$CellType)
p1 <- plotEmbedding(
    gt_proj, 
    colorBy = "cellColData", 
    name = "RNACellType", 
    pal = pal
)
plotPDF(p1, name = "Plot-UMAP-RNA-Integration.pdf", ArchRProj = gt_proj, addDOC = FALSE, width = 5, height = 5)

# Export
write.csv(getCellColData(gt_proj), 'Guttube_endoderm/export/cell_metadata.csv', quote=FALSE)


# Save 
gt_proj <- saveArchRProject(ArchRProj = gt_proj)




# Impute weights 
gt_proj <- addImputeWeights(gt_proj)
imputeWeights <- getImputeWeights(gt_proj)
markerGenes <- c('Nkx2-1', 'Nkx2-5', 'Pax8', 'Nkx2-3', 'Isl1', 'Otx2', 
  'Prrx2', 'Six1', 'Foxg1', 'Irx3', 'Hoxb1', 'Meis2', 'Gata6', 'Foxa3', 
  'Cdx2', 'Hoxa7', 'Hoxb8', 'Hoxc8', 'Hoxc9', 'Tlx2', 'Irx1', 'Ppy', 'Pdx1', 
  'Fabp1', 'Hhex', 'Socs2', 'Wfdc2', 'Krt19', 'Nepn')
markerMat <- ArchR:::.getMatrixValues(
  ArchRProj = gt_proj, 
  name = markerGenes, 
  matrixName = 'GeneScoreMatrix', 
  log2Norm = FALSE, 
)
imputedMat <- imputeMatrix(mat = as.matrix(markerMat), imputeWeights = imputeWeights)
write.csv(t(imputedMat), 'Guttube_endoderm/export/markers_imputed_exprs.csv', quote=FALSE)






























p <- plotBrowserTrack(
    ArchRProj = proj, 
    groupBy = "Clusters", 
    geneSymbol = c('Hoxc9'), 
    upstream = 50000,
    downstream = 50000
)
grid::grid.newpage()
grid::grid.draw(p$Hoxc9)











all_frags <- getFragmentsFromArrow('../data/ArchR/GFPPos_Rep1.arrow', 'chr18')
nfr_frags <- getFragmentsFromArrow('../data/ArchR/GFPPos_Rep1.arrow', 'chr18', maxFragmentLength=147)


