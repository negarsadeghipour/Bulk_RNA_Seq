# Script to perform WGCNA
# setwd("C:/Users/sade0010/OneDrive - University of Oklahoma/University of Oklahoma/bioinformatics practice/Bioinformagician-WGCNA")

library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)

allowWGCNAThreads()#allow multi-threading (optional)


# 1. Fetch Data -----------------------------------------------------------

wd = getwd()

data <- read.delim(paste0(wd,'/GSE152418_p20047_Study1_RawCounts1.txt'), header = T)


# get metadata
geo_id <- "GSE152418"
gse <- getGEO(geo_id, GSEMatrix = TRUE)
phenoData <- pData(phenoData(gse[[1]]))
head(phenoData)
phenoData <- phenoData[,c(1,2,46:50)]

# prepare data
data[1:10, 1:10]

data <- data %>%
  gather(key = "samples", value = "counts", -ENSEMBLID) %>%
  mutate(samples = gsub('\\.', '-', samples)) %>%
  inner_join(., phenoData, by = c('samples' = 'title')) %>%
  select(1, 3, 4) %>%
  spread(key = "geo_accession", value = "counts") %>%
  column_to_rownames(var = "ENSEMBLID")



# 2. QC - outlier detection -----------------------------------------------

# detect outlier genes

gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# remove genes that are detected as outliers
data <- data[gsg$goodGenes == TRUE,]

# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(data)), method = "average")
plot(htree)

# pca - method 2
pca <- prcomp(t(data))
pca.data <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.data <- as.data.frame(pca.data)

ggplot(pca.data, aes(PC1, PC2)) + 
  geom_point() + 
  geom_text(label = rownames(pca.data)) + 
  labs(x = paste0('PC1:', pca.var.percent[1], '%'),
       y = paste0('PC2:', pca.var.percent[2], '%'))


### NOTE: make sure there are no batch effects in the data
samples.to.be.excluded <- c('GSM4615000', 'GSM4614993', 'GSM4624995')
data.subset <- data[,!(colnames(data) %in% samples.to.be.excluded)]


# 3. Normalization --------------------------------------------------------

# create a deseq2 dataset

colData <- phenoData %>%
  filter(!row.names(.) %in% samples.to.be.excluded)

# fixing column names in colData
names(colData)
names(colData) <- gsub(':ch1', '', names(colData))
names(colData) <- gsub('\\s', '_', names(colData))

# making the rownames and column names identical
all(rownames(colData) %in% colnames(data.subset))
all(rownames(colData) == colnames(data.subset))

# create dds

dds <- DESeqDataSetFromMatrix(countData = data.subset,
                              colData = colData,
                              design = ~1) # not specifing model

## remove all genes with counts < 15 in more than 75% of samples (31*0.75 = 23.75)
## suggested by WGCNA on RNAseq FAQ

dds75 <- dds[rowSums(counts(dds) >= 15) >= 24,]
nrow(dds75)

# perform variance stabilization
dds_norm <- vst(dds75)

# get normalized counts
norm.counts <- assay(dds_norm) %>%
  t()


# 4. Network Construction -------------------------------------------------

# choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)

sft.data <- sft$fitIndices

# visualization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) + 
  geom_hline(yintercept = 0.8, color = 'red') + 
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)

# convert matrix to numeric

norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 18
temp_cor <- cor
cor <- WGCNA::cor

# dendrogram --------------------------------------------------------------


# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 14000,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)

cor <- temp_cor

module_eigengenes <- bwnet$MEs

# print out a preview
head(module_eigengenes)


# get number of genes for each module
table(bwnet$colors)

# Plot the dendrogram and the module colors efore and after merging
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang = 0.03,
                    guideHang = 0.05)

# gray module = all genes that don't fall into other modules were assigned to the gray module

traits <- colData %>%
  mutate(disease_state_bin = ifelse(grepl('COVID', disease_state), 1, 0)) %>%
  select(8)

# binarize categorical variables

colData$severity <- factor(colData$severity, levels = c("Healthy", "Convalescent", "ICU", "Moderate", "Severe"))

severity.out <- binarizeCategoricalColumns(colData$severity,
                                           includePairwise = FALSE,
                                           includeLevelVsAll = TRUE,
                                           minCount = 1)

traits <- cbind(traits, severity.out)


# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

# visualize module-trait associations as a heatmap

heatmap.data <- merge(module_eigengenes, traits, by = "row.names")

head(heatmap.data)


heatmap.data <- heatmap.data %>%
  column_to_rownames(var = "Row.names")

CorLevelPlot(heatmap.data, 
             x = names(heatmap.data)[21:25],
             y = names(heatmap.data)[1:20],
             col = c("blue", "skyblue", "white", "pink", "red"))

module.gene.mapping <- as.data.frame(bwnet$colors)
module.gene.mapping %>%
  filter(`bwnet$colors` == 'turquoise') %>%
  rownames()


# 6B. Intramodular analysis: Identifying driger genes ---------------------

# calculate the module membership and the associated p-values

# The module membership/intramodular connectivity is calculated as the correlation of the eigengenes and the genes
# This quantifies the similarity of all genes on the array to every module.

module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)

module.membership.measure.pvals[1:10, 1:10]

# calculate the gene significance and associated p-values

gene.signf.corr <- cor(norm.counts, traits$data.Severe.vs.all, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)

gene.signf.corr.pvals %>%
  as.data.frame() %>%
  arrange(V1) %>%
  head(25)