# script to download, manipulate and visualize gene expression data (GSE183947)
# setwd("C:/Users/sade0010/OneDrive - University of Oklahoma/University of Oklahoma/bioinformatics practice/bioinformagician")


# load libraries ----------------------------------------------------------


# for manipulation
library(dplyr)
library(GEOquery)

# for visualization 
library(ggplot2)

# for both
library(tidyverse)


# read in data ------------------------------------------------------------

# expression data
dat <- read.csv(file = "../bioinformagician/GSE183947_fpkm.csv")
dim(dat)

# get metadata from GEOquery package

gse <- getGEO(GEO = 'GSE183947', GSEMatrix = TRUE)

gse

metadata <- pData(phenoData(gse[[1]]))
head(metadata)
view(metadata)

# modify the metadata

metadata.modified <- metadata %>%
  select(1, 10, 11, 17) %>%
  rename(tissue = characteristics_ch1, metastasis = characteristics_ch1.1) %>%
  mutate(tissue = gsub("tissue:", "", tissue)) %>%
  mutate(metastasis = gsub("metastasis:", "", metastasis))


# reshaping the format of data to long

dat.long <- dat %>%
  rename(gene = X) %>%
  pivot_longer(!gene, names_to = "samples", values_to = "FPKM")
# gather(key = 'samples', value = 'FPKM', -gene) this is an alternative to pivot long

# join dataframe = dat.long + meta.modified

dat.long <- dat.long %>%
  left_join(., metadata.modified, by = c("samples" = "description"))

# Exploring the data ------------------------------------------------------


# explore data
dat.long %>%
  filter(gene =='BRCA1' | gene == 'BRCA2') %>%
  group_by(gene, tissue) %>%
  summarize(mean_FPKM = mean(FPKM),
            median_FPKM = median(FPKM)) %>%
  arrange(-mean_FPKM)




# plots -------------------------------------------------------------------

# 1. barplot

dat.long %>%
  filter(gene == 'BRCA1') %>%
  ggplot(., aes(x = samples, y = FPKM, fill = tissue)) +
  geom_col()

# 2. density

dat.long %>%
  filter(gene == 'BRCA1') %>%
  ggplot(., aes(x = FPKM, fill = tissue)) +
  geom_density(alpha = 0.5)

# 3. boxplot
dat.long %>%
  filter(gene == 'BRCA1') %>%
  ggplot(., aes(x = metastasis, y = FPKM)) +
  #geom_boxplot()
  geom_violin() + 
  geom_point()

# 4. scatterplot
dat.long %>%
  filter(gene == 'BRCA1' | gene == 'BRCA2') %>%
  pivot_wider(names_from = gene, values_from = FPKM) %>%
  ggplot(., aes(x = BRCA1, y = BRCA2, color = tissue)) +
  geom_point() + 
  geom_smooth(method = 'lm', se = FALSE)

# 5. heatmap
genes.of.interest <- c('BRCA1', 'BRCA2', 'TP53', 'ALK', 'MYCN')

pdf("heatmap_save2.pdf")
dat.long %>%
  filter(gene %in% genes.of.interest) %>%
  ggplot(., aes(x = samples, y = gene, fill = FPKM)) +
  geom_tile() +
  scale_fill_gradient(low = 'white', high = 'red')

dev.off()

ggsave(p, filename = 'heatmap_save1.pdf', width = 10, height = 8)
