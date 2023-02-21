# Info -------------------------------------------------------------------------

# Temperature effect on gene expression
# Subset of 54 individuals
# Temperature treatment in the lab was used as categorical factor
# 
# CL - 2023
# 

# Libraries---------------------------------------------------------------------

library(tidyverse)

library(vegan)
library(DESeq2)

library(eulerr)
library(pheatmap)
library(ggplotify)
library(ggpubr)

library(topGO)
library(trinotateR)

library(WGCNA)

library(readxl)
library(here)
library(vcfR)
library(adegenet)

library(dendextend)

library(pcadapt)
library(qvalue)

library(reshape2)

library(RColorBrewer)

library(dartR)

library(lsmeans)
library(effects)

library(ggdendro)
library(ggVennDiagram)

library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)

# Map - sampling sites ---------------------------------------------------------

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

lati <- c(48.5921923, 45.3772867, 50.3060011)
long <- c(-68.5915477, -61.0402811, -54.27895)
sites <- c("St-Lawrence\nEstuary\n(SLE)", 
           "Eastern\nScotian Shelf\n(ESS)", 
           "Northeast\nNewfoundland\nCoast (NNC)")

map.fig <- ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(-72, -51), ylim = c(43, 52), expand = FALSE) +
  theme_bw() +
  annotation_scale(location = "br", width_hint = 0.25, bar_cols = c("grey", "white"),
                   line_col = "grey30") +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(0, "in"), pad_y = unit(0.1, "in"),
                         style = north_arrow_fancy_orienteering(
                           line_col = "grey30", 
                           fill = c("white", "grey"),
                         )) +
  annotate(geom = "point", x = long, y = lati, shape = c(17, 19, 8), 
           color = "black", size = 5, stroke = 1.1) +
  annotate(geom = "label", x = long, y = lati, label = sites, 
           color = "black", size = 3.5, fontface = "bold",
           vjust = 1.3, hjust = 0.5,
           fill = alpha(c("white"),0.75), label.size = NA) +
  xlab("") + ylab("")
map.fig

# Data and Pre-filtering--------------------------------------------------------

## Data 
metaData <- read.table("00_Data/metaData_shrimp.txt", header = T)
metaData$Origin <- factor(metaData$Origin, levels = c("NNC", "SLE", "ESS"))
metaData$Treatment <- factor(metaData$Treatment)

counts <- read.table("00_Data/aa.isoforms.counts.54.samples",
                     header = TRUE,
                     row.names = 1)

metaData <- metaData[with(metaData, order(Origin, Temperature)), ]
counts <- counts[, match(metaData$ID, colnames(counts))]

## Shrimp size

size <- metaData$Size_mm
mean(size, na.rm = T)
sd(size, na.rm = T)

mod.size <- aov(size ~ Treatment + Origin, metaData)
plot(mod.size)
summary(mod.size)

## Pre-filtering

### Library size & cutoff:
min.lib.size <- min(colSums(counts)) / 10^6
cutoff <- 10 / min.lib.size

### Minimum count in at least 6 samples (because each group contains 6 replicates)
keep <- rowSums(edgeR::cpm(counts) > cutoff) >= 6
table(keep)

counts.filtered <- counts[keep, ]

# Output data
write.csv(counts.filtered, file = "00_Data/01_counts.filtered")

# individual order
identical(colnames(counts.filtered), metaData$ID)

# Design
design <- as.formula(~ Origin + Treatment + Origin:Treatment)
ddsObj <- DESeqDataSetFromMatrix(
  countData = counts.filtered,
  colData = metaData,
  design = design
)

# Normalized data
count.vst <- vst(ddsObj, blind = T)
write.csv(assay(count.vst), "02_Results/count_vst_54samples.csv")
write.csv(metaData, "02_Results/metaData_vst_54samples.csv")

# RDA RNA-seq ------------------------------------------------------------------

# Euclidean distance
dist.count <- vegdist(t(assay(count.vst)), method = "euclidean")

# Tank replicate effect:
rda.tank <- capscale(dist.count ~ Tank + Condition(Origin * Treatment), metaData)
RsquareAdj(rda.tank)
anova.cca(rda.tank)

# Temperature as categorical factor:
rda.all <- capscale(dist.count ~ Origin + Treatment + Origin:Treatment, metaData)
RsquareAdj(rda.all)
anova.cca(rda.all)
var.all <- sum(data.frame(anova.cca(rda.all))[,2])

# Marginal effect
anova.cca(rda.all, by="margin", scope="Origin")
anova.cca(rda.all, by="margin", scope="Treatment")
anova.cca(rda.all, by="margin", scope="Origin:Treatment")

# Variation partitioning to extract adjusted R^2

factor.mt <- model.matrix(~ 0+Origin*Treatment, metaData)
pop <- factor.mt[,1:3]
trt <- factor.mt[,4:5]
interaction <- factor.mt[,6:9]
res.varpart <- varpart(t(assay(count.vst)), pop, trt, interaction)
res.varpart

plot(res.varpart)


## Principal Component Analyis - PCA ----
pca.df <-  plotPCA(count.vst,
                   intgroup = c("Treatment", "Origin"),
                   ntop = nrow(ddsObj), returnData = TRUE)

percentVar <- round(100 * attr(pca.df, "percentVar"), 2)

pca.plot <- ggplot(pca.df, aes(PC1, PC2, color = Treatment, shape = Origin)) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", col = "grey") +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_bw() +
  labs(shape = "Origin") +
  labs(colour = "Treatment") +
  # scale_colour_manual(values = c("#006ddb", "#db6d00", "#920000"),
  scale_colour_manual(values = c("#117733", "#db6d00", "#920000"),
                      labels = c("2\u00B0C", "6\u00B0C", "10\u00B0C")) +
  scale_shape_manual(values = c(8, 19, 17)) +
  theme(
    axis.text = element_text(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank()
  )
pca.plot

# ggsave("02_Results/PCA_Temperature.pdf", pca.plot, width = 4, height = 2.5, scale = 1.5)


## Plot RDA ----
rda.s <- scores(rda.all, display="site")
rda.df<- data.frame(CAP1 = rda.s[,1],
                    CAP2 = rda.s[,2],
                    Origin = metaData$Origin,
                    Treatment = metaData$Treatment)

eigen.val <- rda.all$CCA$eig
rda1 <- round((eigen.val[1]/sum(eigen.val))*100,2)
rda2 <- round((eigen.val[2]/sum(eigen.val))*100,2)
percentVar.rda <- c(rda1, rda2)

paxis <- data.frame(anova.cca(rda.all, by="axis"))

rda.plot <- ggplot(rda.df, aes(CAP1, CAP2, color = Treatment, shape = Origin)) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", col = "grey") +
  geom_point(size = 4) +
  xlab(paste0("db-RDA1: ", percentVar.rda[1], "%")) +
  ylab(paste0("db-RDA2: ", percentVar.rda[2], "%")) +
  theme_bw() +
  labs(shape = "Origin") +
  labs(colour = "Temperature") +
  scale_colour_manual(values = c("#117733", "#db6d00", "#920000"), 
                      labels = c("2\u00B0C", "6\u00B0C", "10\u00B0C")) +
  scale_shape_manual(values = c(8, 19, 17)) +
  theme(
    axis.text = element_text(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank()
  )
rda.plot

# DETs identification ----------------------------------------------------------

# Cutoff for ajusted p-value:
padj.cutoff <- 0.001

# DETs according to populations

design(ddsObj) <- formula(~ Origin)
dds.pop <- DESeq(ddsObj, test = "LRT", reduced = ~1)
res.pop <- results(dds.pop, alpha = padj.cutoff)
res.pop.sig <- subset(res.pop, padj < padj.cutoff)

write.csv(res.pop, "02_Results/DESeq2_LRT_Origin.csv")

# DETs according to temperature as categorical factor

design(ddsObj) <- formula(~ Treatment)
dds.temp <- DESeq(ddsObj, test = "LRT", reduced = ~1)
res.temp <- results(dds.temp, alpha = padj.cutoff)
res.temp.sig <- subset(res.temp, padj < padj.cutoff)

write.csv(res.temp, "02_Results/DESeq2_LRT_Treatment.csv")

# DETs according to Origin x Temperature interaction

design(ddsObj) <- formula(~ Origin + Treatment + Origin:Treatment)
dds.pxt <- DESeq(ddsObj, test = "LRT", reduced = ~ Origin + Treatment)
res.pxt <- results(dds.pxt, alpha = padj.cutoff)
res.pxt.sig <- subset(res.pxt, padj < padj.cutoff)

write.csv(res.pxt, "02_Results/DESeq2_LRT_OxT.csv")

## Venn diagram ----

s1 <- list(
  "Treatment" = row.names(res.temp.sig),
  "Origin" = row.names(res.pop.sig),
  "Origin:Treatment" = row.names(res.pxt.sig)
)
venn.plot <- plot(euler(s1, shape = "circle"),
  quantities = T,
  fill = c("gray88", "gray50", "red"),
  alpha = 0.7, edges = F,
  adjust_labels = T
)
venn.plot

## ggsave("02_Results/Venn_DETs_pxt.pdf", venn.plot, width = 3, height = 2, scale = 1.5)


# Module of co-expression ------------------------------------------------------

options(stringsAsFactors = FALSE)

count.vst.mat <- assay(count.vst)

## 1) For DETs among Origin ----

#----- Select expresssion data from normalized count
res.1.pop <- res.pop
res.1.pop.sig <- res.pop.sig

count.pop <- subset(
  count.vst.mat,
  row.names(count.vst.mat) %in% rownames(res.1.pop.sig)
)
datExpr <- t(count.pop)
names(datExpr) <- rownames(datExpr)

datTraits <- model.matrix(~ 0 + metaData$Origin)
rownames(datTraits) <- metaData$ID
names(datTraits) <- levels(factor(metaData$Origin))

sampleTree <- hclust(dist(datExpr), method = "average")
traitColors <- numbers2colors(datTraits, signed = FALSE)
plotDendroAndColors(sampleTree, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

#----- Soft-thresholding powers / network topology analysis function
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1, 2))
cex1 <- 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2",
     type = "n", main = paste("Scale independence")
)
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red") # this line corresponds to using an R^2 cut-off of h

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")

#----- One-step network construction and module detection
net <- blockwiseModules(datExpr,
                        power = 8,
                        TOMType = "unsigned", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = (1 - 0.75), # Merge modules with the criteria: one minus their correlation
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = F,
                        verbose = 3
)

# number of module detected:
table(net$colors)

write.table(data.frame(net$colors),
            file = "02_Results/Module_Population.txt",
            quote = F)


# Module assignment and module eigengene information
moduleLabels <- paste("_", net$colors, sep="")
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs

# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MEs, "Eigengene dendrogram - Origin", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)

# Expression per module

mod.express <- datExpr
identical(row.names(mod.express), metaData$ID)
identical(colnames(mod.express), names(net$colors))

row.names(mod.express) <- metaData$Origin
colnames(mod.express) <- paste("O_", net$colors, sep="")

n.DETs.pop <- data.frame(Name = colnames(mod.express)) %>% 
  group_by(Name) %>% summarise(count = n())
n.DETs.pop$NewName <- "NA"
n.DETs.pop[1,3] <- "ME_0"
for(i in 2:nrow(n.DETs.pop)){
  n.DETs.pop[i,3] <- paste(n.DETs.pop[i,1], " (",
                           n.DETs.pop[i,2], ")", sep = "")
}


mod.eig <- moduleEigengenes(mod.express, net$colors)
mod.express1 <- mod.eig$eigengenes
mod.express1$Origin <- metaData$Origin

pop.express <- reshape2::melt(mod.express1)
colnames(pop.express) <- c("Origin", "Module", "Expression")
levels(pop.express$Module) <- paste("O_", levels(as.factor(net$colors)), sep="")
levels(pop.express$Module) <- n.DETs.pop$NewName
pop.express$Origin <- factor(pop.express$Origin, levels = c("NNC", "SLE", "ESS"))


express.pop.plot <- ggplot(subset(pop.express, !Module == "ME_0"), 
                           aes(x = Origin, y = scale(Expression))) +
  stat_summary(fun.data = mean_se, geom = "pointrange") +
  xlab(paste0("Origin")) +
  ylab("") +
  facet_wrap(~ Module, scales = "free", ncol = 3) + 
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_text(hjust = 0)
  )  +
  ggtitle("Origin")

express.pop.plot


## 2) For DETs among Treatments ----

#----- Select expresssion data from normalized count
res.1.trt <- res.temp
res.1.trt.sig <- res.temp.sig

count.temp <- subset(
  count.vst.mat,
  row.names(count.vst.mat) %in% rownames(res.1.trt.sig)
)
datExpr <- t(count.temp)
names(datExpr) <- rownames(datExpr)

datTraits <- model.matrix(~ 0 + metaData$Treatment)
rownames(datTraits) <- metaData$ID
names(datTraits) <- c("02C", "06C", "10C")

sampleTree <- hclust(dist(datExpr), method = "average")
traitColors <- numbers2colors(datTraits, signed = FALSE)
plotDendroAndColors(sampleTree, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap"
)

#----- Soft-thresholding powers / network topology analysis function
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1, 2))
cex1 <- 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2",
     type = "n", main = paste("Scale independence")
)
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red") # this line corresponds to using an R^2 cut-off of h

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity")
)
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")

#----- One-step network construction and module detection
net <- blockwiseModules(datExpr,
                        power = 7,
                        TOMType = "unsigned", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = (1 - 0.75), # Merge modules with the criteria: one minus their correlation
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = F,
                        verbose = 3
)

# number of module detected:
table(net$colors)

write.table(data.frame(net$colors),
            file = "02_Results/Module_Treatment.txt",
            quote = F)

# Module assignment and module eigengene information
moduleLabels <- paste("_", net$colors, sep="")
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs

# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MEs, "Eigengene dendrogram - Treatment", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)

# Expression per module

mod.express <- datExpr
identical(row.names(mod.express), metaData$ID)
identical(colnames(mod.express), names(net$colors))

row.names(mod.express) <- metaData$Treatment
colnames(mod.express) <- paste("T_", net$colors, sep="")

n.DETs.trt <- data.frame(Name = colnames(mod.express)) %>% 
  group_by(Name) %>% summarise(count = n())
n.DETs.trt$NewName <- "NA"
n.DETs.trt[1,3] <- "ME_0"
for(i in 2:nrow(n.DETs.trt)){
  n.DETs.trt[i,3] <- paste(n.DETs.trt[i,1], " (",
                           n.DETs.trt[i,2], ")", sep = "")
}


mod.eig <- moduleEigengenes(mod.express, net$colors)
mod.express1 <- mod.eig$eigengenes
mod.express1$Treatment <- metaData$Treatment

trt.express <- reshape2::melt(mod.express1)
colnames(trt.express) <- c("Treatment", "Module", "Expression")
levels(trt.express$Module) <- paste("O_", levels(as.factor(net$colors)), sep="")
levels(trt.express$Module) <- n.DETs.trt$NewName
trt.express$Treatment <- factor(trt.express$Treatment, levels = c("02C", "06C", "10C"))


express.trt.plot <- ggplot(subset(trt.express, !Module == "ME_0"), 
                           aes(x = Treatment, y = scale(Expression))) +
  stat_summary(fun.data = mean_se, geom = "pointrange") +
  xlab("Temperature (\u00B0C)") +
  ylab(paste0("Eigengenes")) +
  scale_x_discrete(labels = c("2", "6", "10")) +
  facet_wrap(~ Module, scales = "free", ncol = 1) + 
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(color = "black"),
        strip.background = element_blank(),
        strip.text.x = element_text(hjust = 0)
  ) +
  ggtitle("Temperature") 

express.trt.plot


# Functional analyses ----------------------------------------------------------

# list of genes/transcripts
annot <- read_trinotate("00_Data/trinotate_annotation_report.tsv")
all.genes <- row.names(ddsObj)
pop.genes <- row.names(res.pop.sig)
temp.genes <- row.names(res.temp.sig)
pxt.genes <- row.names(res.pxt.sig)

# Create a list with for each gene, containing vectors with all terms for each gene

GO.pfam <- split_GO(annot, hit = "gene_ontology_pfam")
pfam <- split_pfam(annot)
gene.go <- data.frame(Gene.stable.ID = GO.pfam$transcript,
                      GO.term.accession = GO.pfam$go)
gene2GO <- tapply(gene.go$GO.term.accession, gene.go$Gene.stable.ID, function(x)x)

## Enrichment analyses ----
go2gene.total <- NULL

### For Origin effect ----

res.GO.enrich.pop <- mat.or.vec(0,0)

# Define vector that is 1 if gene belong to a given meta-modules and 0 otherwise
tmp <- ifelse(all.genes %in% pop.genes, 1, 0)
geneList <- tmp

# geneList needs names that match those for GO terms,
names(geneList) <- unlist(all.genes, function(x)x[1])

ontology <- c("BP", "CC", "MF")

for (i in 1:3) {  
  # Create topGOdata object:
  GOdata <- new("topGOdata",
                ontology = ontology[i],
                allGenes = geneList,
                geneSelectionFun = function(x)(x == 1),
                annot = annFUN.gene2GO, gene2GO = gene2GO)
  
  # Run Fisher’s Exact Test:
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  tab <- GenTable(GOdata, raw.p.value = resultFisher, topNodes = length(resultFisher@score),
                  numChar = 1000)
  tab$padj <- p.adjust(tab$raw.p.value, method = "fdr")
  
  # Retrieve transcript id related to a GO.ID
  allGO = genesInTerm(GOdata)
  go2gene.total <- c(go2gene.total, allGO)
  
  ## - create table results
  rm(temp)
  temp <- tab
  #  temp <- subset(tab, padj <= 0.1)
  temp$ontology <- ontology[i]
  
  ## res.GO.enrich.pop <- mat.or.vec(0,0)
  res.GO.enrich.pop <- rbind(res.GO.enrich.pop, temp)
}

write.csv(res.GO.enrich.pop,
          "02_Results/Results_GO_Enrich_Origin.csv")

# Figure for annotation analyses

df.go.pop <- subset(res.GO.enrich.pop, as.numeric(padj) <= 0.05)

df.go.pop$padj_log <- log10(as.numeric(df.go.pop$padj)) * (-1)

df.go.pop <- subset(df.go.pop, padj_log >= 0.1) %>% 
  arrange(ontology, padj_log) %>%
  mutate(order = row_number())

plot.go.pop <- ggplot(df.go.pop, aes(order, padj_log, fill=ontology)) + 
  geom_col()  + 
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Dark2")[2:3]) +
  xlab("") +
  ylab(expression(paste("-log10(", italic("P"),"-value)"))) +
  labs(fill = "Ontology") +
  coord_flip() +
  theme_classic()+
  theme(
    strip.background = element_rect(fill= "gray100"),
    strip.text.y = element_text(angle = 180+90, size = 12),
    legend.position = "right",
    axis.text = element_text(colour="black")
  ) +
  scale_x_continuous(
    breaks = df.go.pop$order,
    labels = df.go.pop$Term,
    expand = c(0,0)
  )

plot.go.pop


### For Temperature effect ----

res.GO.enrich.temp <- mat.or.vec(0,0)

# Define vector that is 1 if gene belong to a given meta-modules and 0 otherwise
tmp <- ifelse(all.genes %in% temp.genes, 1, 0)
geneList <- tmp

# geneList needs names that match those for GO terms,
names(geneList) <- unlist(all.genes, function(x)x[1])

ontology <- c("BP", "CC", "MF")

for (i in 1:3) {  
  # Create topGOdata object:
  GOdata <- new("topGOdata",
                ontology = ontology[i],
                allGenes = geneList,
                geneSelectionFun = function(x)(x == 1),
                annot = annFUN.gene2GO, gene2GO = gene2GO)
  
  # Run Fisher’s Exact Test:
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  tab <- GenTable(GOdata, raw.p.value = resultFisher, topNodes = length(resultFisher@score),
                  numChar = 120)
  tab$padj <- p.adjust(tab$raw.p.value, method = "fdr")
  
  # Retrieve transcript id related to a GO.ID
  allGO = genesInTerm(GOdata)
  go2gene.total <- c(go2gene.total, allGO)
  
  ## - creat table restults
  rm(temp)
  temp <- tab
  #  temp <- subset(tab, padj <= 0.1)
  temp$ontology <- ontology[i]
  
  ## res.GO.enrich.temp <- mat.or.vec(0,0)
  res.GO.enrich.temp <- rbind(res.GO.enrich.temp, temp)
}

write.csv(res.GO.enrich.temp,
          "02_Results/Results_GO_Enrich_Treatment.csv")

# Figure for annotation analyses

df.go.temp <- subset(res.GO.enrich.temp, as.numeric(padj) <= 0.05)

df.go.temp$padj_log <- log10(as.numeric(df.go.temp$padj)) * (-1)

df.go.temp <- subset(df.go.temp, padj_log >= 0.1) %>% 
  arrange(ontology, padj_log) %>%
  mutate(order = row_number())

plot.go.temp <- ggplot(df.go.temp, aes(order, padj_log, fill=ontology)) + 
  geom_col()  + 
  scale_fill_brewer(palette = "Dark2") + 
  xlab("") +
  ylab(expression(paste("-log10(", italic("P"),"-value)"))) +
  labs(fill = "Ontology") +
  coord_flip() +
  theme_classic()+
  theme(
    strip.background = element_rect(fill= "gray100"),
    strip.text.y = element_text(angle = 180+90, size = 12),
    legend.position = "right",
    axis.text = element_text(colour="black")
  ) +
  scale_x_continuous(
    breaks = df.go.temp$order,
    labels = df.go.temp$Term,
    expand = c(0,0)
  )

plot.go.temp

# ### For Population x Temperature interaction ----
# 
# res.GO.enrich.pxt <- mat.or.vec(0,0)
# 
# # Define vector that is 1 if gene belong to a given meta-modules and 0 otherwise
# tmp <- ifelse(all.genes %in% pxt.genes, 1, 0)
# geneList <- tmp
# 
# # geneList needs names that match those for GO terms,
# names(geneList) <- unlist(all.genes, function(x)x[1])
# 
# ontology <- c("BP", "CC", "MF")
# 
# for (i in 1:3) {  
#   # Create topGOdata object:
#   GOdata <- new("topGOdata",
#                 ontology = ontology[i],
#                 allGenes = geneList,
#                 geneSelectionFun = function(x)(x == 1),
#                 annot = annFUN.gene2GO, gene2GO = gene2GO)
#   
#   # Run Fisher’s Exact Test:
#   resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
#   tab <- GenTable(GOdata, raw.p.value = resultFisher, topNodes = length(resultFisher@score),
#                   numChar = 120)
#   tab$padj <- p.adjust(tab$raw.p.value, method = "fdr")
#   
#   # Retrieve transcript id related to a GO.ID
#   allGO = genesInTerm(GOdata)
#   go2gene.total <- c(go2gene.total, allGO)
#   
#   ## - creat table restults
#   rm(temp)
#   temp <- tab
#   #  temp <- subset(tab, padj <= 0.1)
#   temp$ontology <- ontology[i]
#   
#   ## res.GO.enrich.pxt<- mat.or.vec(0,0)
#   res.GO.enrich.pxt<- rbind(res.GO.enrich.pxt, temp)
# }
# 
# write.csv(res.GO.enrich.pxt,
#           "02_Results/Results_GO_Enrich_PxT.csv")
# 
# # Figure for annotation analyses
# 
# df.go.pxt <- subset(res.GO.enrich.pxt, as.numeric(padj) <= 0.05)
# df.go.pxt
# ### no functions detected for population x temperature interaction ###
# 
# 
# ### For ME_2 treatment: genes that are upregulated in warmer temperature ----
# 
# res.GO.enrich.mod2 <- mat.or.vec(0,0)
# 
# # Define vector that is 1 if gene belong to a given meta-modules and 0 otherwise
# mod.temp <- read.table("02_Results/Module_Treatment.txt")
# mod.temp.2 <- row.names(subset(mod.temp, net.colors == 2))
# 
# tmp <- ifelse(all.genes %in% mod.temp.2, 1, 0)
# geneList <- tmp
# 
# # geneList needs names that match those for GO terms,
# names(geneList) <- unlist(all.genes, function(x)x[1])
# 
# ontology <- c("BP", "CC", "MF")
# 
# for (i in 1:3) {  
#   # Create topGOdata object:
#   GOdata <- new("topGOdata",
#                 ontology = ontology[i],
#                 allGenes = geneList,
#                 geneSelectionFun = function(x)(x == 1),
#                 annot = annFUN.gene2GO, gene2GO = gene2GO)
#   
#   # Run Fisher’s Exact Test:
#   resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
#   tab <- GenTable(GOdata, raw.p.value = resultFisher, topNodes = length(resultFisher@score),
#                   numChar = 120)
#   tab$padj <- p.adjust(tab$raw.p.value, method = "fdr")
#   
#   # Retrieve transcript id related to a GO.ID
#   allGO = genesInTerm(GOdata)
#   go2gene.total <- c(go2gene.total, allGO)
#   
#   ## - creat table restults
#   rm(temp)
#   temp <- tab
#   #  temp <- subset(tab, padj <= 0.1)
#   temp$ontology <- ontology[i]
#   
#   ## res.GO.enrich.pxt<- mat.or.vec(0,0)
#   res.GO.enrich.mod2 <- rbind(res.GO.enrich.mod2, temp)
# }
# 
# write.csv(res.GO.enrich.mod2,
#           "02_Results/Results_GO_Enrich_Module2_trt.csv")
# 
# # Figure for annotation analyses
# 
# df.go.mod2<- subset(res.GO.enrich.mod2, as.numeric(padj) <= 0.05)
# ### no functions detected for module 2 ###
# 
# ### For ME_1 and 3 treatment: genes that are upregulated in warmer temperature ----
# 
# res.GO.enrich.mod13 <- mat.or.vec(0,0)
# 
# # Define vector that is 1 if gene belong to a given meta-modules and 0 otherwise
# mod.temp <- read.table("02_Results/Module_Treatment.txt")
# mod.temp.2 <- row.names(subset(mod.temp, net.colors %in% c(1,3)))
# 
# tmp <- ifelse(all.genes %in% mod.temp.2, 1, 0)
# geneList <- tmp
# 
# # geneList needs names that match those for GO terms,
# names(geneList) <- unlist(all.genes, function(x)x[1])
# 
# ontology <- c("BP", "CC", "MF")
# 
# for (i in 1:3) {  
#   # Create topGOdata object:
#   GOdata <- new("topGOdata",
#                 ontology = ontology[i],
#                 allGenes = geneList,
#                 geneSelectionFun = function(x)(x == 1),
#                 annot = annFUN.gene2GO, gene2GO = gene2GO)
#   
#   # Run Fisher’s Exact Test:
#   resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
#   tab <- GenTable(GOdata, raw.p.value = resultFisher, topNodes = length(resultFisher@score),
#                   numChar = 120)
#   tab$padj <- p.adjust(tab$raw.p.value, method = "fdr")
#   
#   # Retrieve transcript id related to a GO.ID
#   allGO = genesInTerm(GOdata)
#   go2gene.total <- c(go2gene.total, allGO)
#   
#   ## - creat table restults
#   rm(temp)
#   temp <- tab
#   #  temp <- subset(tab, padj <= 0.1)
#   temp$ontology <- ontology[i]
#   
#   ## res.GO.enrich.pxt<- mat.or.vec(0,0)
#   res.GO.enrich.mod13 <- rbind(res.GO.enrich.mod13, temp)
# }
# 
# write.csv(res.GO.enrich.mod13,
#           "02_Results/Results_GO_Enrich_Module13_trt.csv")
# 
# # Figure for annotation analyses
# 
# df.go.mod13<- subset(res.GO.enrich.mod13, as.numeric(padj) <= 0.05)
# 
# df.go.mod13$padj_log <- log10(as.numeric(df.go.mod13$padj)) * (-1)
# 
# df.go.mod13 <- subset(df.go.mod13, padj_log >= 0.1) %>% 
#   arrange(ontology, padj_log) %>%
#   mutate(order = row_number())
# 
# plot.go.mod13 <- ggplot(df.go.mod13, aes(order, padj_log, fill=ontology)) + 
#   geom_col()  + 
#   scale_fill_brewer(palette = "Dark2") + 
#   xlab("") +
#   ylab(expression(paste("-log10(", italic("P"),"-value)"))) +
#   labs(fill = "Ontology") +
#   coord_flip() +
#   theme_classic()+
#   theme(
#     strip.background = element_rect(fill= "gray100"),
#     strip.text.y = element_text(angle = 180+90, size = 12),
#     legend.position = "right",
#     axis.text = element_text(colour="black")
#   ) +
#   scale_x_continuous(
#     breaks = df.go.mod13$order,
#     labels = df.go.mod13$Term,
#     expand = c(0,0)
#   )
# 
# plot.go.mod13


# Modules of co-expression and Functionnal analyses ----------------------------


### Figures GO-MOdules 

#### 1- Plot significant GO terms

enrich.pop <- subset(res.GO.enrich.pop, padj <= 0.05) %>%
  mutate(group = "Origin")

enrich.temp <- subset(res.GO.enrich.temp, padj <= 0.05) %>%
  mutate(group = "Temperature")

enrich.all <- rbind(enrich.pop, enrich.temp)
enrich.all$padj_log <- log10(enrich.all$padj)*(-1)
enrich.all$ontology <- as.factor(enrich.all$ontology)
levels(enrich.all$ontology) <- c("Biological Process",
                                 "Cellular Component",
                                 "Molecular Function")

go.subset <- enrich.all %>% 
  arrange(ontology, group, padj_log) %>%
  mutate(order = row_number())

plot.bp.1 <- ggplot(go.subset, aes(order, padj_log, fill=group)) + 
  geom_col()  + 
  scale_fill_brewer(palette = "Set1") +
  xlab("") +
  ylab(expression(paste("-log10(", italic("P"),"-value)"))) +
  labs(fill = "DETs \nacross:") +
  coord_flip() +
  theme_classic() +
  facet_grid(ontology ~ .,
             scales = "free_y",
             space = "free",
             switch = "both") +
  theme(
    strip.text.x = element_text(hjust = 0),
    strip.background = element_rect(fill = "gray85", colour = "gray85"),
    strip.text.y = element_text(size = 10),
    strip.placement = "top",
    legend.position = "right",
    axis.text = element_text(colour="black")
  ) +
  scale_x_continuous(
    breaks = go.subset$order,
    labels = go.subset$Term,
    expand = c(0,0)
  ) 

plot.bp.1


#### 2- Plot for associated modules of co-expressed genes

mod.pop <- read.table("02_Results/Module_Population.txt")
mod.temp <- read.table("02_Results/Module_Treatment.txt")

mod.pop$transcript <- row.names(mod.pop)
mod.temp$transcript <- row.names(mod.temp)

mod.pop$group <- "Origin"
mod.temp$group <- "Temperature"

modules.all <- rbind(mod.pop, mod.temp)
modules.all$net.colors <- paste(modules.all$group, modules.all$net.colors, sep="-")
### Significant GO_ID

list.go <- go2gene.total[enrich.all$GO.ID]

list.det <- lapply(list.go,function(x) x[x %in% c(pop.genes, temp.genes)] )

df.go.det <- mat.or.vec(0,0)

for(i in 1:length(list.det)) {
  go.id <- list.det[i]
  temp <- data.frame(GO = names(go.id),
                     transcript = go.id)
  colnames(temp) <- c("GO", "transcript")
  df.go.det <- rbind(df.go.det, temp)
}

df.go.det <- left_join(df.go.det, modules.all)

modules.go.df <- left_join(df.go.det, go.subset, by = c("GO" = "GO.ID", 
                                                        "group" = "group"))
mod.subset <- subset(modules.go.df, !net.colors %in% c("Origin-0", "Temperature-0"))
mod.subset <- subset(mod.subset, !Term == "NA")

mod.subset$net.colors <- factor(mod.subset$net.colors,
                                levels = c("Origin", 
                                           levels(factor(mod.subset$net.colors))[1:6],
                                           "", "Temperature",
                                           levels(factor(mod.subset$net.colors))[7:9]))

df.onto.mod <- mod.subset %>% group_by(GO, net.colors, Term, ontology) %>% 
  summarise(Number = n())
df.onto.mod$DETs <- substr(df.onto.mod$net.colors, 1, 3)
df.prop <- df.onto.mod %>% group_by(GO, DETs, ontology, Term) %>% 
  summarise(net.colors = net.colors,
            Prop = prop.table(Number))

plot.bp.2 <- ggplot(mod.subset, aes(order)) +
  geom_bar(aes(fill = net.colors), 
           position = "fill") +
  scale_fill_manual(values = c("white",
                               brewer.pal(9, "Spectral")[1:6],
                               "white", "white",
                               brewer.pal(9, "Spectral")[7:9]),
                    drop = FALSE,
                    labels = c("Origin", 
                               paste("O_", 1:6, sep = ""),
                               "",
                               "Temperature", 
                               paste("T_", 1:3, sep = ""))) + 
  xlab("") +
  ylab("Proportion of transcripts") +
  labs(fill = "Module of \nco-expression") +
  coord_flip() +
  theme_classic()+
  scale_x_continuous(
    breaks = mod.subset$order,
    labels = mod.subset$Term,
    expand = c(0,0)
  ) +
  facet_grid(ontology ~ .,
             scales = "free_y",
             space = "free",
             switch = "both") +
  theme(
    strip.text = element_blank(),
    legend.position = "right",
    axis.text = element_text(colour="black"),
    axis.text.y = element_blank()
  ) 
plot.bp.2


fig.enrich <- ggarrange(plot.bp.1, plot.bp.2, 
                        ncol = 2,
                        widths = c(2, 1),
                        labels = c(
                          "A. Gene ontology (GO) term enrichment analysis",
                          "B. Distributiono of co-expressed transcripts modules"
                        ),
                        font.label = list(size = 11),
                        hjust = 0, vjust = -0.5
) +
  theme(plot.margin = margin(.7, 0.1, 0.1, 0.1, "cm"))
fig.enrich


# Disentangling origin effect --------------------------------------------------

## Genetic variation -----------------------------------------------------------

### Data ----
vcf.all <- read.vcfR("00_Data/populations.22227snps.94ind.n2HW.single.recode.vcf")
gl.all  <- vcfR::vcfR2genlight(vcf.all) 

temp.ID <- metaData$ID_GQ
meta.temp <- subset(metaData, !ID_GQ == "TPB085") # bad sequencing quality for genomics

# genlight object
gl.temp  <- gl.all[indNames(gl.all) %in% temp.ID]
gl.temp <- gl.temp[match(meta.temp$ID_GQ, indNames(gl.temp))]

## Genotype matrix according to the number of alternative allele
gen.temp <- as.matrix(gl.temp)

## Change sample names to match the rna-seq data
rownames(gen.temp) <- plyr::mapvalues(rownames(gen.temp), 
                                      from = metaData$ID_GQ, 
                                      to = metaData$ID)

## Re-order according to metaData
gen.temp <- gen.temp[match(meta.temp$ID, row.names(gen.temp)),]

identical(rownames(gen.temp), meta.temp$ID)


identical(indNames(gl.temp), meta.temp$ID_GQ)
pop(gl.temp) <- meta.temp$Origin
pop.fst <- gl.fst.pop(gl.temp, nboots = 100, percent = 95, nclusters = 20)

fst <- pmax(pop.fst$Fsts, 0)
diag(fst) <- 0
fst.matrix <- as.matrix(as.dist(fst))

pvalue <- pop.fst$Pvalues
diag(pvalue) <- 0
pvalue <- as.dist(pvalue)
pvalue <- as.matrix(pvalue)
diag(pvalue) <- 1

stars <- cut(pvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), 
             label=c("***", "**", "*", ""))  

stars.matrix <- matrix(stars, ncol = 3)

rownames(stars.matrix) <- rownames(fst.matrix)
colnames(stars.matrix) <- colnames(fst.matrix)


fst <- reshape2::melt(fst.matrix)
pvalue <- reshape2::melt(stars.matrix)
co <- cbind(fst, pvalue$value)
colnames(co) <- c("var1", "var2", "Fst", "pvalue")

p <- ggplot(aes(x=var1, y=var2, fill=Fst), data=co) +
  geom_tile(col = "black") + 
  scale_fill_gradient2(low="#2C7BB6", mid="white", high="#D7191C") + 
  geom_text(aes(label=stars), color="black", size=6) + 
  labs(y=NULL, x=NULL, fill="Fst") + 
  # geom_vline(xintercept=c(0.5, 1.5, 2.5, 3.5), size=.5, color="black") +
  # geom_hline(yintercept=c(0.5, 1.5, 2.5, 3.5), size=.5, color="black") +
  theme_classic() + theme(axis.text = element_text(angle = 0, hjust = 0.5, size = 12),
                          axis.line = element_blank(),
                          axis.ticks = element_blank(),
                          plot.margin = margin(0.5, .5, .7, 0.5, "cm"))

dend <- hclust(as.dist(fst.matrix))

data <- dendro_data(dend, type = "rectangle")
dendro <- ggplot(segment(data)) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), size = 1.5) +
  theme_classic() + 
  theme(line = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(), 
        plot.margin = margin(t = 0.5, r = 4, b = 0, l = 2.5, "cm")) +
  scale_x_reverse()

plot.fst <- ggarrange(dendro, p, heights = c(.3,1), 
                   ncol = 1)

plot.fst




# ## dartR
# 
gl.temp@other$loc.metrics.flags$monomorphs <- FALSE
gl.temp@ploidy <- rep(as.integer(2), nInd(gl.temp))

he.temp <- gl.report.heterozygosity(gl.temp)

gi.temp <- gl2gi(gl.temp)
summary.loci <- summary(gi.temp)

NNC.loc <- poppr::locus_table(gi.temp, pop = "NNC")
SS.loc <- poppr::locus_table(gi.temp, pop = "ESS")
SLE.loc <- poppr::locus_table(gi.temp, pop = "SLE")

he.pop <- data.frame(Population = rep(c("ESS", "SLE", "NNC"), each = nrow(NNC.loc)),
                     Hexp = c(data.frame(SS.loc)$Hexp,
                              data.frame(SLE.loc)$Hexp,
                              data.frame(NNC.loc)$Hexp))

TukeyHSD(aov(Hexp ~ Population, he.pop))

he.pop$Population <- factor(he.pop$Population, levels = c("NNC", "SLE", "ESS"))

mean.he <- ggplot(he.pop, aes(x = Population, y = Hexp, shape = Population)) +
  stat_summary(fun.data = "mean_se", size = 1, stroke = 1.1) +
  scale_shape_manual(values = c(8, 19, 17)) +
  theme_classic() +
  ylab(expression(Gene~diversity~(italic(H[E]))))  + xlab("") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12),
        plot.margin = margin(2,.5,.3,.5, "cm")) +
  annotate(geom = "text", x = c("NNC", "SLE", "ESS"), y = 0.225,
           label = c("a", "b", "a"), size = 5)
mean.he

he.pop %>% group_by(Population) %>% summarise(mean = mean(Hexp))


## Sampling environmental conditions --------------------------------------------

### Environmental data before sampling ----
env.all <- read_csv("00_Data/Sampling_env_BNAM.csv") %>% 
  dplyr::select(-c(fid, Index, N)) %>%
  data.frame(.)
row.names(env.all) <- c("SLE", "ESS", "NNC")

### Remove correlated variables ----

env <- env.all[,c(6:41)]

cor.env <- cor(env)

library(ggcorrplot)
env.cor.plot <- ggcorrplot(cor.env, method = "circle", 
                           hc.order = F, type = "lower",
                           outline.col = "white", 
                           legend.title = expression(Correlation~italic(r)),
                           tl.cex = 18,
                           show.diag = T) +
  theme(legend.title = element_text(size=18),
        legend.text = element_text(size=15),
        legend.key.size = unit(1.5, "cm"))

env.cor.plot

cor_matrix_rm <- cor.env                 # Modify correlation matrix
cor_matrix_rm[upper.tri(cor_matrix_rm)] <- 0
diag(cor_matrix_rm) <- 0

env.new <- env[ , !apply(cor_matrix_rm,    # Remove highly correlated variables
                         2,
                         function(x) any(abs(x) > 0.90))]
env.new

## Select bottom instead of surface temperature for the two same months
## as more reliable to the species biology

env.new <- env[,c("T_Jul_Bott", "T_Jan_Bott")]

cor(env.new)

# Euclidean distance
vegdist(env.new[,1:2], method = "euclidean")

# PCA plot
pca.env <- rda(env.new, scale = T)

perc <- round(100*(summary(pca.env)$cont$importance[2, 1:2]), 2)

sc_si <- data.frame(scores(pca.env, display="sites", choices=c(1,2)))
sc_si$Origin <- factor(row.names(sc_si), levels = c("NNC", "SLE", "ESS"))
sc_sp <- data.frame(scores(pca.env, display="species", choices=c(1,2)))
rownames(sc_sp) <- c("July\nbottom T\u00B0C", "January\nbottom T\u00B0C")
sc_sp$Origin <- "NNC"

pca.enviro.plot <- ggplot(sc_si, aes(PC1, PC2, shape = Origin)) +
  geom_hline(yintercept = 0, linetype = "dashed", col = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", col = "grey") +
  geom_point(size = 5, stroke = 1.1) +
  geom_text(label = sc_si$Origin,
            vjust = c(-2.1, 0.5, 3),
            hjust = c(0.5, -1, 0.5)) +
  scale_shape_manual(values = c(8, 19, 17)) +
  geom_segment(data = sc_sp, aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow=arrow(length=unit(0.2,"cm")),
               alpha = 0.75, color = 'black') +
  xlab(paste0("PC1 (", perc[1], "%)")) +
  ylab(paste0("PC2 (", perc[2], "%)")) + 
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_blank(),
    legend.position = "none"
  ) +
  annotate("text", x = sc_sp$PC1, y = sc_sp$PC2,
           label = rownames(sc_sp), col = "black",
            vjust = c(0.5,0.2), hjust = c(1,1.2)) +
  ylim(-1.1, 1.5) + xlim(-1.2, 1.2)

pca.enviro.plot

## Shrimp physiological states -------------------------------------------------

df.mm <- read.table("00_Data/Mortality_Molting.txt", header = T)
df.mm$Prop_mort <- 1-(df.mm$n.final/df.mm$n.init)

df.mm$Treatment <- factor(df.mm$Treatment,
                                  levels = c("2C", "6C", "10C"))
df.mm$Origin <- factor(df.mm$Origin,
                               levels = c("NNC", "SLE", "ESS"))

aov.mort <- aov(asin(sqrt(Prop_mort)) ~ Treatment*Origin, df.mm)
anova(aov.mort)
TukeyHSD(aov.mort, "Treatment", group = TRUE, console = T)

aov.molt <- aov(asin(sqrt(molt_rate)) ~ Treatment*Origin, df.mm)
anova(aov.molt)
TukeyHSD(aov.molt)

df <- df.mm
df$Tank <- as.factor(df$Tank)
gpa_mixed = lmerTest::lmer(asin(sqrt(molt_rate)) ~ Treatment*Origin + (1 | Tank), data = df.mm)
anova(gpa_mixed)


df.mm$TxO <- paste(df.mm$Treatment, df.mm$Origin, sep = "_")
aov.molt2 <- aov(asin(sqrt(molt_rate)) ~ TxO, df.mm)
test.res <- agricolae::HSD.test(aov.molt2, "TxO", console = T)


mean.molt <- df.mm %>% group_by(Origin, Treatment) %>%
  summarise(mean_molt = mean(molt_rate),
            se_molt = sd(molt_rate) / sqrt(length(molt_rate)))
mean.molt

mean.mort <- df.mm %>% group_by(Treatment) %>%
  summarise(mean_mort = mean(Prop_mort),
            se_mort = sd(Prop_mort) / sqrt(length(Prop_mort)))
mean.mort


### plot physiology

scaleFUN <- function(x) sprintf("%.3f", x)

df.mm.melted <- reshape::melt(df.mm)
subset.melted <- subset(df.mm.melted, variable %in% c("Prop_mort", "molt_rate"))

subset.melted$variable <- factor(subset.melted$variable,
                               levels = c("Prop_mort", "molt_rate"))

subset.mortality <- subset(subset.melted, variable == "Prop_mort")

plot.mortality <- ggplot(subset.mortality) +
  geom_point(aes(x = Treatment, y = value,
                 color = Treatment, shape = Origin),
             position = position_jitter(0.1), size = 3) +
  scale_colour_manual(values = c("#117733", "#db6d00", "#920000"), 
                      labels = c("2", "6", "10")) +
  scale_shape_manual(values = c(8, 19, 17)) +
  guides(color = "none") +
  theme_classic() +
  scale_x_discrete(labels = c("2", "6", "10")) +
  scale_y_continuous(labels=scaleFUN, n.breaks = 4) +
  ylab("Mortality rate") + xlab("Temperature (\u00B0C)")  +
  annotate(geom = "text", x = c("2C", "6C", "10C"), y = 0.26,
           label = c("a", "ab", "b"), size = 4)

plot.mortality

subset.molting <- subset(subset.melted, variable == "molt_rate")

df.temp <- subset.molting[,c("Origin", "Treatment", "value")] %>%
  group_by(Origin, Treatment) %>% mutate(value = mean(value)) %>% unique()
tukey.res <- test.res$groups
tukey.res$Treatment <- str_split_fixed(row.names(tukey.res), "_", n = 2)[,1]
tukey.res$Origin <- str_split_fixed(row.names(tukey.res), "_", n = 2)[,2]
tukey.res <- left_join(tukey.res, df.temp)

plot.molting <- ggplot() +
  geom_point(data = subset.molting, aes(x = Treatment, y = value,
                 color = Treatment, shape = Origin),
             position = position_jitter(0.1), size = 3) +
  scale_colour_manual(values = c("#117733", "#db6d00", "#920000"), 
                      labels = c("2", "6", "10")) +
  scale_shape_manual(values = c(8, 19, 17)) +
  guides(color = "none") +
  theme_classic() +
  scale_x_discrete(labels = c("2", "6", "10")) +
  scale_y_continuous(labels=scaleFUN, n.breaks = 4) +
  ylab("Molting rate") + xlab("Temperature (\u00B0C)") +
  geom_text(data = tukey.res, aes(x = Treatment, y = value, 
                           label = groups), nudge_x = -0.3)

plot.molting

fig.physio <- ggarrange(plot.mortality, NULL, plot.molting,
                        nrow = 1,
                        labels = c("A.", "", "B."),
                        widths = c(1, 0.1, 1),
                        common.legend = T, legend = "right")
fig.physio


# RNA vs Genet vs Enviro -------------------------------------------------------

## RNA data
rna.mat <- t(assay(count.vst))
dist.rna.mat <- vegdist(rna.mat, method = "euclidean")

## Environment data
env.new$Origin <- row.names(env.new)
tempo <- left_join(metaData, env.new)
env.mat <- as.matrix(tempo[,grep("T_", colnames(tempo))])
row.names(env.mat) <- tempo$ID

# ## Genetic data
fst.pop <- pop.fst$Fsts
diag(fst.pop) <- 0
fst.pop[fst.pop < 0 ] <- 0
fst.pop <- as.dist(fst.pop)

fst.score <- data.frame(ape::pcoa(fst.pop)$vector)
fst.score$Origin <- row.names(fst.score)
tempo <- left_join(metaData, fst.score)
fst.mat <- as.matrix(tempo[,grep("Axis", colnames(tempo))])
row.names(fst.mat) <- tempo$ID

## Variation partitioning 
identical(rownames(rna.mat), rownames(env.mat))
identical(rownames(rna.mat), rownames(fst.mat))

Treatment <- model.matrix(~ 0 + metaData$Treatment)

dist.count <- vegdist(t(assay(count.vst)), method = "euclidean")


## model selection -----
### Model selection

mod0 <- capscale(dist.rna.mat ~ 1)
mod1 <- capscale(dist.rna.mat ~ Treatment + env.mat + Treatment:env.mat +
                   Treatment + fst.mat + Treatment:fst.mat, scale = T)

# res <- ordistep(mod0, mod1, direction = "both", permutations = 999)
# res$anova

res <- ordiR2step(mod0, mod1, permutations = 999)
res$anova

sel.osR2_adj <- res
sel.osR2_adj$anova$`Pr(>F)` <- p.adjust (res$anova$`Pr(>F)`, method = 'holm', n = 6)
sel.osR2_adj$anova


### Var partition

rda.var.part <- capscale(dist.rna.mat ~ Treatment + env.mat + fst.mat)
anova.cca(rda.var.part)
anova.cca(rda.var.part, by = "margin", scope = "Treatment")
anova.cca(rda.var.part, by = "margin", scope = "env.mat")
anova.cca(rda.var.part, by = "margin", scope = "fst.mat")

mm <- model.matrix(~ Treatment + env.mat + fst.mat)[,-1]
colnames(mm)

trt <- mm[,1:3]
env <- mm[,4:5]
gen <- mm[,6]

res.varpart <- varpart(rna.mat, trt, env, gen)
res.varpart

##plot(res.varpart)

anova.cca(capscale(dist.rna.mat ~ Treatment))
anova.cca(capscale(dist.rna.mat ~ env.mat))
anova.cca(capscale(dist.rna.mat ~ fst.mat))



# Figures ----------------------------------------------------------------------

## Fig_1 sampling map
ggsave("02_Results/Fig1_Map.jpeg", map.fig, width = 8, height = 5, scale = 2,
       units = "cm")  


# Fig2 mortality and molting rate
ggsave("02_Results/Fig2_mort_molt_raw_data.jpg", 
       fig.physio + theme(plot.background = element_rect(fill = "white")),
       width = 8, height = 3.5, scale = 2,
       units = "cm")

## Fig_3 gene expression
rda.plot1 <- ggarrange(rda.plot, labels = "A.", hjust = 0.25) +
  theme(plot.margin = margin(0.5,0.5,0, 1, "cm")) 

fig.express <- ggarrange(express.trt.plot, NULL, express.pop.plot, NULL,
                         ncol = 4,
                         labels = c("B.", "",
                                    "C.", ""),
                         widths = c(1.2, 0.1, 3, 0.2),
                         hjust = c(0.25, -0.5, -1, -0.5)) +
  theme(plot.margin = margin(0.5,0.5,0.5,1, "cm")) 

fig3 <- ggarrange(rda.plot1, fig.express,
                  labels = "",
                  nrow = 2, heights = c(1, 1.5))
fig3
ggsave("02_Results/Fig3_Gene_expression.jpg", 
       fig3 + theme(plot.background = element_rect(fill = "white")), 
       width = 8.5, height = 8.5, scale = 2.1,
       units = "cm")

       ##width = 5, height = 5.5, scale = 1.5)  


# Fig4 genetic variation
fig.gen <- ggarrange(mean.he, plot.fst,
          ncol = 2, labels = c("A.", "B."),
          widths = c(1, 1.3),
          vjust = 3)
fig.gen

ggsave("02_Results/Fig4_genetic_diversity.jpg", 
       fig.gen + theme(plot.background = element_rect(fill = "white")),
       width = 8, height = 4, scale = 2,
       units = "cm")

# Fig5 Environment at origin
ggsave("02_Results/Fig5_PCA_env.jpg", pca.enviro.plot,
       width = 7, height = 3.5, scale = 2.3,
       units = "cm")


# Fig_S1 Gene expression and modules
ggsave("02_Results/Supp_Fig1_GO_Modules_all.pdf", fig.enrich,
       width = 8, height = 4, scale = 1.5)

# Fig_S2

library(gridExtra)

env.table <- env.all[,c(6:41)]

env.cor <- cor(env.table)

# Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }

  upper_tri <- get_upper_tri(env.cor)
  upper_tri

  # Melt the correlation matrix
  melted_cormat <- melt(upper_tri, na.rm = TRUE)

  # Heatmap
ggheatmap <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 11, hjust = 1))+
#    coord_fixed() + 
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))+ 
  scale_y_discrete(guide = guide_axis(position = "none"))

ggheatmap  


env.table <- round(env.all[,c(6:41)]) %>% format(nsmall = 2) %>% t() %>% data.frame()

env.reshap <- data.frame(Origin = rep(c("SLE","ESS","NNC"), each = nrow(env.table)),
                         Name = rep(rownames(env.table), 3),
                         Value = c(env.table$SLE, env.table$ESS, env.table$NNC))
env.reshap$Name <- factor(env.reshap$Name,
                          levels = rownames(env.table))

df.table <- ggplot(env.reshap, aes(x = Origin, y = Name,
                           label = Value)) +
  geom_text(size = 3.5, hjust = 1) + 
  theme_minimal() + 
  theme(panel.grid.major = element_blank(), legend.position = "none",
        panel.border = element_blank(), 
        axis.text.x =  element_text(colour = "black", size = 14, hjust = 1),
        axis.text.y = element_text(size = 11, hjust = 0),
        axis.ticks =  element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  scale_x_discrete(guide = guide_axis(position = "top")) +
  scale_y_discrete(guide = guide_axis(position = "right"))
df.table

fig.cor.env <- ggarrange(df.table + 
                           theme(plot.margin = margin(r = 1, l = 36)),
                         ggheatmap +
                           theme(plot.margin = margin(r = 12, l = 0)), 
                         ncol = 2,
          widths = c(1,2.2),
          labels = c("A.","B."),
       #   hjust = c(0, -2),
          align = "h")
fig.cor.env

ggsave("02_Results/Supp_Fig2_Env_correlation.pdf", fig.cor.env,
       width = 5, height = 3, scale = 2)
