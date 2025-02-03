# script for enrichment

########################################################################################### requirements ####

require(DESeq2)
require(edgeR)
require(topGO)
library(data.table)
library(dplyr)
library(formattable)
library(tidyr)
library(htmltools)
library(webshot)    
library(scales)
library(ggplot2)
library(rrvgo)
library(org.Dm.eg.db)

export_formattable <- function(f, file, width = "100%", height = NULL, 
                               background = "white", delay = 0.2)
{
  w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(w, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay)
}

########################################################################################### parameters ####

pvgo <- 0.01        # pvalue for enrichment
pvde <- 0.05        # pvalue for DE
lgfc <- 1.5         # logFold change
ndsz <- 05          # node size for enrichment
algr <- "classic"   # algorythm for enrichment
onto <- "BP"        # GO ontology
threshold <- 1    # clustering threshold
threshold_all <- 0.7    # clustering threshold

########################################################################################### DE analyes AD_AvsB ####

rnaseqMatrix <- read.table("RSEM_crema.filtered.gene.counts.matrix", header=T)

rnaseqMatrix_red <- round(rnaseqMatrix[c(1:5,11:15)])
rnaseqMatrix_red <- rnaseqMatrix_red[rowSums(cpm(rnaseqMatrix_red) > 1) >= 2,]

conditions = data.frame(conditions=factor(c(rep("controllo", 5), rep("trattato", 5))))

rownames(conditions) = colnames(rnaseqMatrix_red)
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = rnaseqMatrix_red,
                                            colData = conditions,
                                            design = ~ conditions)
dds = DESeq(ddsFullCountTable)
contrast=c("conditions","trattato","controllo")
res = results(dds, contrast)
baseMean_controllo <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "controllo"])
baseMean_trattato <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "trattato"])
res = cbind(baseMean_controllo, baseMean_trattato, as.data.frame(res))
res = cbind(sampleA="A", sampleB="B", as.data.frame(res))
res$padj[is.na(res$padj)]  <- 1
res = as.data.frame(res[order(res$pvalue),])

up <- subset(res, padj < pvde & log2FoldChange > +lgfc)
up <- up[order(up$log2FoldChange,decreasing = T),]
write(row.names(up), "AD/AD_AvsB_UP_genes_lst.txt")

dn <- subset(res, padj < pvde & log2FoldChange < -lgfc)
dn <- dn[order(dn$log2FoldChange,decreasing = T),]
write(row.names(dn), "AD/AD_AvsB_DN_genes_lst.txt")

res[row.names(res) %in% c("TRINITY_DN1729_c0_g1",
                          "TRINITY_DN1268_c0_g1",
                          "TRINITY_DN2392_c0_g1",
                          "TRINITY_DN3637_c0_g2",
                          "TRINITY_DN4809_c0_g1"),]

########################################################################################### go enrichment AD_AvsB_UP ####

geneID2GO <- readMappings(file = "GO_crema_geneUniverse_AD")

geneUniverse <- names(geneID2GO)

genesOfInterest <- read.table("AD/AD_AvsB_UP_genes_lst.txt")

genesOfInterest <- as.character(genesOfInterest$V1)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse
myGOdata <- new("topGOdata", description="", ontology=onto, allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize=ndsz)
test <- runTest(myGOdata, algorithm=algr, statistic="fisher")
table <- GenTable(myGOdata,
                  test = test,
                  orderBy = "test", topNodes=1000)
table <- subset(table, test < pvgo)
table <- table[order(table$test,decreasing = F),]

write.table(table,"AD_AvsB_UP.tsv",  row.names=F, quote=F, sep = "\t")

row.names(table) <- NULL

table$color <- -log(as.numeric(table$test))

table$ratio <- table$Significant / table$Expected
table$ratio <- round(table$ratio,2)

colnames(table) <- c("id","term","Annotated","Significant","Expected","p","color","sig/exp")
table[nrow(table) + 1,] = c(" "," "," "," "," ",pvgo,-log(0.01),0)
table[nrow(table) + 1,] = c(" "," "," "," "," ",0,Inf,0)

plot <- formattable(table, list(color = FALSE, Expected = FALSE,
                                `id` = formatter("span", style = ~ style(color = "black",font.weight = "bold")), 
                                `p`= color_tile("orange", "white"), 
                                `sig/exp`= color_bar("orange")),
                    table.attr = 'style="font-size: 20px; font-family: Arial";\"')

export_formattable(plot,"AD_AvsB_UP.png")


#

simMatrix <- calculateSimMatrix(head(table$id, -2),
                                orgdb="org.Dm.eg.db",
                                ont="BP",
                                method="Rel")

scores <- setNames(-log10(as.numeric(head(table$p, -2))), head(table$id, -2))
reducedTerms_AD_AvsB_UP <- reduceSimMatrix(simMatrix,
                                           scores,
                                           threshold=threshold,
                                           orgdb="org.Dm.eg.db")

names(reducedTerms_AD_AvsB_UP)[names(reducedTerms_AD_AvsB_UP) == 'go'] <- 'id'
table$term <- NULL
reducedTerms_AD_AvsB_UP <- merge(reducedTerms_AD_AvsB_UP, table, by ="id")

########################################################################################### go enrichment AD_AvsB_DN ####

geneID2GO <- readMappings(file = "GO_crema_geneUniverse_AD")

geneUniverse <- names(geneID2GO)

genesOfInterest <- read.table("AD/AD_AvsB_DN_genes_lst.txt")

genesOfInterest <- as.character(genesOfInterest$V1)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse
myGOdata <- new("topGOdata", description="", ontology=onto, allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize=ndsz)
test <- runTest(myGOdata, algorithm=algr, statistic="fisher")
table <- GenTable(myGOdata,
                  test = test,
                  orderBy = "test", topNodes=1000)
table <- subset(table, test < pvgo)
table <- table[order(table$test,decreasing = F),]

write.table(table,"AD_AvsB_DN.tsv",  row.names=F, quote=F, sep = "\t")

row.names(table) <- NULL

table$color <- -log(as.numeric(table$test))

table$ratio <- table$Significant / table$Expected
table$ratio <- round(table$ratio,2)

colnames(table) <- c("id","term","Annotated","Significant","Expected","p","color","sig/exp")
table[nrow(table) + 1,] = c(" "," "," "," "," ",pvgo,-log(0.01),0)
table[nrow(table) + 1,] = c(" "," "," "," "," ",0,Inf,0)

plot <- formattable(table, list(color = FALSE, Expected = FALSE,
                                `id` = formatter("span", style = ~ style(color = "black",font.weight = "bold")), 
                                `p`= color_tile("lightblue", "white"), 
                                `sig/exp`= color_bar("lightblue")),
                    table.attr = 'style="font-size: 20px; font-family: Arial";\"')

export_formattable(plot,"AD_AvsB_DN.png")

#

simMatrix <- calculateSimMatrix(head(table$id, -2),
                                orgdb="org.Dm.eg.db",
                                ont="BP",
                                method="Rel")

scores <- setNames(-log10(as.numeric(head(table$p, -2))), head(table$id, -2))
reducedTerms_AD_AvsB_DN <- reduceSimMatrix(simMatrix,
                                           scores,
                                           threshold=threshold,
                                           orgdb="org.Dm.eg.db")

names(reducedTerms_AD_AvsB_DN)[names(reducedTerms_AD_AvsB_DN) == 'go'] <- 'id'
table$term <- NULL
reducedTerms_AD_AvsB_DN <- merge(reducedTerms_AD_AvsB_DN, table, by ="id")

########################################################################################### DE analyes AD_AvsC ####

rnaseqMatrix <- read.table("RSEM_crema.filtered.gene.counts.matrix", header=T)

rnaseqMatrix_red <- round(rnaseqMatrix[c(1:5,21:25)])
rnaseqMatrix_red <- rnaseqMatrix_red[rowSums(cpm(rnaseqMatrix_red) > 1) >= 2,]

conditions = data.frame(conditions=factor(c(rep("controllo", 5), rep("trattato", 5))))

rownames(conditions) = colnames(rnaseqMatrix_red)
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = rnaseqMatrix_red,
                                            colData = conditions,
                                            design = ~ conditions)
dds = DESeq(ddsFullCountTable)
contrast=c("conditions","trattato","controllo")
res = results(dds, contrast)
baseMean_controllo <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "controllo"])
baseMean_trattato <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "trattato"])
res = cbind(baseMean_controllo, baseMean_trattato, as.data.frame(res))
res = cbind(sampleA="A", sampleB="B", as.data.frame(res))
res$padj[is.na(res$padj)]  <- 1
res = as.data.frame(res[order(res$pvalue),])

up <- subset(res, padj < pvde & log2FoldChange > +lgfc)
up <- up[order(up$log2FoldChange,decreasing = T),]
write(row.names(up), "AD/AD_AvsC_UP_genes_lst.txt")

dn <- subset(res, padj < pvde & log2FoldChange < -lgfc)
dn <- dn[order(dn$log2FoldChange,decreasing = T),]
write(row.names(dn), "AD/AD_AvsC_DN_genes_lst.txt")

res[row.names(res) %in% c("TRINITY_DN1729_c0_g1",
                          "TRINITY_DN1268_c0_g1",
                          "TRINITY_DN2392_c0_g1",
                          "TRINITY_DN3637_c0_g2",
                          "TRINITY_DN4809_c0_g1"),]

########################################################################################### go enrichment AD_AvsC_UP ####

geneID2GO <- readMappings(file = "GO_crema_geneUniverse_AD")

geneUniverse <- names(geneID2GO)

genesOfInterest <- read.table("AD/AD_AvsC_UP_genes_lst.txt")

genesOfInterest <- as.character(genesOfInterest$V1)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse
myGOdata <- new("topGOdata", description="", ontology=onto, allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize=ndsz)
test <- runTest(myGOdata, algorithm=algr, statistic="fisher")
table <- GenTable(myGOdata,
                  test = test,
                  orderBy = "test", topNodes=1000)
table <- subset(table, test < pvgo)
table <- table[order(table$test,decreasing = F),]

write.table(table,"AD_AvsC_UP.tsv",  row.names=F, quote=F, sep = "\t")

row.names(table) <- NULL

table$color <- -log(as.numeric(table$test))

table$ratio <- table$Significant / table$Expected
table$ratio <- round(table$ratio,2)

colnames(table) <- c("id","term","Annotated","Significant","Expected","p","color","sig/exp")
table[nrow(table) + 1,] = c(" "," "," "," "," ",pvgo,-log(0.01),0)
table[nrow(table) + 1,] = c(" "," "," "," "," ",0,Inf,0)

plot <- formattable(table, list(color = FALSE, Expected = FALSE,
                                `id` = formatter("span", style = ~ style(color = "black",font.weight = "bold")), 
                                `p`= color_tile("orange", "white"), 
                                `sig/exp`= color_bar("orange")),
                    table.attr = 'style="font-size: 20px; font-family: Arial";\"')

export_formattable(plot,"AD_AvsC_UP.png")


#

simMatrix <- calculateSimMatrix(head(table$id, -2),
                                orgdb="org.Dm.eg.db",
                                ont="BP",
                                method="Rel")

scores <- setNames(-log10(as.numeric(head(table$p, -2))), head(table$id, -2))
reducedTerms_AD_AvsC_UP <- reduceSimMatrix(simMatrix,
                                           scores,
                                           threshold=threshold,
                                           orgdb="org.Dm.eg.db")

names(reducedTerms_AD_AvsC_UP)[names(reducedTerms_AD_AvsC_UP) == 'go'] <- 'id'
table$term <- NULL
reducedTerms_AD_AvsC_UP <- merge(reducedTerms_AD_AvsC_UP, table, by ="id")

########################################################################################### go enrichment AD_AvsC_DN ####

geneID2GO <- readMappings(file = "GO_crema_geneUniverse_AD")

geneUniverse <- names(geneID2GO)

genesOfInterest <- read.table("AD/AD_AvsC_DN_genes_lst.txt")

genesOfInterest <- as.character(genesOfInterest$V1)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse
myGOdata <- new("topGOdata", description="", ontology=onto, allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize=ndsz)
test <- runTest(myGOdata, algorithm=algr, statistic="fisher")
table <- GenTable(myGOdata,
                  test = test,
                  orderBy = "test", topNodes=1000)
table <- subset(table, test < pvgo)
table <- table[order(table$test,decreasing = F),]

write.table(table,"AD_AvsC_DN.tsv",  row.names=F, quote=F, sep = "\t")

row.names(table) <- NULL

table$color <- -log(as.numeric(table$test))

table$ratio <- table$Significant / table$Expected
table$ratio <- round(table$ratio,2)

colnames(table) <- c("id","term","Annotated","Significant","Expected","p","color","sig/exp")
table[nrow(table) + 1,] = c(" "," "," "," "," ",pvgo,-log(0.01),0)
table[nrow(table) + 1,] = c(" "," "," "," "," ",0,Inf,0)

plot <- formattable(table, list(color = FALSE, Expected = FALSE,
                                `id` = formatter("span", style = ~ style(color = "black",font.weight = "bold")), 
                                `p`= color_tile("lightblue", "white"), 
                                `sig/exp`= color_bar("lightblue")),
                    table.attr = 'style="font-size: 20px; font-family: Arial";\"')

export_formattable(plot,"AD_AvsC_DN.png")

#

simMatrix <- calculateSimMatrix(head(table$id, -2),
                                orgdb="org.Dm.eg.db",
                                ont="BP",
                                method="Rel")

scores <- setNames(-log10(as.numeric(head(table$p, -2))), head(table$id, -2))
reducedTerms_AD_AvsC_DN <- reduceSimMatrix(simMatrix,
                                           scores,
                                           threshold=threshold,
                                           orgdb="org.Dm.eg.db")

names(reducedTerms_AD_AvsC_DN)[names(reducedTerms_AD_AvsC_DN) == 'go'] <- 'id'
table$term <- NULL
reducedTerms_AD_AvsC_DN <- merge(reducedTerms_AD_AvsC_DN, table, by ="id")

########################################################################################### DE analyes AD_AvsD ####

rnaseqMatrix <- read.table("RSEM_crema.filtered.gene.counts.matrix", header=T)

rnaseqMatrix_red <- round(rnaseqMatrix[c(1:5,31:35)])
rnaseqMatrix_red <- rnaseqMatrix_red[rowSums(cpm(rnaseqMatrix_red) > 1) >= 2,]

conditions = data.frame(conditions=factor(c(rep("controllo", 5), rep("trattato", 5))))

rownames(conditions) = colnames(rnaseqMatrix_red)
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = rnaseqMatrix_red,
                                            colData = conditions,
                                            design = ~ conditions)
dds = DESeq(ddsFullCountTable)
contrast=c("conditions","trattato","controllo")
res = results(dds, contrast)
baseMean_controllo <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "controllo"])
baseMean_trattato <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "trattato"])
res = cbind(baseMean_controllo, baseMean_trattato, as.data.frame(res))
res = cbind(sampleA="A", sampleB="B", as.data.frame(res))
res$padj[is.na(res$padj)]  <- 1
res = as.data.frame(res[order(res$pvalue),])

up <- subset(res, padj < pvde & log2FoldChange > +lgfc)
up <- up[order(up$log2FoldChange,decreasing = T),]
write(row.names(up), "AD/AD_AvsD_UP_genes_lst.txt")

dn <- subset(res, padj < pvde & log2FoldChange < -lgfc)
dn <- dn[order(dn$log2FoldChange,decreasing = T),]
write(row.names(dn), "AD/AD_AvsD_DN_genes_lst.txt")

res[row.names(res) %in% c("TRINITY_DN1729_c0_g1",
                          "TRINITY_DN1268_c0_g1",
                          "TRINITY_DN2392_c0_g1",
                          "TRINITY_DN3637_c0_g2",
                          "TRINITY_DN4809_c0_g1"),]

########################################################################################### go enrichment AD_AvsD_UP ####

geneID2GO <- readMappings(file = "GO_crema_geneUniverse_AD")

geneUniverse <- names(geneID2GO)

genesOfInterest <- read.table("AD/AD_AvsD_UP_genes_lst.txt")

genesOfInterest <- as.character(genesOfInterest$V1)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse
myGOdata <- new("topGOdata", description="", ontology=onto, allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize=ndsz)
test <- runTest(myGOdata, algorithm=algr, statistic="fisher")
table <- GenTable(myGOdata,
                  test = test,
                  orderBy = "test", topNodes=1000)
table <- subset(table, test < pvgo)
table <- table[order(table$test,decreasing = F),]

write.table(table,"AD_AvsD_UP.tsv",  row.names=F, quote=F, sep = "\t")

row.names(table) <- NULL

table$color <- -log(as.numeric(table$test))

table$ratio <- table$Significant / table$Expected
table$ratio <- round(table$ratio,2)

colnames(table) <- c("id","term","Annotated","Significant","Expected","p","color","sig/exp")
table[nrow(table) + 1,] = c(" "," "," "," "," ",pvgo,-log(0.01),0)
table[nrow(table) + 1,] = c(" "," "," "," "," ",0,Inf,0)

plot <- formattable(table, list(color = FALSE, Expected = FALSE,
                                `id` = formatter("span", style = ~ style(color = "black",font.weight = "bold")), 
                                `p`= color_tile("orange", "white"), 
                                `sig/exp`= color_bar("orange")),
                    table.attr = 'style="font-size: 20px; font-family: Arial";\"')

export_formattable(plot,"AD_AvsD_UP.png")


#

simMatrix <- calculateSimMatrix(head(table$id, -2),
                                orgdb="org.Dm.eg.db",
                                ont="BP",
                                method="Rel")

scores <- setNames(-log10(as.numeric(head(table$p, -2))), head(table$id, -2))
reducedTerms_AD_AvsD_UP <- reduceSimMatrix(simMatrix,
                                           scores,
                                           threshold=threshold,
                                           orgdb="org.Dm.eg.db")

names(reducedTerms_AD_AvsD_UP)[names(reducedTerms_AD_AvsD_UP) == 'go'] <- 'id'
table$term <- NULL
reducedTerms_AD_AvsD_UP <- merge(reducedTerms_AD_AvsD_UP, table, by ="id")

########################################################################################### go enrichment AD_AvsD_DN ####

geneID2GO <- readMappings(file = "GO_crema_geneUniverse_AD")

geneUniverse <- names(geneID2GO)

genesOfInterest <- read.table("AD/AD_AvsD_DN_genes_lst.txt")

genesOfInterest <- as.character(genesOfInterest$V1)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse
myGOdata <- new("topGOdata", description="", ontology=onto, allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize=ndsz)
test <- runTest(myGOdata, algorithm=algr, statistic="fisher")
table <- GenTable(myGOdata,
                  test = test,
                  orderBy = "test", topNodes=1000)
table <- subset(table, test < pvgo)
table <- table[order(table$test,decreasing = F),]

write.table(table,"AD_AvsD_DN.tsv",  row.names=F, quote=F, sep = "\t")

row.names(table) <- NULL

table$color <- -log(as.numeric(table$test))

table$ratio <- table$Significant / table$Expected
table$ratio <- round(table$ratio,2)

colnames(table) <- c("id","term","Annotated","Significant","Expected","p","color","sig/exp")
table[nrow(table) + 1,] = c(" "," "," "," "," ",pvgo,-log(0.01),0)
table[nrow(table) + 1,] = c(" "," "," "," "," ",0,Inf,0)

plot <- formattable(table, list(color = FALSE, Expected = FALSE,
                                `id` = formatter("span", style = ~ style(color = "black",font.weight = "bold")), 
                                `p`= color_tile("lightblue", "white"), 
                                `sig/exp`= color_bar("lightblue")),
                    table.attr = 'style="font-size: 20px; font-family: Arial";\"')

export_formattable(plot,"AD_AvsD_DN.png")

#

simMatrix <- calculateSimMatrix(head(table$id, -2),
                                orgdb="org.Dm.eg.db",
                                ont="BP",
                                method="Rel")

scores <- setNames(-log10(as.numeric(head(table$p, -2))), head(table$id, -2))
reducedTerms_AD_AvsD_DN <- reduceSimMatrix(simMatrix,
                                           scores,
                                           threshold=threshold,
                                           orgdb="org.Dm.eg.db")

names(reducedTerms_AD_AvsD_DN)[names(reducedTerms_AD_AvsD_DN) == 'go'] <- 'id'
table$term <- NULL
reducedTerms_AD_AvsD_DN <- merge(reducedTerms_AD_AvsD_DN, table, by ="id")

########################################################################################### DE analyes CT_AvsB ####

rnaseqMatrix <- read.table("RSEM_crema.filtered.gene.counts.matrix", header=T)

rnaseqMatrix_red <- round(rnaseqMatrix[c(6:10,16:20)])
rnaseqMatrix_red <- rnaseqMatrix_red[rowSums(cpm(rnaseqMatrix_red) > 1) >= 2,]

conditions = data.frame(conditions=factor(c(rep("controllo", 5), rep("trattato", 5))))

rownames(conditions) = colnames(rnaseqMatrix_red)
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = rnaseqMatrix_red,
                                            colData = conditions,
                                            design = ~ conditions)
dds = DESeq(ddsFullCountTable)
contrast=c("conditions","trattato","controllo")
res = results(dds, contrast)
baseMean_controllo <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "controllo"])
baseMean_trattato <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "trattato"])
res = cbind(baseMean_controllo, baseMean_trattato, as.data.frame(res))
res = cbind(sampleA="A", sampleB="B", as.data.frame(res))
res$padj[is.na(res$padj)]  <- 1
res = as.data.frame(res[order(res$pvalue),])

up <- subset(res, padj < pvde & log2FoldChange > +lgfc)
up <- up[order(up$log2FoldChange,decreasing = T),]
write(row.names(up), "CT/CT_AvsB_UP_genes_lst.txt")

dn <- subset(res, padj < pvde & log2FoldChange < -lgfc)
dn <- dn[order(dn$log2FoldChange,decreasing = T),]
write(row.names(dn), "CT/CT_AvsB_DN_genes_lst.txt")

res[row.names(res) %in% c("TRINITY_DN1729_c0_g1",
                          "TRINITY_DN1268_c0_g1",
                          "TRINITY_DN2392_c0_g1",
                          "TRINITY_DN3637_c0_g2",
                          "TRINITY_DN4809_c0_g1"),]

########################################################################################### go enrichment CT_AvsB_UP ####

geneID2GO <- readMappings(file = "GO_crema_geneUniverse_CT")

geneUniverse <- names(geneID2GO)

genesOfInterest <- read.table("CT/CT_AvsB_UP_genes_lst.txt")

genesOfInterest <- as.character(genesOfInterest$V1)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse
myGOdata <- new("topGOdata", description="", ontology=onto, allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize=ndsz)
test <- runTest(myGOdata, algorithm=algr, statistic="fisher")
table <- GenTable(myGOdata,
                  test = test,
                  orderBy = "test", topNodes=1000)
table <- subset(table, test < pvgo)
table <- table[order(table$test,decreasing = F),]

write.table(table,"CT_AvsB_UP.tsv",  row.names=F, quote=F, sep = "\t")

row.names(table) <- NULL

table$color <- -log(as.numeric(table$test))

table$ratio <- table$Significant / table$Expected
table$ratio <- round(table$ratio,2)

colnames(table) <- c("id","term","Annotated","Significant","Expected","p","color","sig/exp")
table[nrow(table) + 1,] = c(" "," "," "," "," ",pvgo,-log(0.01),0)
table[nrow(table) + 1,] = c(" "," "," "," "," ",0,Inf,0)

plot <- formattable(table, list(color = FALSE, Expected = FALSE,
                                `id` = formatter("span", style = ~ style(color = "black",font.weight = "bold")), 
                                `p`= color_tile("orange", "white"), 
                                `sig/exp`= color_bar("orange")),
                    table.attr = 'style="font-size: 20px; font-family: Arial";\"')

export_formattable(plot,"CT_AvsB_UP.png")


#

simMatrix <- calculateSimMatrix(head(table$id, -2),
                                orgdb="org.Dm.eg.db",
                                ont="BP",
                                method="Rel")

scores <- setNames(-log10(as.numeric(head(table$p, -2))), head(table$id, -2))
reducedTerms_CT_AvsB_UP <- reduceSimMatrix(simMatrix,
                                           scores,
                                           threshold=threshold,
                                           orgdb="org.Dm.eg.db")

names(reducedTerms_CT_AvsB_UP)[names(reducedTerms_CT_AvsB_UP) == 'go'] <- 'id'
table$term <- NULL
reducedTerms_CT_AvsB_UP <- merge(reducedTerms_CT_AvsB_UP, table, by ="id")

########################################################################################### go enrichment CT_AvsB_DN ####

geneID2GO <- readMappings(file = "GO_crema_geneUniverse_CT")

geneUniverse <- names(geneID2GO)

genesOfInterest <- read.table("CT/CT_AvsB_DN_genes_lst.txt")

genesOfInterest <- as.character(genesOfInterest$V1)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse
myGOdata <- new("topGOdata", description="", ontology=onto, allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize=ndsz)
test <- runTest(myGOdata, algorithm=algr, statistic="fisher")
table <- GenTable(myGOdata,
                  test = test,
                  orderBy = "test", topNodes=1000)
table <- subset(table, test < pvgo)
table <- table[order(table$test,decreasing = F),]

write.table(table,"CT_AvsB_DN.tsv",  row.names=F, quote=F, sep = "\t")

row.names(table) <- NULL

table$color <- -log(as.numeric(table$test))

table$ratio <- table$Significant / table$Expected
table$ratio <- round(table$ratio,2)

colnames(table) <- c("id","term","Annotated","Significant","Expected","p","color","sig/exp")
table[nrow(table) + 1,] = c(" "," "," "," "," ",pvgo,-log(0.01),0)
table[nrow(table) + 1,] = c(" "," "," "," "," ",0,Inf,0)

plot <- formattable(table, list(color = FALSE, Expected = FALSE,
                                `id` = formatter("span", style = ~ style(color = "black",font.weight = "bold")), 
                                `p`= color_tile("lightblue", "white"), 
                                `sig/exp`= color_bar("lightblue")),
                    table.attr = 'style="font-size: 20px; font-family: Arial";\"')

export_formattable(plot,"CT_AvsB_DN.png")

#

simMatrix <- calculateSimMatrix(head(table$id, -2),
                                orgdb="org.Dm.eg.db",
                                ont="BP",
                                method="Rel")

scores <- setNames(-log10(as.numeric(head(table$p, -2))), head(table$id, -2))
reducedTerms_CT_AvsB_DN <- reduceSimMatrix(simMatrix,
                                           scores,
                                           threshold=threshold,
                                           orgdb="org.Dm.eg.db")

names(reducedTerms_CT_AvsB_DN)[names(reducedTerms_CT_AvsB_DN) == 'go'] <- 'id'
table$term <- NULL
reducedTerms_CT_AvsB_DN <- merge(reducedTerms_CT_AvsB_DN, table, by ="id")

########################################################################################### DE analyes CT_AvsC ####

rnaseqMatrix <- read.table("RSEM_crema.filtered.gene.counts.matrix", header=T)

rnaseqMatrix_red <- round(rnaseqMatrix[c(6:10,26:30)])
rnaseqMatrix_red <- rnaseqMatrix_red[rowSums(cpm(rnaseqMatrix_red) > 1) >= 2,]

conditions = data.frame(conditions=factor(c(rep("controllo", 5), rep("trattato", 5))))

rownames(conditions) = colnames(rnaseqMatrix_red)
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = rnaseqMatrix_red,
                                            colData = conditions,
                                            design = ~ conditions)
dds = DESeq(ddsFullCountTable)
contrast=c("conditions","trattato","controllo")
res = results(dds, contrast)
baseMean_controllo <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "controllo"])
baseMean_trattato <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "trattato"])
res = cbind(baseMean_controllo, baseMean_trattato, as.data.frame(res))
res = cbind(sampleA="A", sampleB="B", as.data.frame(res))
res$padj[is.na(res$padj)]  <- 1
res = as.data.frame(res[order(res$pvalue),])

up <- subset(res, padj < pvde & log2FoldChange > +lgfc)
up <- up[order(up$log2FoldChange,decreasing = T),]
write(row.names(up), "CT/CT_AvsC_UP_genes_lst.txt")

dn <- subset(res, padj < pvde & log2FoldChange < -lgfc)
dn <- dn[order(dn$log2FoldChange,decreasing = T),]
write(row.names(dn), "CT/CT_AvsC_DN_genes_lst.txt")

res[row.names(res) %in% c("TRINITY_DN1729_c0_g1",
                          "TRINITY_DN1268_c0_g1",
                          "TRINITY_DN2392_c0_g1",
                          "TRINITY_DN3637_c0_g2",
                          "TRINITY_DN4809_c0_g1"),]

########################################################################################### go enrichment CT_AvsC_UP ####

geneID2GO <- readMappings(file = "GO_crema_geneUniverse_CT")

geneUniverse <- names(geneID2GO)

genesOfInterest <- read.table("CT/CT_AvsC_UP_genes_lst.txt")

genesOfInterest <- as.character(genesOfInterest$V1)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse
myGOdata <- new("topGOdata", description="", ontology=onto, allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize=ndsz)
test <- runTest(myGOdata, algorithm=algr, statistic="fisher")
table <- GenTable(myGOdata,
                  test = test,
                  orderBy = "test", topNodes=1000)
table <- subset(table, test < pvgo)
table <- table[order(table$test,decreasing = F),]

write.table(table,"CT_AvsC_UP.tsv",  row.names=F, quote=F, sep = "\t")

row.names(table) <- NULL

table$color <- -log(as.numeric(table$test))

table$ratio <- table$Significant / table$Expected
table$ratio <- round(table$ratio,2)

colnames(table) <- c("id","term","Annotated","Significant","Expected","p","color","sig/exp")
table[nrow(table) + 1,] = c(" "," "," "," "," ",pvgo,-log(0.01),0)
table[nrow(table) + 1,] = c(" "," "," "," "," ",0,Inf,0)

plot <- formattable(table, list(color = FALSE, Expected = FALSE,
                                `id` = formatter("span", style = ~ style(color = "black",font.weight = "bold")), 
                                `p`= color_tile("orange", "white"), 
                                `sig/exp`= color_bar("orange")),
                    table.attr = 'style="font-size: 20px; font-family: Arial";\"')

export_formattable(plot,"CT_AvsC_UP.png")


#

simMatrix <- calculateSimMatrix(head(table$id, -2),
                                orgdb="org.Dm.eg.db",
                                ont="BP",
                                method="Rel")

scores <- setNames(-log10(as.numeric(head(table$p, -2))), head(table$id, -2))
reducedTerms_CT_AvsC_UP <- reduceSimMatrix(simMatrix,
                                           scores,
                                           threshold=threshold,
                                           orgdb="org.Dm.eg.db")

names(reducedTerms_CT_AvsC_UP)[names(reducedTerms_CT_AvsC_UP) == 'go'] <- 'id'
table$term <- NULL
reducedTerms_CT_AvsC_UP <- merge(reducedTerms_CT_AvsC_UP, table, by ="id")
########################################################################################### go enrichment CT_AvsC_DN ####

geneID2GO <- readMappings(file = "GO_crema_geneUniverse_CT")

geneUniverse <- names(geneID2GO)

genesOfInterest <- read.table("CT/CT_AvsC_DN_genes_lst.txt")

genesOfInterest <- as.character(genesOfInterest$V1)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse
myGOdata <- new("topGOdata", description="", ontology=onto, allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize=ndsz)
test <- runTest(myGOdata, algorithm=algr, statistic="fisher")
table <- GenTable(myGOdata,
                  test = test,
                  orderBy = "test", topNodes=1000)
table <- subset(table, test < pvgo)
table <- table[order(table$test,decreasing = F),]

write.table(table,"CT_AvsC_DN.tsv",  row.names=F, quote=F, sep = "\t")

row.names(table) <- NULL

table$color <- -log(as.numeric(table$test))

table$ratio <- table$Significant / table$Expected
table$ratio <- round(table$ratio,2)

colnames(table) <- c("id","term","Annotated","Significant","Expected","p","color","sig/exp")
table[nrow(table) + 1,] = c(" "," "," "," "," ",pvgo,-log(0.01),0)
table[nrow(table) + 1,] = c(" "," "," "," "," ",0,Inf,0)

plot <- formattable(table, list(color = FALSE, Expected = FALSE,
                                `id` = formatter("span", style = ~ style(color = "black",font.weight = "bold")), 
                                `p`= color_tile("lightblue", "white"), 
                                `sig/exp`= color_bar("lightblue")),
                    table.attr = 'style="font-size: 20px; font-family: Arial";\"')

export_formattable(plot,"CT_AvsC_DN.png")

#

simMatrix <- calculateSimMatrix(head(table$id, -2),
                                orgdb="org.Dm.eg.db",
                                ont="BP",
                                method="Rel")

scores <- setNames(-log10(as.numeric(head(table$p, -2))), head(table$id, -2))
reducedTerms_CT_AvsC_DN <- reduceSimMatrix(simMatrix,
                                           scores,
                                           threshold=threshold,
                                           orgdb="org.Dm.eg.db")

names(reducedTerms_CT_AvsC_DN)[names(reducedTerms_CT_AvsC_DN) == 'go'] <- 'id'
table$term <- NULL
reducedTerms_CT_AvsC_DN <- merge(reducedTerms_CT_AvsC_DN, table, by ="id")

########################################################################################### DE analyes CT_AvsD ####

rnaseqMatrix <- read.table("RSEM_crema.filtered.gene.counts.matrix", header=T)

rnaseqMatrix_red <- round(rnaseqMatrix[c(6:10,36:40)])
rnaseqMatrix_red <- rnaseqMatrix_red[rowSums(cpm(rnaseqMatrix_red) > 1) >= 2,]

conditions = data.frame(conditions=factor(c(rep("controllo", 5), rep("trattato", 5))))

rownames(conditions) = colnames(rnaseqMatrix_red)
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = rnaseqMatrix_red,
                                            colData = conditions,
                                            design = ~ conditions)
dds = DESeq(ddsFullCountTable)
contrast=c("conditions","trattato","controllo")
res = results(dds, contrast)
baseMean_controllo <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "controllo"])
baseMean_trattato <- rowMeans(counts(dds, normalized=TRUE)[,colData(dds)$conditions == "trattato"])
res = cbind(baseMean_controllo, baseMean_trattato, as.data.frame(res))
res = cbind(sampleA="A", sampleB="B", as.data.frame(res))
res$padj[is.na(res$padj)]  <- 1
res = as.data.frame(res[order(res$pvalue),])

up <- subset(res, padj < pvde & log2FoldChange > +lgfc)
up <- up[order(up$log2FoldChange,decreasing = T),]
write(row.names(up), "CT/CT_AvsD_UP_genes_lst.txt")

dn <- subset(res, padj < pvde & log2FoldChange < -lgfc)
dn <- dn[order(dn$log2FoldChange,decreasing = T),]
write(row.names(dn), "CT/CT_AvsD_DN_genes_lst.txt")

res[row.names(res) %in% c("TRINITY_DN1729_c0_g1",
                          "TRINITY_DN1268_c0_g1",
                          "TRINITY_DN2392_c0_g1",
                          "TRINITY_DN3637_c0_g2",
                          "TRINITY_DN4809_c0_g1"),]

########################################################################################### go enrichment CT_AvsD_UP ####

geneID2GO <- readMappings(file = "GO_crema_geneUniverse_CT")

geneUniverse <- names(geneID2GO)

genesOfInterest <- read.table("CT/CT_AvsD_UP_genes_lst.txt")

genesOfInterest <- as.character(genesOfInterest$V1)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse
myGOdata <- new("topGOdata", description="", ontology=onto, allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize=ndsz)
test <- runTest(myGOdata, algorithm=algr, statistic="fisher")
table <- GenTable(myGOdata,
                  test = test,
                  orderBy = "test", topNodes=1000)
table <- subset(table, test < pvgo)
table <- table[order(table$test,decreasing = F),]

write.table(table,"CT_AvsD_UP.tsv",  row.names=F, quote=F, sep = "\t")

row.names(table) <- NULL

table$color <- -log(as.numeric(table$test))

table$ratio <- table$Significant / table$Expected
table$ratio <- round(table$ratio,2)

colnames(table) <- c("id","term","Annotated","Significant","Expected","p","color","sig/exp")
table[nrow(table) + 1,] = c(" "," "," "," "," ",pvgo,-log(0.01),0)
table[nrow(table) + 1,] = c(" "," "," "," "," ",0,Inf,0)

plot <- formattable(table, list(color = FALSE, Expected = FALSE,
                                `id` = formatter("span", style = ~ style(color = "black",font.weight = "bold")), 
                                `p`= color_tile("orange", "white"), 
                                `sig/exp`= color_bar("orange")),
                    table.attr = 'style="font-size: 20px; font-family: Arial";\"')

export_formattable(plot,"CT_AvsD_UP.png")


#

simMatrix <- calculateSimMatrix(head(table$id, -2),
                                orgdb="org.Dm.eg.db",
                                ont="BP",
                                method="Rel")

scores <- setNames(-log10(as.numeric(head(table$p, -2))), head(table$id, -2))
reducedTerms_CT_AvsD_UP <- reduceSimMatrix(simMatrix,
                                           scores,
                                           threshold=threshold,
                                           orgdb="org.Dm.eg.db")

names(reducedTerms_CT_AvsD_UP)[names(reducedTerms_CT_AvsD_UP) == 'go'] <- 'id'
table$term <- NULL
reducedTerms_CT_AvsD_UP <- merge(reducedTerms_CT_AvsD_UP, table, by ="id")

########################################################################################### go enrichment CT_AvsD_DN ####

geneID2GO <- readMappings(file = "GO_crema_geneUniverse_CT")

geneUniverse <- names(geneID2GO)

genesOfInterest <- read.table("CT/CT_AvsD_DN_genes_lst.txt")

genesOfInterest <- as.character(genesOfInterest$V1)
geneList <- factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) <- geneUniverse
myGOdata <- new("topGOdata", description="", ontology=onto, allGenes=geneList,  annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize=ndsz)
test <- runTest(myGOdata, algorithm=algr, statistic="fisher")
table <- GenTable(myGOdata,
                  test = test,
                  orderBy = "test", topNodes=1000)
table <- subset(table, test < pvgo)
table <- table[order(table$test,decreasing = F),]

write.table(table,"CT_AvsD_DN.tsv",  row.names=F, quote=F, sep = "\t")

row.names(table) <- NULL

table$color <- -log(as.numeric(table$test))

table$ratio <- table$Significant / table$Expected
table$ratio <- round(table$ratio,2)

colnames(table) <- c("id","term","Annotated","Significant","Expected","p","color","sig/exp")
table[nrow(table) + 1,] = c(" "," "," "," "," ",pvgo,-log(0.01),0)
table[nrow(table) + 1,] = c(" "," "," "," "," ",0,Inf,0)

plot <- formattable(table, list(color = FALSE, Expected = FALSE,
                                `id` = formatter("span", style = ~ style(color = "black",font.weight = "bold")), 
                                `p`= color_tile("lightblue", "white"), 
                                `sig/exp`= color_bar("lightblue")),
                    table.attr = 'style="font-size: 20px; font-family: Arial";\"')

export_formattable(plot,"CT_AvsD_DN.png")

#

simMatrix <- calculateSimMatrix(head(table$id, -2),
                                orgdb="org.Dm.eg.db",
                                ont="BP",
                                method="Rel")

scores <- setNames(-log10(as.numeric(head(table$p, -2))), head(table$id, -2))
reducedTerms_CT_AvsD_DN <- reduceSimMatrix(simMatrix,
                                           scores,
                                           threshold=threshold,
                                           orgdb="org.Dm.eg.db")

names(reducedTerms_CT_AvsD_DN)[names(reducedTerms_CT_AvsD_DN) == 'go'] <- 'id'
table$term <- NULL
reducedTerms_CT_AvsD_DN <- merge(reducedTerms_CT_AvsD_DN, table, by ="id")

########################################################################################### merge ####

order <- read.table("order.txt", sep = "\t")
colnames(order) <- c("term","category")

reducedTerms_AD_AvsB_UP$comparison <- "UP_AD_AvsB"
reducedTerms_AD_AvsB_DN$comparison <- "DN_AD_AvsB"
reducedTerms_AD_AvsC_UP$comparison <- "UP_AD_AvsC"
reducedTerms_AD_AvsC_DN$comparison <- "DN_AD_AvsC"
reducedTerms_AD_AvsD_UP$comparison <- "UP_AD_AvsD"
reducedTerms_AD_AvsD_DN$comparison <- "DN_AD_AvsD"
reducedTerms_CT_AvsB_UP$comparison <- "UP_CT_AvsB"
reducedTerms_CT_AvsB_DN$comparison <- "DN_CT_AvsB"
reducedTerms_CT_AvsC_UP$comparison <- "UP_CT_AvsC"
reducedTerms_CT_AvsC_DN$comparison <- "DN_CT_AvsC"
reducedTerms_CT_AvsD_UP$comparison <- "UP_CT_AvsD"
reducedTerms_CT_AvsD_DN$comparison <- "DN_CT_AvsD"

df <- rbind(reducedTerms_AD_AvsB_UP,
            reducedTerms_AD_AvsB_DN,
            reducedTerms_AD_AvsC_UP,
            reducedTerms_AD_AvsC_DN,
            reducedTerms_AD_AvsD_UP,
            reducedTerms_AD_AvsD_DN,
            reducedTerms_CT_AvsB_UP,
            reducedTerms_CT_AvsB_DN,
            reducedTerms_CT_AvsC_UP,
            reducedTerms_CT_AvsC_DN,
            reducedTerms_CT_AvsD_UP,
            reducedTerms_CT_AvsD_DN)

df <- df[,c(1,7,9,10,11,12,14,15)]

simMatrix <- calculateSimMatrix(df[,1],
                                orgdb="org.Dm.eg.db",
                                ont="BP",
                                method="Rel")

reducedSimMatrix <- reduceSimMatrix(simMatrix, scores=NULL,
                                    threshold=threshold_all,
                                    orgdb="org.Dm.eg.db")

order <- unique(reducedSimMatrix[order(reducedSimMatrix[,2],decreasing=FALSE),][,3])

df <- df %>% pivot_longer(cols = contains("comparison"))

df$direction <- ifelse(grepl("UP", df$value), "0", "1")
df$tiss <- ifelse(grepl("AD", df$value), "ABD","H+T")
#df$value <- gsub('.{3}$', '', df$value)
#df$value <- gsub('^.{3}', '', df$value)

cluster <- reducedSimMatrix[,c(1,2)]
rownames(cluster) <- NULL
colnames(cluster) <- c("id","cluster")
df <- merge(df, cluster, by = "id")

df <- merge(df, order, by = "term")

figure <- ggplot(df, aes(x = value, y = term, colour=as.factor(direction), size = -log(as.numeric(p)), alpha = as.numeric(`sig/exp`))) + 
  theme_bw() +
  geom_point() +
  facet_grid(category ~ tiss, scales="free", space="free") +
  theme(panel.grid.minor = element_blank()) +
  scale_color_manual(values=c("red", "blue")) + 
  theme(axis.text.x = element_text(angle=90)) +
  theme(strip.background = element_rect(color="black", fill="white")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(axis.text=element_text(size=9)) 

ggsave("crema_de_panel.pdf", 
       width = 15,
       height = 30,
       figure)

