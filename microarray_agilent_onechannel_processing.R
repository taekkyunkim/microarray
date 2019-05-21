############ 1. Preparation ############
# 1. Install "limma" if it is not installed yet.
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma")
BiocManager::install("corrplot")
BiocManager::install("ComplexHeatmap")
BiocManager::install("circlize")

############ 2. Normalization ############
# 1. Remove all stuffs in the R environment
rm(list=ls())

# 2. Set working directory
setwd('XXXXXXXXX')  #######  specify the working directory which includes input data files
path <- getwd()

# 3. Call libraries to be needed for analysis
library(limma)

# 4. List up file names
filelist <- list.files(path = ".", pattern = "txt")

# 5. Quantile normalization of gMeanSignals
gMean <- read.maimages(files=filelist, source="agilent", columns=list(G="gMeanSignal"))
annot <- gMean$genes
gMeanQnorm <- normalizeBetweenArrays(gMean, method="quantile")
write.table(gMeanQnorm, 'gMeanQuantileNormamlizedData.txt', sep="\t", row.names=FALSE)

# 6. Quantile normalization of gProcessedSignals
gProcessed <- read.maimages(files=filelist, source="agilent", columns=list(G="gProcessedSignal"))
annot <- gProcessed$genes
gProcessedQnorm <- normalizeBetweenArrays(gProcessed, method="quantile")
write.table(gProcessedQnorm, 'gProcesseduantileNormamlizedData.txt', sep="\t", row.names=FALSE)

# 7. Present & Absent Call
data <- gProcessedQnorm  ####### specify dataset
present_absent <- annot$ControlType == 0 & rowMeans(data$E) > mean(data$E)

# 8. Correlation coefficient
library(corrplot)
data <- gProcessedQnorm  ####### specify dataset
M<-cor(data$E)
col1 <- colorRampPalette(c("red", "white", "blue")) 
corrplot(M, method="color", outline = FALSE, diag = TRUE, is.corr = FALSE, tl.cex = 0.5, col = col1(100))

# 9. Scatter plot (draw one pair only at a time due to the running time issue)
data <- gProcessedQnorm  ####### specify dataset
ids_1 <- 1  ####### specify the column number of the first sample you want
ids_2 <- 2  ####### specify the column number of the second sample you want
plot(data$E[,ids_1], data$E[,ids_2], xlab=filelist[ids_1], ylab=filelist[ids_2], pch=1)


############ 3. Differential expression analysis ############
# 1. Set groups
data <- gProcessedQnorm  #######  choose which dataset you will use, gMeanQnorm or gProcessQnorm
case_id <- c(1,2)  #######  choose which columns you will use for case group
control_id <- c(3,4)  #######  choose which columns you will use for control group
expr <- data$E[, c(case_id, control_id)]
case_group <- c(rep(1, length(case_id)), rep(0, length(control_id)))
control_group <- c(rep(0, length(case_id)), rep(1, length(control_id)))

# 2. Generate design matrix
design <- cbind(Case=case_group,Control=control_group)
fit<-lmFit(expr, design);
cont.matrix <- makeContrasts(CasevsControl=Case-Control, levels=design)
fit <- contrasts.fit(fit, cont.matrix)
fit <- eBayes(fit);

# 3. Save results
result <- topTable(fit, number="all", adjust="BH", sort="none")
annot_result <- cbind(annot, result, present_absent)
write.table(annot_result, "LIMMA_output.txt", sep="\t", row.names=FALSE)

# 4. Volcano plot (p-value and fold change)
PValue <- annot_result$P.Value
#FDR <- annot_result$adj.P.Val
log2Fold <- annot_result$logFC
plot(log2Fold, -log10(PValue))

# 5. Pvalue distribution
PValue <- annot_result$P.Value
hist(PValue, breaks=100)

# 6. Differentially expressed genes
pvalue_cutoff <- 0.01 ####### specify cutoff
log2_fold_cutoff <- 0.585  #######  specify cutoff
data <- gProcessedQnorm  #######  specify dataset

present_id <- which(annot_result$ControlType == 0 & rowMeans(data$E) > mean(data$E))
PValue <- annot_result$P.Value
log2Fold <- annot_result$logFC
upid <- which(PValue < pvalue_cutoff & log2Fold > log2_fold_cutoff)
upid <- intersect(upid, present_id)
dnid <- which(PValue < pvalue_cutoff & log2Fold < -log2_fold_cutoff)
dnid <- intersect(dnid, present_id)

summary_table <- cbind(length(upid),
                       length(unique(annot_result$GeneName[upid])),
                       length(dnid),
                       length(unique(annot_result$GeneName[dnid])))
colnames(summary_table) <- c('No. of Up-regulated Probes', 
                             'No. of Up-regulated Genes',
                             'No. of Down-regulated Probes',
                             'No. of Down-regulated Genes')

barplot(summary_table[,c(1,3)], ylab='No. of probes')
barplot(summary_table[,c(2,4)], ylab='No. of genes')

deg_table <- annot_result[c(upid, dnid), ]
write.table(deg_table, "LIMMA_DEGs.txt", sep="\t", row.names=FALSE)

# 7. Generate heatmap for DEGs
library(ComplexHeatmap)
library(circlize)

mncn <- function(x) {
  n = nrow(x)
  ones = rep(1, n)
  H = diag(n) - (1/n) * (ones %*% t(ones))
  H %*% x
}

expr <- data$E[, c(case_id, control_id)]
degids <- c(upid, dnid)

tm <- expr[degids, ]
normMat <- t(mncn(t(tm)))
type <- gsub("s\\d+_", "", filelist[c(case_id,control_id)])
ha <- HeatmapAnnotation(df = data.frame(type = type))
Heatmap(normMat, name = "expression", 
        km = 1, 
        col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")),
        top_annotation = ha, 
        top_annotation_height = unit(4, "mm"), 
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 10),
        cluster_columns = FALSE,
        cluster_rows = TRUE,
        show_row_names = FALSE, 
        show_row_dend = TRUE,
        show_column_names = FALSE)



