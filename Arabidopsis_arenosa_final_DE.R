# Import of libraries ----

library(limma)
library(DESeq2)
library(edgeR)
library(topGO)
library(dplyr)
library(xlsx)
library(org.At.tair.db)
library(jsonlite)
library(outliers)
library(readr)


quart <- function(x) {
        x <- sort(x)
        n <- length(x)
        m <- (n+1)/2
        if (floor(m) != m) {
                l <- m-1/2; u <- m+1/2
        } else {
                l <- m-1; u <- m+1
        }
        c(Q1=median(x[1:l]), Q3=median(x[u:n]))
}

setwd("D:/MCB/PIOTR_R/LEGACY/Arenosa_proper/")
# Loading generated count tables and parsing to one big dataframe ----

counts_root <- read.csv("counts/arenosa_root_counts.csv")
rownames(counts_root) <- counts_root$X
counts_root <- counts_root[2:10]
counts_shoot <- read.csv("counts/arenosa_shoot_counts.csv")
rownames(counts_shoot) <- counts_shoot$X
counts_shoot <- counts_shoot[2:11]

# Splitting data into Shoot and Root ----

Shoot_df <- counts_shoot
Shoot_class <- DataFrame(condition = factor(rep(c("E","TM"), c(5,5))))
head(Shoot_df)
print(Shoot_class)

Root_df <- counts_root
Root_class <- DataFrame(condition = factor(rep(c("E","TM"), c(5,4))))
head(Root_df)
print(Root_class)

samples_shoot <- substr(colnames(Shoot_df), 7, 7)
samples_root <- substr(colnames(Root_df), 7, 7)


design_shoot <-  data.frame(Es = ifelse(samples_shoot=="E",1,0),
                            Ts=ifelse(samples_shoot=="T",1,0))

design_root <-  data.frame(Es = ifelse(samples_root=="E",1,0),
                           Ts=ifelse(samples_root=="T",1,0))

rownames(design_shoot) <- colnames(Shoot_df)
rownames(design_root) <- colnames(Root_df)

# keep = filterByExpr(as.matrix(Shoot_df), design_shoot)
# keep["__no_feature"] <- FALSE
# keep["__ambiguous"] <- FALSE
# keep["__too_low_aQual"] <- FALSE
# keep["__not_aligned"] <- FALSE
# keep["__alignment_not_unique"] <- FALSE
# 
# Shoot_df <- Shoot_df[keep,]
# outlier(Shoot_df)
# 
# 
# keep = filterByExpr(as.matrix(Root_df), design_root)
# keep["__no_feature"] <- FALSE
# keep["__ambiguous"] <- FALSE
# keep["__too_low_aQual"] <- FALSE
# keep["__not_aligned"] <- FALSE
# keep["__alignment_not_unique"] <- FALSE
# 
# Root_df <- Root_df[keep,]
# outlier(Root_df)
dim(Shoot_df)
dim(Root_df)


# DESeq2 shoot ----------

coldata_shoot <- data.frame(condition = samples_shoot)
coldata_root <- data.frame(condition = samples_root)
library(DESeq2)
dds_shoot <-  DESeqDataSetFromMatrix(countData = Shoot_df,
                                     colData = coldata_shoot,
                                     design = formula(~condition))
dds_shoot <- DESeq(dds_shoot)
resultsNames(dds_shoot)
results_shoot <- results(dds_shoot)
results_shoot <- results_shoot[!is.na(results_shoot$padj),]
significant_shoot_DESeq2 <- results_shoot[results_shoot$padj <= 0.05,]
dim(significant_shoot_DESeq2)

pchs <- rep('.', dim(Shoot_df)[1])
colors <- rep('black', dim(Shoot_df)[1])
names(pchs) <- names(colors) <- rownames(Shoot_df)
selected_transcripts <- rownames(significant_shoot_DESeq2)
pchs[selected_transcripts] <- "x"
colors[selected_transcripts] <- "purple"
plot(results_shoot$baseMean, results_shoot$log2FoldChange,
     xlab = 'Average Expression', ylab = "logFC",
     main = "DESeq2 shoot MA plot", pch = pchs, col = colors,
     xlim = c(0,2e03))
# DESeq2 root ----
dds_root <-  DESeqDataSetFromMatrix(countData = Root_df,
                                    colData = coldata_root,
                                    design = formula(~condition))
dds_root <- DESeq(dds_root)
resultsNames(dds_root)
results_root <- results(dds_root)
results_root <- results_root[!is.na(results_root$padj),]
significant_root_DESeq2 <- results_root[results_root$padj <= 0.05,]
dim(significant_root_DESeq2)

pchs <- rep('.', dim(Root_df)[1])
colors <- rep('black', dim(Root_df)[1])
names(pchs) <- names(colors) <- rownames(Root_df)
selected_transcripts <- rownames(significant_root_DESeq2)
pchs[selected_transcripts] <- "x"
colors[selected_transcripts] <- "purple"
plot(results_root$baseMean, results_root$log2FoldChange,
     xlab = 'Average Expression', ylab = "logFC",
     main = "DESeq2 root MA plot", pch = pchs, col = colors,
     xlim = c(0,2e03))

# edgeR shoot ----
shoot_dge <- DGEList(counts = Shoot_df)
shoot_dge <- calcNormFactors(shoot_dge)
shoot_dge <- estimateDisp(shoot_dge, design_shoot)
shoot_fit <- glmQLFit(shoot_dge, design_shoot)
shoot_lrt <- glmLRT(shoot_fit)
shoot_contrasts <- makeContrasts(EvsT = Es-Ts,
                                 levels = design_shoot)
shoot_qlf <- glmQLFTest(shoot_fit, contrast = shoot_contrasts)
topTags(shoot_qlf, n = 30)


shoot_results_edgeR <- as.data.frame(shoot_qlf$table)
shoot_results_edgeR

significant_shoot_edgeR <- shoot_results_edgeR[shoot_results_edgeR$PValue <= 0.05,]

pchs <- rep('.', dim(Shoot_df)[1])
colors <- rep('black', dim(Shoot_df)[1])
names(pchs) <- names(colors) <- rownames(Shoot_df)
selected_transcripts <- rownames(significant_shoot_edgeR)
pchs[selected_transcripts] <- "x"
colors[selected_transcripts] <- "purple"
plot(shoot_results_edgeR$logCPM, shoot_results_edgeR$logFC,
     xlab = 'Average Expression', ylab = "logFC",
     main = "edgeR shoot MA plot", pch = pchs, col = colors)
# edgeR root ----
root_dge <- DGEList(counts = Root_df)
root_dge <- calcNormFactors(root_dge)
root_dge <- estimateDisp(root_dge, design_root)
root_fit <- glmQLFit(root_dge, design_root)
root_lrt <- glmLRT(root_fit)
root_contrasts <- makeContrasts(EvsT = Es-Ts,
                                levels = design_root)
root_qlf <- glmQLFTest(root_fit, contrast = root_contrasts)
topTags(root_qlf, n = 30)


root_results_edgeR <- as.data.frame(root_qlf$table)
root_results_edgeR

significant_root_edgeR <- root_results_edgeR[root_results_edgeR$PValue <= 0.05,]

pchs <- rep('.', dim(Root_df)[1])
colors <- rep('black', dim(Root_df)[1])
names(pchs) <- names(colors) <- rownames(Root_df)
selected_transcripts <- rownames(significant_root_edgeR)
pchs[selected_transcripts] <- "x"
colors[selected_transcripts] <- "purple"
plot(root_results_edgeR$logCPM, root_results_edgeR$logFC,
     xlab = 'Average Expression', ylab = "logFC",
     main = "edgeR root MA plot", pch = pchs, col = colors)

# limma shoot ----
shoot_dge_limma <- DGEList(counts = Shoot_df)
shoot_dge_limma <- calcNormFactors(shoot_dge_limma)
v_shoot <- voom(shoot_dge_limma, design_shoot, plot = T)
v_shoot_fit <- lmFit(v_shoot, design_shoot)
cf_shoot <- contrasts.fit(v_shoot_fit, shoot_contrasts)
eBayes_shoot <- eBayes(cf_shoot, proportion = 0.01)
limma_countsTMMvoom_shoot <- topTable(eBayes_shoot,
                                      coef = "Es",
                                      number = Inf,
                                      adjust.method = "BH",
                                      sort.by="none",
                                      confint = T)



significant_shoot_limma <- limma_countsTMMvoom_shoot[limma_countsTMMvoom_shoot$adj.P.Val <= 0.05,]
LogFC_boolvec <- abs(significant_shoot_limma$logFC)>=2
significant_shoot_limma <- significant_shoot_limma[LogFC_boolvec,]
AveExp_Q1 <- quart(significant_shoot_limma$AveExpr)[1]
AveExp_boolvec <- significant_shoot_limma$AveExpr >= AveExp_Q1
significant_shoot_limma <- significant_shoot_limma[AveExp_boolvec,]
dim(significant_shoot_limma)

pchs <- rep('.', dim(Shoot_df)[1])
colors <- rep('black', dim(Shoot_df)[1])
names(pchs) <- names(colors) <- rownames(Shoot_df)
selected_transcripts <- rownames(significant_shoot_limma)
pchs[selected_transcripts] <- "x"
colors[selected_transcripts] <- "red"
plot(limma_countsTMMvoom_shoot$AveExpr, limma_countsTMMvoom_shoot$logFC,
     xlab = 'Average Expression', ylab = "logFC",
     main = "AAE+ vs AAE- MA plot for shoot", pch = pchs, col = colors)




diffex <- limma_countsTMMvoom_shoot
diffex <- diffex[,c("AveExpr", "logFC", "adj.P.Val")]
name <- rownames(diffex)
diffex <- cbind(name, diffex)
colnames(diffex) <- c("name", "baseMeanLog2", "log2FoldChange", "padj")


library(ggpubr)
p <- ggpubr::ggmaplot(data = diffex,
                 size = 1, fdr = 0.05, fc = 0,
                 legend = "top", top = 0)
p

p_data <- p$data

shoot_name_vec <- p_data[p_data$sig != "NS","name"]

limma_countsTMMvoom_shoot <- limma_countsTMMvoom_shoot[rownames(limma_countsTMMvoom_shoot) %in% name_vec,]


# limma root ----
root_dge_limma <- DGEList(counts = Root_df)
root_dge_limma <- calcNormFactors(root_dge_limma)
v_root <- voom(root_dge_limma, design_root, plot = T)
v_root_fit <- lmFit(v_root, design_root)
cf_root <- contrasts.fit(v_root_fit, root_contrasts)
eBayes_root <- eBayes(cf_root, proportion = 0.01)
limma_countsTMMvoom_root <- topTable(eBayes_root,
                                     number = Inf,
                                     adjust.method = "BH",
                                     sort.by="none",
                                     confint = T)
vp <- volcanoplot(eBayes_root, coef = 1, highlight = T)
# write.csv(file = "Aarenosa_voom.csv", x = limma_countsTMMvoom_root)

significant_root_limma <- limma_countsTMMvoom_root[limma_countsTMMvoom_root$adj.P.Val <= 0.05,]
LogFC_boolvec <- abs(significant_root_limma$logFC)>=2
significant_root_limma <- significant_root_limma[LogFC_boolvec,]
AveExp_Q1 <- quart(significant_root_limma$AveExpr)[1]
AveExp_boolvec <- significant_root_limma$AveExpr >= AveExp_Q1
significant_root_limma <- significant_root_limma[AveExp_boolvec,]
dim(significant_root_limma)

pchs <- rep('.', dim(Root_df)[1])
colors <- rep('black', dim(Root_df)[1])
names(pchs) <- names(colors) <- rownames(Root_df)
selected_transcripts <- rownames(significant_root_limma)
pchs[selected_transcripts] <- "x"
colors[selected_transcripts] <- "red"
plot(limma_countsTMMvoom_root$AveExpr, limma_countsTMMvoom_root$logFC,
     xlab = 'Average Expression', ylab = "logFC",
     main = "AAE+ vs AAE- MA plot for root", pch = pchs, col = colors)


write.csv(file = "significant_root_limma.csv", x = significant_root_limma)
write.csv(file = "significant_shoot_limma.csv", x = significant_shoot_limma)

potential <- limma_countsTMMvoom_root[limma_countsTMMvoom_root$adj.P.Val <= 0.05,]
LogFC_boolvec <- abs(potential$logFC)>=2
potential <- potential[LogFC_boolvec,]

dim(potential)

AveExp_Q1 <- quart(potential$AveExpr)[1]
AveExp_boolvec <- potential$AveExpr < AveExp_Q1
potential <- potential[AveExp_boolvec,]
dim(potential)
write.csv(file = "significant_low_copy_root_limma.csv", x = potential)



diffex <- limma_countsTMMvoom_root
diffex <- diffex[,c("AveExpr", "logFC", "adj.P.Val")]
name <- rownames(diffex)
rownames(diffex) <- name
colnames(diffex) <- c("baseMeanLog2", "log2FoldChange", "padj")

diffex
library(ggpubr)
p <- ggpubr::ggmaplot(data = diffex,
                 fdr = 0.05, fc = 0, size = 2,
                 legend = "top", top = 0)
p
p_data <- p$data

root_name_vec <- p_data[p_data$sig != "NS","name"]

limma_countsTMMvoom_root <- limma_countsTMMvoom_root[rownames(limma_countsTMMvoom_root) %in% name_vec,]

ggpubr::ggmaplot(diff_express, main = expression("Group 1" %->% "Group 2"),
         fdr = 0.05, fc = 2, size = 0.4,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         # genenames = as.vector(diff_express$name),
         legend = "top", top = 0,
         font.label = c("bold", 11),
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_minimal())
data(diff_express)


View(diffex)
View(diff_express)
