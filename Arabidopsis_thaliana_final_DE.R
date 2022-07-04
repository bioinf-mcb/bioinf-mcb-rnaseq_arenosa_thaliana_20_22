# memory.limit(size=16000000)
# options(expressions = 5e5)

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

setwd("D:/MCB/PIOTR_R/Thaliana_proper/")
library(limma)
library(DESeq2)
library(edgeR)
library(topGO)
library(dplyr)
library(xlsx)
library(org.At.tair.db)
library(jsonlite)
library(glmnetUtils)

Experiment <- read.csv("./counts/thaliana_root_counts.csv")
rownames(Experiment) <- Experiment$X
Experiment <- Experiment[,2:10]
View(Experiment)


samples <- substr(colnames(Experiment), 9, 9)
design_limma <- data.frame(Es = ifelse(samples == "E", 1, 0),
                           Rs = ifelse(samples == "R", 1, 0))

rownames(design_limma) <- colnames(Experiment)
design_limma
# Counts normalization and MVE computation ----
# keep = filterByExpr(as.matrix(Experiment), design_limma)
# keep["__no_feature"] <- FALSE
# keep["__ambiguous"] <- FALSE
# keep["__too_low_aQual"] <- FALSE
# keep["__not_aligned"] <- FALSE 
# keep["__alignment_not_unique"] <- FALSE
# transposed <- data.table::transpose(Experiment[keep,])
# colnames(transposed) <- rownames(Experiment[keep,])
# rownames(transposed) <- c("niE1", "niE2", "niE3",
#                           "niE4", "niE5", "niR2",
#                           "niR3", "niR4", "niR5")
# 
# library(outliers)
# outlier(transposed)
# mod <- lm(Type ~ ., data = transposed)
# mod <- glmnet(Type ~ ., data = transposed, use.model.frame=TRUE)
# cooksd <- cooks.distance(mod$glmnet)
# plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")
# abline(h = 4*mean(cooksd, na.rm=T), col="red")
# text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")

# outlier(Experiment[keep,])


DGElimma <- DGEList(counts = Experiment)
DGElimma <- calcNormFactors(DGElimma)
voom_limma <- voom(DGElimma, design_limma, plot = T)

# Fitting and p-value adjustment ----
contrasts <- makeContrasts(EvsR = Es-Rs,
                           levels = design_limma)
voom_fit <- lmFit(voom_limma, contrast = contrasts)
cf_limma <- contrasts.fit(voom_fit, contrasts)
eBayes_limma <- eBayes(cf_limma, proportion = 0.01)

TMM_voom_counts <- topTable(eBayes_limma,
                            number = Inf,
                            adjust.method = "BH",
                            sort.by = "none",
                            confint = T)
write.csv(x = TMM_voom_counts, file = "Athaliana_voom.csv")

# Plotting of significant genes
significant_limma <- TMM_voom_counts[TMM_voom_counts$adj.P.Val <= 0.05,]
dim(significant_limma)
LogFC_boolvec <- abs(significant_limma$logFC) >= 1
significant_limma <- significant_limma[LogFC_boolvec,]
AveExp_Q1 <- quart(significant_limma$AveExpr)[1]
AveExp_boolvec <- significant_limma$AveExpr >= AveExp_Q1
significant_limma <- significant_limma[AveExp_boolvec,]
dim(significant_limma)
write.csv(file = "significant_limma.csv", x = significant_limma)



pchs <- rep('.', dim(TMM_voom_counts)[1])
colors <- rep('black', dim(TMM_voom_counts)[1])
names(pchs) <- names(colors) <- rownames(TMM_voom_counts)
selected_transcripts <- rownames(significant_limma)

pchs[selected_transcripts] <- "x"

colors[selected_transcripts] <- "purple"

plot(TMM_voom_counts$AveExpr, TMM_voom_counts$logFC,
     xlab = 'Average Expression', ylab = "logFC",
     main = "Limma MA plot", pch = pchs, col = colors)

write.csv(voom_limma$E, "counts_thaliana.csv")

potential <- TMM_voom_counts[TMM_voom_counts$adj.P.Val <= 0.05,]
LogFC_boolvec <- abs(potential$logFC)>=1
potential <- potential[LogFC_boolvec,]
AveExp_Q1 <- quart(potential$AveExpr)[1]
AveExp_boolvec <- potential$AveExpr < AveExp_Q1
potential <- potential[AveExp_boolvec,]
write.csv(file = "significant_low_copy_root_limma.csv", x = potential)






