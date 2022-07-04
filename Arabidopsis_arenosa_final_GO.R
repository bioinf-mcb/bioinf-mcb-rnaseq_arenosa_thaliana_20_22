# Sourcing Differential Expression script ----
source("Arabidopsis_arenosa_final_DE.R")
setwd("D:/MCB/PIOTR_R/LEGACY/Arenosa_proper/")
library(org.At.tair.db)
library(tidyr)
library(ggplot2)
library(ggpubr)
options(scipen = 0, digits = 2)
# Definition of gene filtering function for topGO ----
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


topDiffGenes <- function(allGene, thr=0.05)
{
  return(allGene <= thr);
}

IDtoTerm <- function(GOlist){
  GO_names <- names(GOlist)
  for(i in 1:length(GO_names)){
    GO_names[i] <- AnnotationDbi::Term(GOTERM[GO_names[i]])
  }
  names(GOlist) <- GO_names
  return(GOlist)
}


TermConversion <- function(GOvec){
  for(i in 1:length(GOvec)){
    GOvec[i] <- AnnotationDbi::Term(GOTERM[GOvec[i]])
  }
  return(GOvec)
}




# topGO Biological Process for limma (shoot) ----
genelist_shoot_limma <- limma_countsTMMvoom_shoot$adj.P.Val
names(genelist_shoot_limma) <- rownames(limma_countsTMMvoom_shoot)
names(genelist_shoot_limma)
GO_shoot_limma <- new("topGOdata",
                      description = "Shoot limma BP", ontology = "BP",
                      allGenes = genelist_shoot_limma, geneSel = topDiffGenes,
                      nodeSize = 1, #1000
                      annot = annFUN.org, mapping = "org.At.tair.db")

limma_shoot_Fischer <- runTest(GO_shoot_limma, algorithm = "elim", statistic = "fisher")

limma_res_shoot <- GenTable(GO_shoot_limma, elimFisher = limma_shoot_Fischer,
                            orderBy = "elimFisher", ranksOf = "elimFisher",
                            topNodes = 500, numChar=1000) #20

View(limma_res_shoot)
# # 
# res_shoot <- gather(limma_res_shoot, key = "Legend", value = "Value", Significant, Expected)
# res_shoot$elimFisher <- sapply(X = res_shoot$elimFisher, as.numeric)
# res_shoot["elimFisher"] <- sapply(res_shoot["elimFisher"], format, scientific = T)
# res_shoot[21:40,"elimFisher"] <- c(rep("", 20))
# # res_shoot$Term <- sapply(strwrap(res_shoot$Term, width = 20, simplify = FALSE), paste, collapse="\n")
# # View(res_shoot)
# colnames(res_shoot)[1] <- "GO Term ID"
# colnames(res_shoot)[2] <- "GO Term"
# 
# ggbarplot(res_shoot, "GO Term", "Value",
#                   fill = "Legend",
#                   color = "Legend",
#                   palette = "Paired",
#                   # label = res_shoot[1:10,"elimFisher"],
#                   position = position_dodge2(0.5)) +
#   geom_text(aes(label=elimFisher), check_overlap = T, nudge_y = 4, na.rm = T) +
#   coord_flip()


write.csv(limma_res_shoot, "Arenosa_shoot_BP.csv")

write.xlsx(limma_res_shoot, file = "./results/BiologicalProcess.xlsx",
           sheetName = "limma_shoot", append = TRUE)

sig_limma_shoot <- significant_shoot_limma
genes_limma_shoot <- genesInTerm(GO_shoot_limma)



genes_limma_shoot <- IDtoTerm(genes_limma_shoot)



limma_shoot_json <- toJSON(lapply(genes_limma_shoot, function(x) x[x %in% rownames(sig_limma_shoot)]),pretty = TRUE, auto_unbox = TRUE)
fileConn<-file("./results/limma_shoot_BP.json")
write(limma_shoot_json, fileConn)
close(fileConn)

# topGO Biological Process for limma (root) ----
genelist_root_limma <- limma_countsTMMvoom_root$adj.P.Val
names(genelist_root_limma) <- rownames(limma_countsTMMvoom_root)
GO_root_limma <- new("topGOdata",
                     description = "Root limma BP", ontology = "BP",
                     allGenes = genelist_root_limma, geneSel = topDiffGenes,
                     nodeSize = 1, #2000
                     annot = annFUN.org, mapping = "org.At.tair.db")

limma_root_Fischer <- runTest(GO_root_limma, algorithm = "elim", statistic = "fisher")

limma_res_root <- GenTable(GO_root_limma, elimFisher = limma_root_Fischer,
                           orderBy = "elimFisher", ranksOf = "elimFisher",
                           topNodes = 500, numChar=1000) #25

# 
View(limma_res_root)
# 
# res_shoot <- gather(limma_res_root, key = "Legend", value = "Value", Significant, Expected)
# res_shoot[c(1,26), "elimFisher"] <- "1e-30"
# res_shoot$elimFisher <- sapply(X = res_shoot$elimFisher, as.numeric)
# res_shoot["elimFisher"] <- sapply(res_shoot["elimFisher"], format, scientific = T)
# res_shoot[26:50,"elimFisher"] <- c(rep("", 25))
# # res_shoot$Term <- sapply(strwrap(res_shoot$Term, width = 20, simplify = FALSE), paste, collapse="\n")
# # View(res_shoot)
# colnames(res_shoot)[1] <- "GO Term ID"
# colnames(res_shoot)[2] <- "GO Term"
# ggbarplot(res_shoot, "GO Term", "Value",
#           fill = "Legend",
#           color = "Legend",
#           palette = "Paired",
#           # label = res_shoot[1:10,"elimFisher"],
#           position = position_dodge2(0.5)) +
#   geom_text(aes(label=elimFisher), check_overlap = T, nudge_y = 200 , na.rm = T) +
#   coord_flip()
# # 
# # 








View(limma_res_root)










write.csv(limma_res_root, "Arenosa_root_BP.csv")

write.xlsx(limma_res_root, file = "./results/BiologicalProcess.xlsx",
           sheetName = "limma_root", append = TRUE)

sig_limma_root <- significant_root_limma
rownames(sig_limma_root) <- rownames(sig_limma_root)
genes_limma_root <- genesInTerm(GO_root_limma)
genes_limma_root <- IDtoTerm(genes_limma_root)
limma_root_json <- toJSON(lapply(genes_limma_root, function(x) x[x %in% rownames(sig_limma_root)]),pretty = TRUE, auto_unbox = TRUE)
fileConn<-file("./results/limma_root_BP.json")
write(limma_root_json, fileConn)
close(fileConn)
# Loading libraries for plotting of heatmaps and defining GO terms of interest ----
library(pheatmap)
library(RColorBrewer)




auksyna <- c("GO:0009733", "GO:0009850")
auksyna <- TermConversion(auksyna)
korzen <- c("GO:0048765", "GO:0010054", "GO:0048767", "GO:0048364")
korzen <- TermConversion(korzen)
pozostale <- c("GO:0006833", "GO:0009828", "GO:0071555", "GO:0009913", "GO:0048527")
pozostale <- TermConversion(pozostale)
nutrienty <- c("GO:0016036", "GO:0006995", "GO:0006817", "GO:0015706")
nutrienty <- TermConversion(nutrienty)
metale <- c("GO:0006826", "GO:0000041", "GO:0046686", "GO:0010043", "GO:0006865")
metale <- TermConversion(metale) 
etylen <- c("GO:0009723", "GO:0009693", "GO:0071369")
etylen <- TermConversion(etylen)
inne <- c("GO:0006950", "GO:0042221", "GO:0009628")
inne <- TermConversion(inne)
wszystkie <- union(names(genes_limma_root), names(genes_limma_shoot))
etylen
wszystkie


roi <- limma_res_root$GO.ID[c(1,3,5,8,18)]
roi <- TermConversion(roi)

# Definition of plotting function of heatmaps ----
CreateSetOfPlots <- function(normcounts, golist, deresult, gogenes, savedir, root = T){
  for(i in 1:length(golist)){
    print(golist[i])
    all_genes <- as.vector(unlist(gogenes[golist[i]]))
    sig_genes <- intersect(rownames(deresult), all_genes)
    print(sig_genes)
    sig_counts <- normcounts[sig_genes,]
    goname <- gsub(":","",golist[i])
    if(root){
      fname <- paste0("./heatmaps/",savedir, "/root/", goname, ".png")
      csv_name <- paste0("./heatmaps/",savedir, "/root/", goname, ".csv")
    }else{
      fname <- paste0("./heatmaps/",savedir, "/shoot/", goname, ".png")
      csv_name <- paste0("./heatmaps/",savedir, "/shoot/", goname, ".csv")
    }
    print(fname)
    out <- tryCatch(
      {
        pheatmap(sig_counts, display_numbers = F, legend = F, angle_col = 0,
                 # border_color = "gray",
                 main = goname, filename = fname, cell_height = 4, height = 6,
                 fontsize = 4,  number_color = "black", cell_width = 3, width = 5,
                 color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100))
        write.csv(sig_counts, csv_name)
        dev.off()
    }, error = function(cond){
      print("Not enough significant genes")
    }, warning = function(cond){
      print("Not enough significant genes")
    }, finally = function(cond){
      print("Not enough significant genes")
    }
    )
    if(length(dev.list()) > 0){
      for (i in dev.list()[1]:dev.list()[length(dev.list())]) {
        # dev.off()
      }
    }
  }
}

# Data transformation ----
normalized_shoot <- v_shoot$E
Mean_E_shoot <- as.vector(rowMeans(normalized_shoot[,1:5]))
Mean_TM_shoot <- as.vector(rowMeans(normalized_shoot[,6:10]))
normalized_mean_shoot <- data.frame(Mean_E_shoot, Mean_TM_shoot)
rownames(normalized_mean_shoot) <- rownames(normalized_shoot)

normalized_root <- v_root$E
Mean_E_root <- as.vector(rowMeans(normalized_root[,1:5]))
Mean_TM_root <- as.vector(rowMeans(normalized_root[,6:9]))
normalized_mean_root <- data.frame(Mean_E_shoot, Mean_TM_shoot)
rownames(normalized_mean_root) <- rownames(normalized_shoot)

colnames(normalized_mean_shoot) <- colnames(normalized_mean_root) <- c("E+", "E-")

rootDE <- significant_root_limma
shootDE <- significant_shoot_limma

gogenes_shoot <- genes_limma_shoot
sink("gogenes_shoot.txt")
print(gogenes_shoot)
sink()
gogenes_root <- genes_limma_root
sink("gogenes_root.txt")
print(gogenes_root)
sink()
# Plotting ----
setwd("./results/GO/")

significant_vec <- rootDE$adj.P.Val <= 0.05
rootDE <- rootDE[significant_vec,]
fold_vec <- abs(rootDE$logFC) >= 2
rootDE <- rootDE[fold_vec,]
q1_filter <- quart(rootDE$AveExpr)[1]
q1_vec <- rootDE$AveExpr >= q1_filter
rootDE <- rootDE[q1_vec,]

significant_vec <- shootDE$adj.P.Val <= 0.05
shootDE <- shootDE[significant_vec,]
fold_vec <- abs(shootDE$logFC) >= 2
shootDE <- shootDE[fold_vec,]
q1_filter <- quart(shootDE$AveExpr)[1]
q1_vec <- shootDE$AveExpr >= q1_filter
shootDE <- shootDE[q1_vec,]

mapIds(org.At.tair.db, keys = rownames(rootDE), column = "PMID", keytype = "TAIR", mutiVals = "first")

CreateSetOfPlots(normalized_mean_shoot, roi, shootDE, gogenes_shoot, "terms_of_interest", root = F)
CreateSetOfPlots(normalized_mean_root, roi, rootDE, gogenes_root, "terms_of_interest")




CreateSetOfPlots(normalized_mean_shoot, auksyna, shootDE, gogenes_shoot, "auksyna", root = F)
CreateSetOfPlots(normalized_mean_root, auksyna, rootDE, gogenes_root, "auksyna")

CreateSetOfPlots(normalized_mean_shoot, korzen, shootDE, gogenes_shoot, "korzen", root = F)
CreateSetOfPlots(normalized_mean_root, korzen, rootDE, gogenes_root, "korzen")

CreateSetOfPlots(normalized_mean_shoot, pozostale, shootDE, gogenes_shoot, "pozostale", root = F)
CreateSetOfPlots(normalized_mean_root, pozostale, rootDE, gogenes_root, "pozostale")

CreateSetOfPlots(normalized_mean_shoot, nutrienty, shootDE, gogenes_shoot, "nutrienty", root = F)
CreateSetOfPlots(normalized_mean_root, nutrienty, rootDE, gogenes_root, "nutrienty")

CreateSetOfPlots(normalized_mean_shoot, metale, shootDE, gogenes_shoot, "metale", root = F)
CreateSetOfPlots(normalized_mean_root, metale, rootDE, gogenes_root, "metale")

CreateSetOfPlots(normalized_mean_shoot, etylen, shootDE, gogenes_shoot, "etylen", root = F)
CreateSetOfPlots(normalized_mean_root, etylen, rootDE, gogenes_root, "etylen")

CreateSetOfPlots(normalized_mean_shoot, wszystkie, shootDE, gogenes_shoot, "wspolne", root = F)
CreateSetOfPlots(normalized_mean_root, wszystkie, rootDE, gogenes_root, "wspolne")

CreateSetOfPlots(normalized_mean_shoot, inne, shootDE, gogenes_shoot, "inne", root = F)
CreateSetOfPlots(normalized_mean_root, inne, rootDE, gogenes_root, "inne")

