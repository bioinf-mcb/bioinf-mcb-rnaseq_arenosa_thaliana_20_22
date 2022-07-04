# Sourcing DE script ----
source("Arabidopsis_thaliana_final_DE.R")

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
#Setting up essentials ----
genelist_limma <- TMM_voom_counts$adj.P.Val
names(genelist_limma) <- rownames(TMM_voom_counts)

topDiffGenes <- function(allGene, thr = 0.05)
{
  return(allGene <= thr);
}

# topGO Biological process ----
topGO_limma <- new("topGOdata",
                   description = "Limma BP", ontology = "BP",
                   allGenes = genelist_limma, geneSel = topDiffGenes,
                   nodeSize = 1,
                   annot = annFUN.org, mapping = "org.At.tair.db")

limma_Fisher <- runTest(topGO_limma, algorithm = "elim", statistic = "fisher")

GO_res <- GenTable(topGO_limma, elimFisher = limma_Fisher,
                   orderBy = "elimFisher", ranksOf = "elimFisher",
                   topNodes = 500)

tail(GO_res)

write.csv(GO_res, "Thaliana_BP.csv")

genes_GO_limma <- genesInTerm(topGO_limma)
genes_GO_limma <- IDtoTerm(genes_GO_limma)
save(genes_GO_limma,file ="genes_GO_limma.RData")

sink("genes_limma.txt")
print(genes_GO_limma)
sink()

# Plotting ----
library(pheatmap)
library(RColorBrewer)
auxin <- c("GO:0009733", "GO:0009850")
auxin <- TermConversion(auxin)
root <- c("GO:0048765", "GO:0010054", "GO:0048767", "GO:0048364")
root <- TermConversion(root)
other <- c("GO:0006833", "GO:0009828", "GO:0071555", "GO:0009913", "GO:0048527")
other <- TermConversion(other)
nutrients <- c("GO:0016036", "GO:0006995", "GO:0006817", "GO:0015706")
nutrients <- TermConversion(nutrients)
metals <- c("GO:0006826", "GO:0000041", "GO:0046686", "GO:0010043", "GO:0006865")
metals <- TermConversion(metals)
ethylen <- c("GO:0009723", "GO:0009693", "GO:0071369")
ethylen <- TermConversion(ethylen)
wszystkie <- names(genes_GO_limma)

normalized_counts <- voom_limma$E
MeanE <- as.vector(rowMeans(normalized_counts[,1:5]))
MeanR <- as.vector(rowMeans(normalized_counts[,6:9]))
normalized_mean <- data.frame(MeanE, MeanR)
colnames(normalized_mean) <- c("Mean normalized expression E", "Mean normalized expression R")
rownames(normalized_mean) <- rownames(normalized_counts)
head(normalized_mean)

CreateSetOfPlots <- function(normcounts, golist, deresult, gogenes, savedir){
  for(i in 1:length(golist)){
    print(golist[i])
    all_genes <- as.vector(unlist(gogenes[golist[i]]))
    sig_genes <- intersect(rownames(deresult), all_genes)
    sig_counts <- normcounts[sig_genes,]
    goname <- gsub(":","",golist[i])
    #     rownames(sig_counts) <- gsub("gene:", "", rownames(sig_counts))
    fname <- paste0(savedir, goname, ".png")
    csv_name <- paste0(savedir, goname, ".csv")
    out <- tryCatch(
      {
        pheatmap(sig_counts, display_numbers = T,
                 main = goname, filename = fname,
                 fontsize = 5, number_color = "black",  color = colorRampPalette(rev(brewer.pal(n = 7, name = "YlGnBu")))(100))
        write.csv(sig_counts, csv_name)
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
        dev.off()
      }
    }
  }
}

CreateSetOfPlots(normalized_mean, auxin, significant_limma, genes_GO_limma, savedir = "./results/GO/heatmaps/auksyna/")
CreateSetOfPlots(normalized_mean, root, significant_limma, genes_GO_limma, savedir = "./results/GO/heatmaps/korzen/")
CreateSetOfPlots(normalized_mean, other, significant_limma, genes_GO_limma, savedir = "./results/GO/heatmaps/pozostale/")
CreateSetOfPlots(normalized_mean, nutrients, significant_limma, genes_GO_limma, savedir = "./results/GO/heatmaps/nutrienty/")
CreateSetOfPlots(normalized_mean, metals, significant_limma, genes_GO_limma, savedir = "./results/GO/heatmaps/metale/")
CreateSetOfPlots(normalized_mean, ethylen, significant_limma, genes_GO_limma, savedir = "./results/GO/heatmaps/etylen/")
CreateSetOfPlots(normalized_mean, wszystkie, significant_limma, genes_GO_limma, savedir = "./results/GO/heatmaps/wszystkie/")

