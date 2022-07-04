library(pheatmap)
library(RColorBrewer)

Arenosa_voom <- significant_root_limma
Thaliana_voom <- read.csv("D:/MCB/PIOTR_R/LEGACY/Thaliana_proper/results/DE/significant_limma.csv")

Arenosa_counts <- v_root$E
Thaliana_counts <- read.csv('D:/MCB/PIOTR_R/LEGACY/Thaliana_proper/counts_thaliana.csv')

rownames(Thaliana_counts) <- Thaliana_counts$X
Thaliana_counts <- Thaliana_counts[,2:dim(Thaliana_counts)[2]]

rownames(Thaliana_voom) <- Thaliana_voom$X
Thaliana_voom <- Thaliana_voom[,2:dim(Thaliana_voom)[2]]

Mean_Endophyte_Arenosa <- as.vector(rowMeans(Arenosa_counts[,1:5]))
Mean_Wild_Arenosa <- as.vector(rowMeans(Arenosa_counts[,6:9]))
normalized_mean_root <- data.frame(Mean_Endophyte_Arenosa, Mean_Wild_Arenosa)
rownames(normalized_mean_root) <- rownames(Arenosa_counts)
Arenosa_counts <- normalized_mean_root

Mean_Endophyte_Thaliana <- as.vector(rowMeans(Thaliana_counts[,1:5]))
Mean_Wild_Thaliana <- as.vector(rowMeans(Thaliana_counts[,6:9]))
normalized_mean_root <- data.frame(Mean_Endophyte_Thaliana, Mean_Wild_Thaliana)
rownames(normalized_mean_root) <- rownames(Thaliana_counts)
Thaliana_counts <- normalized_mean_root


unified_genes <- intersect(rownames(Arenosa_counts), rownames(Thaliana_counts))



counts_merged <- cbind(Arenosa_counts[unified_genes,],
                       Thaliana_counts[unified_genes,])


load("D:/MCB/PIOTR_R/LEGACY/Thaliana_proper/genes_GO_limma.RData")

Thaliana_GO <- genes_GO_limma
Arenosa_GO <- gogenes_root
ETR_GO
WT_GO
arenosa_x_thaliana_intersection <- intersect(names(Arenosa_GO), names(Thaliana_GO))
ETR_x_WT_intersection <- intersect(names(ETR_GO), names(WT_GO))
merged_GO <- taRifx::merge.list(Thaliana_GO, Arenosa_GO)

setwd("D:/MCB/PIOTR_R/LEGACY/Arenosa_proper/results/GO/")


counts_merged <- counts_merged[,c(1,3,2,4)]
colnames(counts_merged) <- c(
  "Mean Expression A.arenosa E+",
  "Mean Expression A.thaliana Ni E+",
  "Mean Expression A.arenosa E-",
  "Mean Expression A.thaliana Ni E-"
)


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

CreateSetOfPlots <- function(normcounts, golist, deresult_arenosa, deresult_thaliana, gogenes, savedir, root = T){
  for(i in 1:length(golist)){
    print(golist[i])
    all_genes <- as.vector(unlist(gogenes[golist[i]]))
    sig_genes_arenosa <- intersect(rownames(deresult_arenosa), all_genes)
    sig_genes_thaliana <- intersect(rownames(deresult_thaliana), all_genes)
    sig_genes <- intersect(sig_genes_arenosa, sig_genes_thaliana)
    print(sig_genes)
    sig_counts <- normcounts[sig_genes,]
    goname <- gsub(":","",golist[i])
    if(root){
      fname <- paste0("./heatmaps/",savedir, "/Arenosa_X_Thaliana/", goname, ".png")
      csv_name <- paste0("./heatmaps/",savedir, "/Arenosa_X_Thaliana/", goname, ".csv")
    }else{
      fname <- paste0("./heatmaps/",savedir, "/shoot/", goname, ".png")
    }
    print(fname)
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

CreateSetOfPlots(counts_merged, auksyna, Arenosa_voom, Thaliana_voom, merged_GO, "auksyna")
CreateSetOfPlots(counts_merged, korzen, Arenosa_voom, Thaliana_voom, merged_GO, "korzen")
CreateSetOfPlots(counts_merged, pozostale, Arenosa_voom, Thaliana_voom, merged_GO, "pozostale")
CreateSetOfPlots(counts_merged, nutrienty, Arenosa_voom, Thaliana_voom, merged_GO, "nutrienty")
CreateSetOfPlots(counts_merged, metale, Arenosa_voom, Thaliana_voom, merged_GO, "metale")
CreateSetOfPlots(counts_merged, etylen, Arenosa_voom, Thaliana_voom, merged_GO, "etylen")
CreateSetOfPlots(counts_merged, arenosa_x_thaliana_intersection, Arenosa_voom, Thaliana_voom, merged_GO, "wspolne")