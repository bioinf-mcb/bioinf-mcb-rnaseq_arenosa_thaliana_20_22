# Sourcing GO script ----
setwd("D:/MCB/PIOTR_R/LEGACy/Thaliana_proper/")
source("Arabidopsis_thaliana_final_GO.R")

# Importing KEGG libraries ----
library(KEGGprofile)

# KEGG testing ----
ni_genes = rownames(significant_limma)
ni_KEGG <- find_enriched_pathway(ni_genes, species = "ath", returned_genenumber = 1, download_latest = F, returned_pvalue = 0.05)
write.csv(ni_KEGG$stastic, "./results/KEGG/KEGG_thaliana_statistic.csv")

expr <- normalized_counts
rownames(expr) <- rownames(expr)
colnames(expr) <- c("E.1", "E.2" ,"E.3", "E.4", "E.5",
                    "R.1", "R.2", "R.3", "R.4")

# KEGG correlation plotting ----
plot_pathway_cor(gene_expr = expr[,c(1,2,3,4,6,7,8,9)], kegg_enriched_pathway = ni_KEGG)


# KEGG root pathway map plotting ----
E_expr <- expr[,1:5]
temp <- apply(E_expr, 1, function(x) length(which(is.na(x))))
E_expr <- E_expr[which(temp == 0), ]
col <- col_by_value(E_expr, col = colorRampPalette(c("blue", "violet", "red"))(1024),
                    range = c(-6, 6))
col

E_expr["AT5G03280",]

KEGGdf <- ni_KEGG$stastic

for(i in 1){
  out <- tryCatch(
    {
      KEGG_n <- "04016"
      KEGG_id <- "MAPK"
      download_KEGGfile(pathway_id = KEGG_n, species = "ath", target_dir = getwd())
      KEGG_database <- parse_XMLfile(KEGG_n, species = "ath", database_dir = getwd())
      print(KEGG_database)
      plot_profile(gene_expr = E_expr, pathway_name = paste0("ath",KEGG_n), result_name = paste0("results/KEGG/",KEGG_id,".png"),
                   KEGG_database = KEGG_database, text_col = "white", type = "bg", bg_col = col,
                   species = "ath", magnify = 1.2, pathway_min = 0)
    }, error = function(cond){
      print("Problem with XML or too global perspective")
    }, warning = function(cond){
      print("Problem with XML or too global perspective")
    }, finally = function(cond){
      print("Problem with XML or too global perspective")
    }
  )
  if(length(dev.list()) > 0){
    for (i in dev.list()[1]:dev.list()[length(dev.list())]) {
      dev.off()
    }
  }
}


for(i in 1:dim(ni_KEGG$stastic)[1]){
  out <- tryCatch(
    {
      KEGG_n <- rownames(KEGGdf)[i]
      KEGG_id <- KEGGdf[i,1]
      download_KEGGfile(pathway_id = KEGG_n, species = "ath", target_dir = getwd())
      KEGG_database <- parse_XMLfile(KEGG_n, species = "ath", database_dir = getwd())
      plot_profile(gene_expr = E_expr, pathway_name = paste0("ath",KEGG_n), result_name = paste0("results/KEGG/",KEGG_id,".png"),
                   KEGG_database = KEGG_database, text_col = "white", type = "bg", bg_col = col,
                   species = "ath", magnify = 1.2, pathway_min = 0)
    }, error = function(cond){
      print("Problem with XML or too global perspective")
    }, warning = function(cond){
      print("Problem with XML or too global perspective")
    }, finally = function(cond){
      print("Problem with XML or too global perspective")
    }
  )
}

