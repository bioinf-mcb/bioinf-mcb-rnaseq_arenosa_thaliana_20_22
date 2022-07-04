# Sourcing GO script ----


setwd("D:/MCB/PIOTR_R/LEGACY/Arenosa_proper/")
source("Arabidopsis_arenosa_final_GO.R")
# Importing KEGG libraries ----


library(KEGGprofile)
# KEGG testing ----
genes_shoot <- rownames(significant_shoot_limma)
genes_root <- rownames(significant_root_limma)
shoot_KEGG <- find_enriched_pathway(genes_shoot, species = "ath", download_latest = T)
root_KEGG <- find_enriched_pathway(genes_root, species = "ath", download_latest = T)

write.csv(shoot_KEGG$stastic, "D:/MCB/PIOTR_R/LEGACY/Arenosa_proper/results/KEGG/KEGG_shoot_statistic.csv")
write.csv(root_KEGG$stastic, "D:/MCB/PIOTR_R/LEGACY/Arenosa_proper/results/KEGG/KEGG_root_statistic.csv")

shoot_expr <- normalized_shoot[genes_shoot,]
rownames(shoot_expr)
root_expr <- normalized_root[genes_root,]
rownames(root_expr)


colnames(shoot_expr) <- c("E.1", "E.2", "E.3", "E.4", "E.5",
                          "TM.1", "TM.2", "TM.3", "TM.4", "TM.5")
colnames(root_expr) <- c("E.1", "E.2", "E.3", "E.4", "E.5",
                         "TM.1", "TM.2", "TM.3", "TM.4")
root_expr <- root_expr[,c(1,2,3,4,6,7,8,9)]

# KEGG correlation plotting ----
plot_pathway_cor(gene_expr = shoot_expr, kegg_enriched_pathway = shoot_KEGG)
plot_pathway_cor(gene_expr = root_expr, kegg_enriched_pathway = root_KEGG)

# KEGG root pathway map plotting ----

root_expr <- normalized_root[genes_root,]

E_expr <- root_expr[,1:5]
temp <- apply(E_expr, 1, function(x) length(which(is.na(x))))
E_expr <- E_expr[which(temp == 0), ]
col <- col_by_value(E_expr, col = colorRampPalette(c("blue", "violet", "red"))(1024),
                    range = c(-6, 6))


setwd("D:/MCB/PIOTR_R/LEGACY/Arenosa_proper")
KEGGdf <- root_KEGG$stastic
for(i in 1:dim(root_KEGG$stastic)[1]){
  out <- tryCatch(
    {
      KEGG_n <- rownames(KEGGdf)[i]
      KEGG_id <- KEGGdf[i,1]
      download_KEGGfile(pathway_id = KEGG_n, species = "ath", target_dir = getwd())
      KEGG_database <- parse_XMLfile(KEGG_n, species = "ath", database_dir = getwd())
      plot_profile(gene_expr = E_expr, pathway_name = paste0("ath",KEGG_n), result_name = paste0("results/KEGG/root/",KEGG_id,".png"),
                   KEGG_database = KEGG_database, text_col = "white", type = "bg", bg_col = col,
                   species = "ath", magnify = 1.2, pathway_min = 1)
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
        #dev.off()
    }
  }
}

# KEGG shoot pathway map plotting ----
E_expr <- shoot_expr[,1:5]
temp <- apply(E_expr, 1, function(x) length(which(is.na(x))))
E_expr <- E_expr[which(temp == 0), ]
col <- col_by_value(E_expr, col = colorRampPalette(c("blue", "violet", "red"))(1024),
                    range = c(-6, 6))


KEGGdf <- shoot_KEGG$stastic
for(i in 1:dim(shoot_KEGG$stastic)[1]){
  out <- tryCatch(
    {
      KEGG_n <- rownames(KEGGdf)[i]
      KEGG_id <- KEGGdf[i,1]
      download_KEGGfile(pathway_id = KEGG_n, species = "ath", target_dir = getwd())
      KEGG_database <- parse_XMLfile(KEGG_n, species = "ath", database_dir = getwd())
      plot_profile(gene_expr = E_expr, pathway_name = paste0("ath",KEGG_n), result_name = paste0("results/KEGG/shoot/",KEGG_id,".png"),
                   KEGG_database = KEGG_database, text_col = "white", type = "bg", bg_col = col,
                   species = "ath", magnify = 1.2, pathway_min = 1)
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
      # dev.off()
    }
  }
}
