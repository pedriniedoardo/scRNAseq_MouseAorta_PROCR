# library -----------------------------------------------------------------
library(tidyverse)
library(fgsea)
library(msigdbr)
library(GSEABase)
library(patchwork)
library(ComplexHeatmap)

# prepare the dataset with all the annoration needed ----------------------
# locate the files
# file <- dir("out/") %>%
#   str_subset(pattern = "res_") %>%
#   str_subset(pattern = ".txt") %>%
#   str_subset(pattern = "shr",negate = T)
# file <- dir("../../out/table/") %>%
#   str_subset(pattern = "res_") %>%
#   str_subset(pattern = ".txt") %>%
#   str_subset(pattern = "HUVEC") %>%
#   str_subset(pattern = "shr",negate = F)
# file

# load the results
# results <- lapply(paste0("../../out/table/",file),function(x){
#   read_tsv(x)
# }) %>%
#   setNames(str_remove_all(file,pattern = ".txt"))
results <- read_tsv("../../out/table/DE_pseudobulk_ENDO_shr_filterExp.tsv") %>%
  split(f = .$conditionVsCtrl)

# GSEA --------------------------------------------------------------------
# use the FC dataset to create the ranked list of genes
# Symbol or Entrez?
# x <- results$res_MutvsWT_shr
list_ranks <- lapply(results, function(x){
  
  x <- dplyr::filter(x,!is.na(symbol)) %>%
    # average logFC in case of duplicated genenames
    group_by(symbol) %>%
    summarise(logFC = mean(log2FoldChange))
  
  ranks <- setNames(x$logFC, x$symbol)
  ranks
})
glimpse(list_ranks)

# score all the signatures in MsigDB from C2 category ---------------------
# library("msigdbr")
msigdbr_collections() %>%
  print(n=30)
#
reactome_gene_sets <- msigdbr(species = "Mus musculus", category = "C5",subcategory = "GO:BP")
# filter for the signautres of interest
reactome_gene_sets_subset <- reactome_gene_sets %>%
  filter(str_detect(gs_name,"POSITIVE_REGULATION_OF_VASCULAR_ENDOTHELIAL_CELL_PROLIFERATION"))
# reactome_gene_sets <- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "CP:KEGG")
# get all the C2 terms
# cgp_gene_sets <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CGP")
# C2_gene_sets = msigdbr(species = "Mus musculus", category = "C2")
# h_gene_sets <- msigdbr(species = "Mus musculus", category = "H")
head(reactome_gene_sets)

# format in order to be accepted by GSEA
pathways <- split(x = reactome_gene_sets_subset$gene_symbol, f = reactome_gene_sets_subset$gs_name)
# head(pathways)

# RUN GSEA ----------------------------------------------------------------
df_tables_GSEA_all <- lapply(list_ranks, function(x){
  fgsea(pathways, x, minSize=10, maxSize=500)
}) %>%
  bind_rows(.id = "dataset") %>%
  # the ladingEdge columns has to be re-arranged in order to save the file as a table (originally is a list)
  mutate(leadingEdge = unlist(lapply(.$leadingEdge, function(x){
    paste0(x,collapse = "|")
  }))) %>%
  arrange(dataset,padj,-abs(NES))

dim(df_tables_GSEA_all)
head(df_tables_GSEA_all,n=20)

# save the whole table
df_tables_GSEA_all %>%
  write_tsv("../../out/table/GSEA_Custom_DE_pseudobulk_ENDO_shr.tsv")


# library(ggrepel)
df_tables_GSEA_all %>%
  # shorten the label of the pathway
  mutate(pathway2 = str_remove(pathway,pattern = "REACTOME_") %>%
           str_sub(start = 1,end = 35)) %>%
  # mutate(min_log10_padj = -log10(padj)) %>%
  ggplot(aes(y = -log10(padj),x = NES,label = pathway2)) + geom_point(aes(size = size),alpha = 0.2) + facet_wrap(~dataset) + theme_bw() +
  theme(strip.background = element_blank())+
  geom_hline(yintercept = -log10(0.05),linetype="dashed",col="gray")+
  geom_text_repel(size = 2,box.padding = 0.5,segment.alpha = 0.6,max.overlaps = 10)
ggsave("../../out/plot/GSEA_REACTOME_DE_pseudobulk_ENDO_shr.pdf",width = 15,height = 10)

# PLOT PROFILE ------------------------------------------------------------
# library("patchwork")
list_mainPathways <- names(pathways)

# plot a single pathway across conditions
plot_pathway <- pmap(list(list_ranks,
                          names(list_ranks)),function(x,y){
  # name <- str_sub("POS_REG_OF_VASC_ENDO_CELL_PROL",start = 6,end = -4)
  name <- "POS_REG_OF_VASC_ENDO_CELL_PROL"
  plotEnrichment(pathways$GOBP_POSITIVE_REGULATION_OF_VASCULAR_ENDOTHELIAL_CELL_PROLIFERATION,
                 x) + labs(title = y,subtitle = name)
})

wrap_plots(plot_pathway)
ggsave(filename = "../../out/plot/GSEA_plot_profile_ENDO_shr_nonredundant_custom.pdf",width = 15,height = 5)
# wrap_plots(plot_pathway) + ggsave(filename = "image/GSEA_plot_profile_CGP.pdf",width = 15,height = 10)
# wrap_plots(plot_pathway) + ggsave(filename = "image/GSEA_plot_profile_H.pdf",width = 15,height = 10)

# PLOT RANK ---------------------------------------------------------------
# get the top and the bottom genes of the signatures in the rank 
df_rank <- results$res_ENDO.D3_shr %>% 
  arrange(desc(log2FoldChange)) %>% 
  mutate(rank = nrow(.):1) %>% 
  mutate(rank2 = 1:nrow(.))

df_rank 

# subset the rank for the terms of interest
id_UP <- df_rank %>% 
  filter(symbol %in% pathways$GOBP_POSITIVE_REGULATION_OF_VASCULAR_ENDOTHELIAL_CELL_PROLIFERATION)

id_UP 

# produce the heatmap
m <- matrix(df_rank$rank,ncol = 1) 
ha_up <- rowAnnotation(foo = anno_mark(at = id_UP$rank2, labels = id_UP$symbol)) 
hm_up <- Heatmap(m, name = "mat", cluster_rows = FALSE, right_annotation = ha_up,column_title = "res_ENDO.D3_shr")

draw(hm_up,heatmap_legend_side = "left",annotation_legend_side = "left")

