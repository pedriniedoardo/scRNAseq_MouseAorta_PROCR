# AIM ---------------------------------------------------------------------
# the aim is to plot the data after integration.
# this integration was run by skipping the seurat integration (by merging the matrices) and running Harmony

# libraries ---------------------------------------------------------------
library(harmony)
library(Seurat)
library(dplyr)
library(cowplot)
library(tidyverse)
library(ggrepel)
library(scales)
library(RColorBrewer)
library(SeuratWrappers)
library(dittoSeq)
library(clustree)
library(pals)
library(patchwork)
library(magick)
library(homologene)

# read in the data --------------------------------------------------------
data.combined <- readRDS("../../out/object/data.combined_harmonySkipIntegration.rds")
DimPlot(data.combined)
# Idents(data.combined)<-"seurat_clusters"

# read in the LUT for homologenes
# library(homologene)
# read in the updated reference table
homologeneData2_240528 <- read_tsv("../../data/Homologene_240528") %>%
  as.data.frame()

# plots -------------------------------------------------------------------
# plot the tree of the cluster dependencies. this will justify the choice of the resolution, not too granula for the moment.
# library(clustree)
# define the columns to use
id_columns <- str_subset(colnames(data.combined@meta.data),pattern = "RNA_snn_res") %>%
  str_subset(negate = T,pattern = ".paper")

clustree::clustree(data.combined@meta.data[,id_columns],
                   prefix = "RNA_snn_res.")
ggsave("../../out/plot/02_UMAPCluster_tree.pdf",width = 10,height = 10)

# plot the UMAP with all the resolutions runs
id_resolution <- str_subset(colnames(data.combined@meta.data),pattern = "RNA_snn_res") %>%
  str_subset(negate = T,pattern = ".paper") %>%
  sort()

list_plot <- lapply(id_resolution,function(x){
  plot <- DimPlot(data.combined,
                  reduction = "umap",
                  group.by = x,
                  label = T,
                  raster = F)
  return(plot)
})

wrap_plots(list_plot)
ggsave("../../out/plot/02_UMAPCluster_resolutions.pdf",width = 20,height = 15)

# save the panel martina uses for the annotation
shortlist_features_list_long_mouse <- list(
  ENDO = c("Pecam1","Cdh5","Clu","Cytl1"),
  FIBRO = c("Igfbp5","Dcn","Col1a2","Gsn"),
  RBc = c("Hbb-bt","Hbb-bs","Hba-a1","Hba-a2"),
  MAC = c("Pld4","C1qa","C1qb"),
  DC = c("Pld4","C1qa","C1qb","Tnip3"),
  Bc = c("Cd79a","Cd79b","Ms4a1","Ly6d"),
  Tc = c("Cd3g","Ms4a4b","Cd3d","Nkg7")
  
)

# if needed conver the panel of genes
# shortlist_features_list_long_mouse <- lapply(shortlist_features_list_long,function(x){
#   df <- homologene(x, inTax = 9606, outTax = 10090,db = homologeneData2_240528)
#   df %>%
#     pull(`10090`)
# })

# make the genes unique to avoid duplication issue
GOI <- unique(unlist(shortlist_features_list_long_mouse))
GOI

# define the grouping variable
Idents(data.combined) <- "RNA_snn_res.0.4"

test <- DotPlot(data.combined, features = GOI,cluster.idents = T)

df_test <- lapply(shortlist_features_list_long_mouse,function(x){
  test$data %>% 
    filter(features.plot %in% x)
}) %>% 
  bind_rows(.id = "cell_type")

head(df_test)

# plot the shortlisted feature per cluster
# notice that this is done only on the subset of the young (control) cells
test_long01 <- df_test %>%
  # force the order
  mutate(id = factor(id,levels = c(9,6,5,3,2,8,4,7,0,1))) %>% 
  mutate(cell_type = factor(cell_type,levels = c("Bc","Tc","RBc","FIBRO","MAC","DC","ENDO"))) %>% 
  ggplot(aes(x = features.plot,y = id)) +
  geom_point(aes(size = pct.exp, col = avg.exp.scaled))+
  scale_radius(range = c(0, 8)) +
  facet_grid(~cell_type,scales = "free",space = "free")+
  theme_cowplot()+
  theme(strip.background = element_blank(),
        axis.text.x = element_text(hjust = 1,angle = 90),
        strip.text.x = element_text(angle = 0))+
  scale_color_gradient(low = "lightgrey",high = "blue",limits = c(-1,2),oob = scales::squish)
ggsave(plot=test_long01,"../../out/plot/02_DotplotLong_res0.4.pdf",width = 10,height = 5)

# try to convey the the classification using also the UMAP and the score module
# seu.int <- AddModuleScore(seu.int, features = sig.list, name = "_score")
# names(seu.int@meta.data)[grep("_score", names(seu.int@meta.data))] <- names(sig.list)
# FeaturePlot(seu.int, features = names(sig.list))

# Rename using a named vector and `all_of()`
data.combined <- Seurat::AddModuleScore(data.combined,
                                        features = shortlist_features_list_long_mouse,
                                        name = "_score")

df_rename <- data.frame(names = data.combined@meta.data %>%
                          colnames() %>%
                          str_subset("_score"),
                        rename = paste0("score_",names(shortlist_features_list_long_mouse)))

lookup <- df_rename$names
names(lookup) <- df_rename$rename

# rename the columns
data.combined@meta.data <- dplyr::rename(data.combined@meta.data,all_of(lookup))

# plot the scores from AddModuleScore
list_plot_02 <- lapply(df_rename$rename,function(x){
  plot <- FeaturePlot(data.combined,features = x,order = T,
                      reduction = "umap",
                      raster = F) + scale_color_viridis_c(option = "inferno")
  return(plot)
})

wrap_plots(list_plot_02)
ggsave("../../out/plot/02_UMAPCluster_ExpertAnnotaiton.pdf",width = 15,height = 15)
# 
# # same as above but as violin plot
# list_plot <- lapply(df_rename$rename, function(x){ 
#   test <- VlnPlot(object = data.combined,features = x, group.by = "RNA_snn_res.0.1",raster = T)
#   return(test)
# })
# 
# # make it a dataframe
# # x <- list_plot[[1]]
# df_violin <- lapply(list_plot,function(x){ 
#   df <- x[[1]]$data 
#   
#   # extract the name of the gene 
#   feature <- colnames(df)[1] 
#   
#   df %>% 
#     mutate(feature = feature) %>% 
#     setNames(c("value","ident","feature")) 
# }) %>% 
#   bind_rows()
# 
# head(df_violin) 
# 
# # how many cells per ident
# df_violin %>%
#   group_by(ident,feature) %>%
#   summarise(n = n()) %>%
#   filter(feature %in% c("score_ASTRO"))
# 
# # plot at maximum 1000 cells per group
# set.seed(123)
# df_plot_violin <- df_violin %>% 
#   group_by(ident,feature) %>%
#   sample_n(size = 1000,replace = F) %>%
#   ungroup()
# 
# df_plot_violin_summary <- df_plot_violin %>%
#   group_by(feature) %>%
#   summarise(med_score = median(value))
# 
# df_plot_violin %>%
#   ggplot(aes(x = ident, y = value)) + 
#   geom_violin(scale = "width")+ 
#   #geom_boxplot(outlier.shape = NA,position = position_dodge(width=0.9),width=0.05) + 
#   geom_point(position = position_jitter(width = 0.2),alpha = 0.05,size = 0.5) + 
#   facet_wrap(~feature,ncol = 1,scales = "free") + 
#   theme_bw() + 
#   geom_hline(data = df_plot_violin_summary,aes(yintercept = med_score),linetype="dashed",col="red") +
#   theme(strip.background = element_blank(),
#         axis.text.x = element_text(hjust = 1,angle = 45))
# ggsave("../../out/plot/ViolinCluster_ExpertAnnotaiton.pdf",width = 7,height = 21)

# use a similar approach to score potential technical clusters
# plot the scores from AddModuleScore
list_plot_technical <- lapply(c("nFeature_RNA","percent.mt","percent.ribo","percent.globin"),function(x){
  plot <- FeaturePlot(data.combined,features = x,order = T,
                      reduction = "umap",
                      raster = F) + scale_color_viridis_c(option = "inferno")
  return(plot)
})

wrap_plots(list_plot_technical)
ggsave("../../out/plot/02_UMAPCluster_technical.pdf",width = 10,height = 8)

list_plot_technical <- lapply(c("nFeature_RNA","percent.mt","percent.ribo","percent.globin"), function(x){ 
  test <- VlnPlot(object = data.combined,features = x, group.by = "RNA_snn_res.0.1",raster = T)
  return(test)
})

# make it a dataframe
# x <- list_plot[[1]]
df_violin_technical <- lapply(list_plot_technical,function(x){ 
  df <- x[[1]]$data 
  
  # extract the name of the gene 
  feature <- colnames(df)[1] 
  
  df %>% 
    mutate(feature = feature) %>% 
    setNames(c("value","ident","feature")) 
}) %>% 
  bind_rows()

head(df_violin_technical) 

# plot at maximum 1000 cells per group
set.seed(123)
df_plot_violin_technical <- df_violin_technical %>% 
  group_by(ident,feature) %>%
  # sample_n(size = 300,replace = F) %>%
  ungroup()

df_plot_violin_technical_summary <- df_plot_violin_technical %>%
  group_by(feature) %>%
  summarise(med_score = median(value))

df_plot_violin_technical %>%
  ggplot(aes(x = ident, y = value)) + 
  geom_violin(scale = "width")+ 
  #geom_boxplot(outlier.shape = NA,position = position_dodge(width=0.9),width=0.05) + 
  geom_point(position = position_jitter(width = 0.2),alpha = 0.05,size = 0.5) + 
  facet_wrap(~feature,ncol = 1,scales = "free") + 
  theme_bw() + 
  geom_hline(data = df_plot_violin_technical_summary,aes(yintercept = med_score),linetype="dashed",col="red") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/plot/02_ViolinCluster_technical_res0.4.pdf",width = 7,height = 12)

# pick one resolution onward 0.4

# main umap
plot03 <- DimPlot(data.combined, reduction = "umap", group.by = "RNA_snn_res.0.4",label = T,raster = F)
ggsave(plot=plot03,"../../out/plot/02_UMAPCluster_data.combined_harmonySkipIntegration_res0.4.pdf",width = 6,height = 5)

plot03a <- DimPlot(data.combined, reduction = "umap", group.by = "RNA_snn_res.0.4",label = T,raster = F,split.by = "orig.ident",ncol = 5)
ggsave(plot=plot03a,"../../out/plot/02_UMAPClusterSplit_data.combined_harmonySkipIntegration_res0.4.pdf",width = 8,height = 4)

# main umap cell cycle
plot04 <- DimPlot(data.combined, reduction = "umap", group.by = "Phase",raster = F,order = T)
ggsave(plot=plot04,"../../out/plot/02_UMAPPhase_data.combined_harmonySkipIntegration.pdf",width = 5,height = 4)

# split by sample
Idents(data.combined) <- "RNA_snn_res.0.4"
df_meta <- data.combined@meta.data %>%
  rownames_to_column("rowname")
df_UMAP <- data.combined@reductions$umap@cell.embeddings %>%
  data.frame() %>%
  rownames_to_column("rowname")

# plot the proporition for the phase per cluster
df_summary_phase <- df_meta %>%
  group_by(RNA_snn_res.0.4,Phase) %>%
  summarise(n = n()) %>%
  group_by(RNA_snn_res.0.4) %>%
  mutate(tot=sum(n)) %>%
  ungroup() %>%
  mutate(prop=n/tot)

df_summary_phase %>%
  mutate(RNA_snn_res.0.4 = fct_relevel(RNA_snn_res.0.4,as.character(0:9))) %>%
  ggplot() +
  geom_col(aes(x=RNA_snn_res.0.4,y=prop,fill=Phase))+theme_cowplot()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))
ggsave("../../out/plot/02_BarplotPhase_summary_data.combined_harmonySkipIntegration_res0.4.pdf",width = 7,height = 6)

# proportion of cell per cluster
df_summary <- df_meta %>%
  group_by(orig.ident,RNA_snn_res.0.4) %>%
  summarise(n = n()) %>%
  group_by(orig.ident) %>%
  mutate(tot=sum(n)) %>%
  ungroup() %>%
  mutate(prop=n/tot)
write_tsv(df_summary,"../../out/table/02_summary_data.combined_harmonySkipIntegration_res0.4.tsv")

color_id <- alphabet(length(unique(df_summary$RNA_snn_res.0.4)))
# check the colors
show_col(color_id)

df_summary %>%
  ggplot() +
  geom_col(aes(x=orig.ident,y=prop,fill=RNA_snn_res.0.4))+theme_cowplot()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90)) +
  scale_fill_manual(values = unname(color_id))
ggsave("../../out/plot/02_Barplot_summary_data.combined_harmonySkipIntegration_res0.4.pdf",width = 7,height = 6)

# render the same plot as an heatmap
sample_prop_wide <- df_summary %>%
  # scale by rows
  group_by(RNA_snn_res.0.4) %>%
  mutate(zscore = (prop-mean(prop))/sd(prop)) %>%
  # make it long
  dplyr::select(orig.ident,RNA_snn_res.0.4,zscore) %>%
  pivot_wider(names_from = orig.ident,values_from = zscore,values_fill = 0) %>%
  # dplyr::select(orig.ident,RNA_snn_res.0.4,prop) %>%
  # pivot_wider(names_from = orig.ident,values_from = prop,values_fill = 0) %>%
  column_to_rownames("RNA_snn_res.0.4")

rowSums(sample_prop_wide)
colSums(sample_prop_wide)

# plot the data as heatmap
# meta_sample_prop <- data.frame(colname = colnames(sample_prop_wide)) %>%
#   left_join(df_meta %>%
#               group_by(original_sample_name,sex,diagnosis,location) %>%
#               summarise(),by=c("colname" = "original_sample_name"))

# column_meta_sample_prop <- HeatmapAnnotation(gender = meta_sample_prop$sex,
#                                              location = meta_sample_prop$location,
#                                              diagnosis = meta_sample_prop$diagnosis,
#                                              col = list(gender = c("m" = "blue",
#                                                                    "f" = "pink"),
#                                                         location = c("lumbar" = "gray90",
#                                                                      "thoracic" = "gray50",
#                                                                      "cervical" = "black"),
#                                                         diagnosis = c("Non-demented control" = "green",
#                                                                       "Multiple sclerosis" = "red")))

ht2_shr_MG2 <- Heatmap(sample_prop_wide, show_column_names = T,raster_by_magick = T,show_row_dend = F, use_raster = T,
                       name = "prop_cell_type \nscale cluster",
                       column_title = "sample",
                       col = viridis::viridis(option = "turbo",n = 10),
                       
                       # row_names_gp = gpar(fontsize = 3),
                       # top_annotation = column_meta_sample_prop,
                       show_row_names = T
                       # cluster_rows = F,
                       # right_annotation = row_ha,
                       # row_split = rep(c(1,2,3,4),c(2,3,4,7))
                       
)
# pdf("../../out/plot/HeatmapCluster_summary_data.combined_harmonySkipIntegration_manualClean_res0.1.pdf",width = 3,height = 6)
draw(ht2_shr_MG2,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(30,2,2, 2), "mm"))
# dev.off()

# render the same plot as an heatmap
sample_prop_wide2 <- df_summary %>%
  # scale by rows
  group_by(orig.ident) %>%
  mutate(zscore = (prop-mean(prop))/sd(prop)) %>%
  # make it long
  dplyr::select(orig.ident,RNA_snn_res.0.4,zscore) %>%
  pivot_wider(names_from = orig.ident,values_from = zscore,values_fill = 0) %>%
  column_to_rownames("RNA_snn_res.0.4")

rowSums(sample_prop_wide2)
colSums(sample_prop_wide2)

# # plot the data as heatmap
# meta_sample_prop2 <- data.frame(colname = colnames(sample_prop_wide2)) %>%
#   left_join(df_meta %>%
#               group_by(original_sample_name,sex,diagnosis,location) %>%
#               summarise(),by=c("colname" = "original_sample_name"))
# 
# column_meta_sample_prop2 <- HeatmapAnnotation(gender = meta_sample_prop2$sex,
#                                               location = meta_sample_prop2$location,
#                                               diagnosis = meta_sample_prop2$diagnosis,
#                                               col = list(gender = c("m" = "blue",
#                                                                     "f" = "pink"),
#                                                          location = c("lumbar" = "gray90",
#                                                                       "thoracic" = "gray50",
#                                                                       "cervical" = "black"),
#                                                          diagnosis = c("Non-demented control" = "green",
#                                                                        "Multiple sclerosis" = "red")))

ht2_shr_MG22 <- Heatmap(sample_prop_wide2, show_column_names = T,raster_by_magick = T,show_row_dend = F, use_raster = T,
                        name = "zscore \nprop_cell_type \nscale sample",
                        column_title = "sample",
                        col = viridis::viridis(option = "turbo",n = 10),
                        
                        # row_names_gp = gpar(fontsize = 3),
                        # top_annotation = column_meta_sample_prop2,
                        show_row_names = T
                        # cluster_rows = F,
                        # right_annotation = row_ha,
                        # row_split = rep(c(1,2,3,4),c(2,3,4,7))
                        
)
# pdf("../../out/plot/HeatmapSample_summary_data.combined_harmonySkipIntegration_manualClean_res0.1.pdf",width = 3,height = 6)
draw(ht2_shr_MG22,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(30,2,2, 2), "mm"))
# dev.off()

# Identify conserved cell type markers ------------------------------------
# data
DefaultAssay(data.combined) <- "RNA"
# notice that in this case the data object of the RNA slot is already filled with the normalzied data, therefore in this case (following the Normalize workfloe for the integration) there is no need to run the NormalizeData on the RNA slot of the integrated object
# sobj_total_h@assays$RNA@data[20:50,1:10]
# dim(sobj_total_h@assays$RNA@data)
# 
# # scale the data see the note on evernote on why this can be done also at this point. this is needed because the scale.data is empty
# sobj_total_h@assays$RNA@scale.data
# all.genes <- rownames(sobj_total_h)
# sobj_total_h <- ScaleData(sobj_total_h,vars.to.regress = c("percent.mt","nCount_RNA","S.Score","G2M.Score"))
# # confirm now the slot is loaded
# sobj_total_h@assays$RNA@scale.data[1:10,1:10]

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
Idents(data.combined) <- "RNA_snn_res.0.4"
sobj_total_h.markers <- RunPrestoAll(data.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# save the table of all markers
sobj_total_h.markers %>%
  write_tsv("../../out/table/02_FindAllMarkers_data.combined_harmonySkipIntegration_res0.4.tsv")

# # save the top 100
# sobj_total_h.markers %>%
#   group_by(cluster) %>%
#   dplyr::slice(1:100) %>%
#   write_tsv("../../out/table/FindAllMarkers_harmonySkipIntegration_AllSoupX_00500_07000_05_top100.tsv")

# sobj_total_h.markers <- read_tsv("../../out/table/FindAllMarkers_harmonySkipIntegration_AllSoupX_01000_06000_15.tsv")

# try plotting the top markers
top_specific_markers <- sobj_total_h.markers %>%
  # filter ribosomal and mt genes
  filter(str_detect(gene,pattern = "^mt-",negate=T)) %>%
  filter(str_detect(gene,pattern = "^Rp[sl][[:digit:]]|^Rplp[[:digit:]]|^Rpsa",negate=T)) %>%
  filter(str_detect(gene,pattern = "^Hb[^(p)]",negate=T)) %>%
  group_by(cluster) %>%
  top_n(5, avg_log2FC)

# And generate e.g. a dotplot:
dittoSeq::dittoDotPlot(data.combined,
                       vars = unique(top_specific_markers$gene), 
                       group.by = "RNA_snn_res.0.4")+scale_color_viridis_c(option = "turbo",name="relative \nexpression")
ggsave("../../out/plot/02_TopMarkersDitto_data.combined_harmonySkipIntegration_res0.4.pdf",width = 12,height = 4)

# attempt annotation based on the results ---------------------------------
data.combined@meta.data$expertAnno.l1 <- 
  data.combined@meta.data %>%
  mutate(expertAnno.l1 = case_when(RNA_snn_res.0.4 %in% c(9) ~ "Bc",
                                   RNA_snn_res.0.4 %in% c(6) ~ "Tc",
                                   RNA_snn_res.0.4 %in% c(5) ~ "RBc",
                                   RNA_snn_res.0.4 %in% c(2,3) ~ "FIBRO",
                                   RNA_snn_res.0.4 %in% c(8) ~ "MAC",
                                   RNA_snn_res.0.4 %in% c(4) ~ "DC",
                                   RNA_snn_res.0.4 %in% c(0,1,7) ~ "ENDO",
                                   TRUE ~ "NA")) %>%
  pull(expertAnno.l1)

# main umap
plot032 <- DimPlot(data.combined, reduction = "umap", group.by = "expertAnno.l1",label = T,raster = F)
ggsave(plot=plot032,"../../out/plot/02_UMAPCluster_data.combined_harmonySkipIntegration_cellid.pdf",width = 6,height = 5)

plot032a <- DimPlot(data.combined, reduction = "umap", group.by = "expertAnno.l1",label = T,raster = F,split.by = "orig.ident",ncol = 5)
ggsave(plot=plot032a,"../../out/plot/02_UMAPClusterSplit_data.combined_harmonySkipIntegration_cellid.pdf",width = 8,height = 4)

# split by sample
Idents(data.combined) <- "expertAnno.l1"
df_meta2 <- data.combined@meta.data %>%
  rownames_to_column("rowname")
df_UMAP2 <- data.combined@reductions$umap@cell.embeddings %>%
  data.frame() %>%
  rownames_to_column("rowname")

# plot the proporition for the phase per cluster
df_summary_phase2 <- df_meta2 %>%
  group_by(expertAnno.l1,Phase) %>%
  summarise(n = n()) %>%
  group_by(expertAnno.l1) %>%
  mutate(tot=sum(n)) %>%
  ungroup() %>%
  mutate(prop=n/tot)

df_summary_phase2 %>%
  mutate(expertAnno.l1 = fct_relevel(expertAnno.l1)) %>%
  ggplot() +
  geom_col(aes(x=expertAnno.l1,y=prop,fill=Phase))+theme_cowplot()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))
ggsave("../../out/plot/02_BarplotPhase_summary_data.combined_harmonySkipIntegration_cellid.pdf",width = 7,height = 6)

# proportion of cell per cluster
df_summary2 <- df_meta2 %>%
  group_by(orig.ident,expertAnno.l1) %>%
  summarise(n = n()) %>%
  group_by(orig.ident) %>%
  mutate(tot=sum(n)) %>%
  ungroup() %>%
  mutate(prop=n/tot)
write_tsv(df_summary2,"../../out/table/02_summary_data.combined_harmonySkipIntegration_cellid.tsv")

color_id2 <- alphabet(length(unique(df_summary2$expertAnno.l1)))
# check the colors
show_col(color_id2)

df_summary2 %>%
  ggplot() +
  geom_col(aes(x=orig.ident,y=prop,fill=expertAnno.l1))+theme_cowplot()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90)) +
  scale_fill_manual(values = unname(color_id))
ggsave("../../out/plot/02_Barplot_summary_data.combined_harmonySkipIntegration_cellid.pdf",width = 7,height = 6)

# calculate the ratio of the proportions per sample
df_summary3 <- df_meta2 %>%
  group_by(orig.ident,expertAnno.l1) %>%
  summarise(n = n()) %>%
  group_by(orig.ident) %>%
  mutate(tot=sum(n)) %>%
  ungroup() %>%
  mutate(nCell_100k=n/tot*100000) %>%
  ungroup() %>%
  # pivot_wider(names_from = orig.ident,values_from = nCell_100k,id_cols = -c(n,tot),values_fill = 0) %>%
  group_by(expertAnno.l1) %>%
  mutate(sumCell_100k = sum(nCell_100k)) %>%
  mutate(ratioCell_100k = nCell_100k/sumCell_100k) %>%
  ungroup()
write_tsv(df_summary3,"../../out/table/02_summary_data.combined_harmonySkipIntegration_cellid2.tsv")

color_id3 <- c("black","gray")
# check the colors
show_col(color_id3)

df_summary3 %>%
  ggplot() +
  geom_col(aes(x=expertAnno.l1,y=ratioCell_100k,fill=orig.ident))+theme_cowplot()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90)) +
  scale_fill_manual(values = unname(color_id3))
ggsave("../../out/plot/02_Barplot_summary_data.combined_harmonySkipIntegration_cellid3.pdf",width = 7,height = 6)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
Idents(data.combined) <- "expertAnno.l1"
sobj_total_h.markers2 <- RunPrestoAll(data.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# save the table of all markers
sobj_total_h.markers2 %>%
  write_tsv("../../out/table/02_FindAllMarkers_data.combined_harmonySkipIntegration_cellid.tsv")

# save the object ---------------------------------------------------------
# save the object with the full annotation
saveRDS(data.combined,"../../out/object/02_data.combined_harmonySkipIntegration_Annotation.rds")
