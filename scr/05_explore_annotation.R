# AIM ---------------------------------------------------------------------
# compare the partial macro annotation with the one suggested by the authors

# libraries ---------------------------------------------------------------
library(Seurat)
library(cowplot)
library(tidyverse)
library(ComplexHeatmap)

# read in the data --------------------------------------------------------
# manual re-annotated dataset
data.combined <- readRDS("../../out/object/data.combined_harmonySkipIntegration_manualClean_AnnotationSCType.rds")

# reference dataset shard by the authors
ref <- readRDS("../../data/mouse_global_high_res.rds")

# wrangling ---------------------------------------------------------------
barcodes_GEO <- data.frame(barcode = data.combined@meta.data$barcode,
                           ref_GEO = 1)
barcodes_new <- data.frame(barcode = rownames(ref@meta.data),
                           ref_new = 1)

df_join <- full_join(barcodes_GEO,barcodes_new)

dim(df_join)

df_join %>%
  group_by(ref_GEO,ref_new) %>%
  summarise(n = n())

data.combined$missing <- data.combined@meta.data %>%
  mutate(missing = case_when(barcode %in% rownames(ref@meta.data)~"present",
                             T ~ "missing")) %>%
  pull(missing)

data.combined@meta.data %>%
  filter(missing == "missing") %>%
  group_by(Annotation.paper,Sample.paper) %>%
  summarise(n = n())

DimPlot(data.combined,group.by = "RNA_snn_res.0.1",split.by = "missing")

table(ref@meta.data$CellType)

# add the annotation ------------------------------------------------------
ref_meta <- ref@meta.data %>%
  rownames_to_column("barcode") %>%
  select(barcode,cluster_high_res,CellType.l2 = CellType) %>%
  # add a group annotation
  mutate(CellType.l1 = case_when(CellType.l2 %in% c("gCap-Trans.",
                                                    "gCap-Prolif.",
                                                    "gCap",
                                                    "Arterial",
                                                    "gCap-EGR",
                                                    "aCap",
                                                    "Venous",
                                                    "gCap-Cytoprotective")~"Endo",
                                 CellType.l2 %in% c("Clara cell",
                                                    "AT2 cell",
                                                    "AT1 cell")~"Epi",
                                 CellType.l2 %in% c("B-cell 1",
                                                    "NK cell",
                                                    "Cd8+ T cell",
                                                    "B-cell 2",
                                                    "Mixed T cell",
                                                    "Pre-T cell",
                                                    "Cd4+ T cell",
                                                    "B-cell 3")~"Lymphoid",
                                 CellType.l2 %in% c("Trans. Mono",
                                                    "Alv. Mac",
                                                    "NC Mono",
                                                    "Neutrophils",
                                                    "Int. Mac",
                                                    "Classical Mono",
                                                    "Mast cell",
                                                    "Cd103+ cDC",
                                                    "Mreg+ DC",
                                                    "Mono-DC")~"Myeloid",
                                 CellType.l2 %in% c("Pericytes",
                                                    "Matrix Fib",
                                                    "Col13a1+ Fib",
                                                    "Mesothelial",    
                                                    "Col14a1+ Fib",
                                                    "SMC")~"Stromal"
                                 ))






# add the new annotation (partial) to the complete object
df_meta_full <- left_join(data.combined@meta.data,ref_meta,by = "barcode")

# what is the top annotation per each cluster of RNA_snn_res.0.1
df_CellType.l1 <- df_meta_full %>%
  group_by(RNA_snn_res.0.1,CellType.l1) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  ungroup() %>%
  # recode the NA as missing
  mutate(CellType.l1_recode = case_when(is.na(CellType.l1) ~ "missing",
                                        T~CellType.l1))

df_CellType.l1_wide <- df_CellType.l1 %>%
  select(-c(n,tot,CellType.l1)) %>%
  pivot_wider(names_from = RNA_snn_res.0.1,values_from = "prop",values_fill = 0) %>%
  column_to_rownames("CellType.l1_recode")

# plot as heatmap
pdf("../../out/plot/heatmap_CellType.l1.pdf",width = 5,height = 3)
Heatmap(df_CellType.l1_wide,
        name = "Proportion",
        col = viridis::viridis(option = "turbo",n = 10),
        # col = colorRamp2(c(-1, 0, 1), colors = c("blue", "white", "red")),
        row_names_side = "right",
        row_names_gp = gpar(fontsize = 8),
        column_names_side = "bottom",
        column_names_gp = gpar(fontsize = 8),
        row_dend_reorder = FALSE,
        column_dend_reorder = FALSE,
        row_title_gp = gpar(fontsize = 10, fontface = "bold"),
        column_title_gp = gpar(fontsize = 10, fontface = "bold"),
        show_column_names = T,
        show_row_names = T)
dev.off()

# print the top annotation per cluster
df_CellType.l1 %>%
  group_by(RNA_snn_res.0.1) %>%
  arrange(RNA_snn_res.0.1,desc(prop)) %>%
  slice(1)

# use the l2 annotation
df_CellType.l2 <- df_meta_full %>%
  group_by(RNA_snn_res.0.1,CellType.l2) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  ungroup() %>%
  # recode the NA as missing
  mutate(CellType.l2_recode = case_when(is.na(CellType.l2) ~ "missing",
                                        T~CellType.l2))

df_CellType.l2_wide <- df_CellType.l2 %>%
  select(-c(n,tot,CellType.l2)) %>%
  pivot_wider(names_from = RNA_snn_res.0.1,values_from = "prop",values_fill = 0) %>%
  column_to_rownames("CellType.l2_recode")

# plot as heatmap
pdf("../../out/plot/heatmap_CellType.l2.pdf",width = 5,height = 6)
Heatmap(df_CellType.l2_wide,
        name = "Proportion",
        col = viridis::viridis(option = "turbo",n = 10),
        # col = colorRamp2(c(-1, 0, 1), colors = c("blue", "white", "red")),
        row_names_side = "right",
        row_names_gp = gpar(fontsize = 8),
        column_names_side = "bottom",
        column_names_gp = gpar(fontsize = 8),
        row_dend_reorder = FALSE,
        column_dend_reorder = FALSE,
        row_title_gp = gpar(fontsize = 10, fontface = "bold"),
        column_title_gp = gpar(fontsize = 10, fontface = "bold"),
        show_column_names = T,
        show_row_names = T)
dev.off()

# print the top annotation per cluster
df_CellType.l1 %>%
  group_by(RNA_snn_res.0.1) %>%
  arrange(RNA_snn_res.0.1,desc(prop)) %>%
  slice(1)

# for cluster 10 and 12 the annotatio is mostly missing. Use the partial inofrmation from the annotated cells to infer the annotation of the whole cluster
df_LUT_clusters_01 <- df_meta_full %>%
  group_by(RNA_snn_res.0.1,CellType.l1,annotation_confident) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(RNA_snn_res.0.1) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  ungroup() %>%
  # recode the NA as missing
  mutate(CellType.l1_recode = case_when(is.na(CellType.l1) ~ "missing",
                                        T~CellType.l1)) %>%
  group_by(RNA_snn_res.0.1) %>%
  arrange(RNA_snn_res.0.1,desc(prop)) %>%
  slice(1)

df_LUT_clusters_02 <- df_meta_full %>%
  group_by(RNA_snn_res.0.1,CellType.l1,annotation_confident) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(RNA_snn_res.0.1) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  ungroup() %>%
  # recode the NA as missing
  mutate(CellType.l1_recode_cons = case_when(is.na(CellType.l1) ~ "missing",
                                        T~CellType.l1)) %>%
  group_by(RNA_snn_res.0.1) %>%
  arrange(RNA_snn_res.0.1,desc(prop)) %>%
  # filter out the missing CellType.l1
  filter(CellType.l1_recode_cons != "missing") %>%
  slice(1)

df_LUT_clusters <- df_LUT_clusters_01 %>%
  left_join(df_LUT_clusters_02 %>% select(RNA_snn_res.0.1,CellType.l1_recode_cons),by = "RNA_snn_res.0.1")

# according to the new annotation cluster 10 is mostly Epi and cluster 12 is mostly Myeloid. this is in keepign tieth teh depiction form the heatmap

# update the annotation in the main object
data.combined@meta.data <- data.combined@meta.data %>%
  # make sure not to lose the rownames
  rownames_to_column() %>%
  # add the l1 consensus annotation
  left_join(df_LUT_clusters %>% select(RNA_snn_res.0.1,CellType.l1_recode,CellType.l1_recode_cons),by = "RNA_snn_res.0.1") %>%
  # add the l2 annotation
  left_join(df_meta_full %>% select(barcode,CellType.l2),by = "barcode") %>%
  column_to_rownames()

# save the final version of the object
saveRDS(data.combined,"../../out/object/data.combined_harmonySkipIntegration_manualClean_AnnotationSCType_manual.rds")
