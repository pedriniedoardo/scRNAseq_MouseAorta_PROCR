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
library(SeuratWrappers)
library(presto)

# read in the data --------------------------------------------------------
data.combined <- readRDS("../../out/object/data.combined_harmonySkipIntegration_manualClean_AnnotationSCType_manual.rds")
Idents(data.combined)<-"CellType.l1_recode_cons"

GOI <- c("Procr")

# EDA ---------------------------------------------------------------------
# check the global Procr expression
# VlnPlot(object = data.combined,features = GOI)
# ggsave("../../out/plot/VlnProcr_cluster.pdf",width = 6,height = 4)

# default plotting --------------------------------------------------------
show_col(brewer.pal(5,"Set1"))
# plot the UMAP with the final annotation
DimPlot(data.combined,group.by = "CellType.l1_recode_cons",cols = brewer.pal(5,"Set1"))
ggsave("../../out/plot/UMAP_CellType_l1_recode_cons.pdf",width = 6,height = 5)
FeaturePlot(data.combined,features = GOI,order = T)
ggsave("../../out/plot/UMAP_CellType_l1_recode_cons_Procr.pdf",width = 6,height = 5)

# Reinhold requested to try to split by timpoints.
FeaturePlot(data.combined,features = GOI,order = T,split.by = "Sample.paper")
ggsave("../../out/plot/UMAP_CellType_l1_recode_cons_Procr_split.pdf",width = 20,height = 5)

# manual plotting ---------------------------------------------------------
# extract the expression data
df_exp <- FetchData(data.combined, vars = GOI,layer = "data") |> 
  rownames_to_column("barcodes") |> 
  pivot_longer(names_to = "gene",values_to = "count",-barcodes) |> 
  # try to min/max normalize the count varaible per gene in order to rescale the difference in terms of expression
  group_by(gene) %>%
  # threshold of positiveness is based on the distriubtion of the expression of the signal in tihs case
  mutate(norm_min_max = ((count - min(count))/(max(count)-min(count))),
         exp_cat = case_when(count > 0~"pos",
                             T~"neg")) %>%
  ungroup() %>%
  mutate(count_fix = count + rnorm(nrow(.))/100000) %>%
  separate(barcodes,into = c("barcode","barcode_id"),sep = "-",remove = F)

head(df_exp)

# plot all the cells ------------------------------------------------------
# split the expresion by cell cycle phase

# average expression
# build the grouping variable
group_id <- paste(data.combined@meta.data$Annotation.paper,
                  data.combined@meta.data$Sample.paper,
                  data.combined@meta.data$CellType.l1_recode_cons,
                  sep = "|")
# add it to the meta
data.combined$group_id <- group_id

# set the idents
Idents(data.combined) <- "group_id"

df_average <- AggregateExpression(data.combined,features = GOI) %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(names_to = "group_id",values_to = "avg_exp",-gene)

# build the pattern to extract the metadata
pattern_Annotation.paper <- paste0(unique(data.combined@meta.data$Annotation.paper),collapse = "|")
pattern_Sample.paper <- paste0(unique(data.combined@meta.data$Sample.paper),collapse = "|") %>% str_replace_all(pattern = "_",replacement = "\\.")
pattern_CellType <- paste0(unique(data.combined@meta.data$CellType.l1_recode_cons),collapse = "|") %>% str_replace_all(pattern = "_",replacement = "\\.")

df_average %>%
  mutate(Annotation.paper = str_extract_all(group_id,pattern = pattern_Annotation.paper) %>% unlist()) %>%
  mutate(Sample.paper = str_extract_all(group_id,pattern = pattern_Sample.paper) %>% unlist()) %>%
  mutate(CellType = str_extract_all(group_id,pattern = pattern_CellType) %>% unlist()) %>%
  # group_by(RNA_snn_res.0.1) %>%
  # mutate(group_avg = mean(avg_exp)) %>%
  # ungroup() %>%
  # mutate(RNA_snn_res.0.1 = fct_reorder(RNA_snn_res.0.1,group_avg,.desc = T)) %>%
  ggplot(aes(x = Sample.paper,y=avg_exp)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitter(width = 0.2),shape = 1)+
  theme_bw() +
  facet_wrap(CellType~gene,nrow = 1)+
  theme(strip.background = element_blank())
ggsave("../../out/plot/AvgExpProcr_CellType_l1_recode_cons_Procr_01.pdf",width = 12,height = 3)

df_average %>%
  mutate(Annotation.paper = str_extract_all(group_id,pattern = pattern_Annotation.paper) %>% unlist()) %>%
  mutate(Sample.paper = str_extract_all(group_id,pattern = pattern_Sample.paper) %>% unlist()) %>%
  mutate(CellType = str_extract_all(group_id,pattern = pattern_CellType) %>% unlist()) %>%
  # group_by(RNA_snn_res.0.1) %>%
  # mutate(group_avg = mean(avg_exp)) %>%
  # ungroup() %>%
  # mutate(RNA_snn_res.0.1 = fct_reorder(RNA_snn_res.0.1,group_avg,.desc = T)) %>%
  ggplot(aes(x = Sample.paper,y=avg_exp,color = CellType)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2,dodge.width = 0.8),shape = 1)+
  theme_bw() +
  facet_wrap(~gene,nrow = 1)+
  theme(strip.background = element_blank()) +
  scale_color_manual(values = brewer.pal(5,"Set1"))
ggsave("../../out/plot/AvgExpProcr_CellType_l1_recode_cons_Procr_02.pdf",width = 5,height = 3)

# test if it is significant
df_average %>%
  mutate(Annotation.paper = str_extract_all(group_id,pattern = pattern_Annotation.paper) %>% unlist()) %>%
  mutate(Sample.paper = str_extract_all(group_id,pattern = pattern_Sample.paper) %>% unlist()) %>%
  mutate(CellType = str_extract_all(group_id,pattern = pattern_CellType) %>% unlist()) %>%
  lm(data = .,avg_exp~CellType*Sample.paper) %>%
  summary()

# single cell
# data.combined$test <- paste0(data.combined$Sample.paper,"|",data.combined$Phase)
VlnPlot(data.combined,features = GOI,group.by = "Sample.paper",split.by = "CellType.l1_recode_cons",split.plot = T,,cols = brewer.pal(5,"Set1"))

# manual plotting of the sc data
# extract the expression data
df_exp <- FetchData(data.combined, vars = GOI,layer = "data") |> 
  rownames_to_column("barcodes") |> 
  pivot_longer(names_to = "gene",values_to = "count",-barcodes) |> 
  # try to min/max normalize the count varaible per gene in order to rescale the difference in terms of expression
  group_by(gene) %>%
  # threshold of positiveness is based on the distriubtion of the expression of the signal in tihs case
  mutate(norm_min_max = ((count - min(count))/(max(count)-min(count))),
         exp_cat = case_when(count > 0~"pos",
                             T~"neg")) %>%
  ungroup() %>%
  mutate(count_fix = count + rnorm(nrow(.))/100000) %>%
  separate(barcodes,into = c("barcode","barcode_id"),sep = "-",remove = F)

head(df_exp)

# violin global expression
left_join(data.combined@meta.data,df_exp,by = c("barcode")) %>%
  ggplot(aes(x=Sample.paper,y = count_fix)) +
  geom_violin()+
  geom_point(position = position_jitter(width = 0.2),size=0.1)+
  theme_bw()+scale_y_sqrt()+
  facet_wrap(gene~CellType.l1_recode_cons)+
  theme(strip.background = element_blank())

# alternative plot, tailor the order of the rows
# build the grouping variable
group_id2 <- paste(data.combined@meta.data$Sample.paper,
                   data.combined@meta.data$CellType.l1_recode_cons,
                   sep = "|")
# add it to the meta
data.combined$group_id2 <- group_id2

# set the idents
Idents(data.combined) <- "group_id2"
test_dotplot <- DotPlot(data.combined,features = GOI)

# generate the ggplot object
# force the order
df_test <- lapply(GOI,function(x){
  test_dotplot$data %>% 
    filter(features.plot %in% x)
}) %>% 
  bind_rows(.id = "ref")

glimpse(df_test)

# plot mapping radius
df_test %>%
  separate(id,into = c("timepoint","CellType"),sep = "\\|",remove = F) %>%
  # force the order
  mutate(timepoint = factor(timepoint,level = c("Ctrl","DPT_3d","DPT_5d","DTP_7d"))) %>%
  # mutate(phase = factor(phase,level = c("G1","G2M","S"))) %>%
  ggplot(aes(x = CellType,y = timepoint)) +
  geom_point(aes(size = pct.exp, col = avg.exp.scaled))+
  scale_radius(range = c(0, 8)) +
  # facet_grid(~cell_type,scales = "free",space = "free")+
  theme_cowplot()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 90))+
  scale_color_gradient(low = "lightgrey",high = "blue")+
  ggtitle("Procr expression")
ggsave("../../out/plot/DotplotProcr_CellType_l1_recode_cons_Procr.pdf",width = 5,height = 3)

# reinhold requested to invert the axis of the plot
df_test %>%
  separate(id,into = c("timepoint","CellType"),sep = "\\|",remove = F) %>%
  # force the order
  mutate(timepoint = factor(timepoint,level = c("Ctrl","DPT_3d","DPT_5d","DTP_7d"))) %>%
  # mutate(phase = factor(phase,level = c("G1","G2M","S"))) %>%
  mutate(CellType = factor(CellType,level = c("Stromal","Myeloid","Lymphoid","Epi","Endo"))) %>%
  ggplot(aes(x = timepoint,y = CellType)) +
  geom_point(aes(size = pct.exp, col = avg.exp.scaled))+
  scale_radius(range = c(0, 8)) +
  # facet_grid(~cell_type,scales = "free",space = "free")+
  theme_cowplot()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 90))+
  scale_color_gradient(low = "lightgrey",high = "blue")+
  ggtitle("Procr expression")
ggsave("../../out/plot/DotplotProcr_CellType_l1_recode_cons_Procr2.pdf",width = 5,height = 3)

# set the ident
Idents(data.combined) <- "CellType.l1_recode_cons"

# vector of the comparisons
test_id2 <- unique(Idents(data.combined))

# run the test for the two cell types vs all the other Endo cells
# ident.2: A second identity class for comparison; if NULL, use all other cells for comparison
df_res_treat2 <- lapply(test_id2,function(test){
  ident_1 <- test
  # ident_2 <- "Ctrl"
  res <- SeuratWrappers::RunPresto(object = data.combined,
                                   ident.1 = ident_1,
                                   # ident.2 = ident_2,
                                   verbose = T,
                                   logfc.threshold = 0,only.pos = F,min.pct = 0.01) %>%
    rownames_to_column("gene")
  return(res)
}) %>%
  setNames(test_id2) %>%
  bind_rows(.id = "comparision") %>%
  mutate(comparison = paste0(comparision,"_vs_All"))

# filter Procr stats
df_res_treat2 %>%
  filter(gene %in% GOI)

# ENDO subset -------------------------------------------------------------
# subset the data
# enumerate the different cell types of endo origin
data.combined@meta.data %>%
  filter(CellType.l1_recode_cons %in% c("Endo")) %>%
  group_by(CellType.l2) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  filter(CellType.l2 %in% c("gCap-Trans.",
                            "gCap-Prolif.",
                            "gCap",
                            "Arterial",
                            "gCap-EGR",
                            "aCap",
                            "Venous",
                            "gCap-Cytoprotective"))

data.combined@meta.data %>%
  filter(CellType.l1_recode_cons %in% c("Endo")) %>%
  group_by(CellType.l2) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  filter(CellType.l2 %in% c("gCap-Trans.",
                            "gCap-Prolif.",
                            "gCap",
                            "Arterial",
                            "gCap-EGR",
                            "aCap",
                            "Venous",
                            "gCap-Cytoprotective")
  ) %>%
  summarise(tot = sum(n))

# the congruent cell types are the following
# c("gCap-Trans.",
#   "gCap-Prolif",
#   "gCap",
#   "Arterial",
#   "gCap-EGR",
#   "aCap",
#   "Venous",
#   "gCap-Cytoprotective")
# exclude all the others

data.combined@meta.data %>%
  group_by(CellType.l2,CellType.l1_recode_cons) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  filter(CellType.l2 %in% c("gCap-Trans.",
                            "gCap-Prolif.",
                            "gCap",
                            "Arterial",
                            "gCap-EGR",
                            "aCap",
                            "Venous",
                            "gCap-Cytoprotective")
  )

# define the filtering criteria
data.combined$endo_cell_id <- data.combined@meta.data %>%
  mutate(endo_cell_id = case_when(CellType.l1_recode_cons %in% c("Endo") & CellType.l2 %in% c("gCap-Trans.",
                                                                                              "gCap-Prolif.",
                                                                                              "gCap",
                                                                                              "Arterial",
                                                                                              "gCap-EGR",
                                                                                              "aCap",
                                                                                              "Venous",
                                                                                              "gCap-Cytoprotective")~1,
                                  T~0)) %>%
  pull(endo_cell_id)

# subset the data
scobj_ENDO <- subset(data.combined,subset = endo_cell_id == 1)
# confirm the size is congruent
dim(scobj_ENDO)

scobj_ENDO@meta.data %>%
  group_by(CellType.l2) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  ggplot(aes(x = reorder(CellType.l2,-n),y = n))+geom_col()+theme_bw()+theme(axis.text.x = element_text(angle = 45,hjust = 1))
ggsave("../../out/plot/Barplot_EndoCellTypes.pdf",width = 5,height = 3)

# Check also the relative proportion of each cell type across time points
scobj_ENDO@meta.data %>%
  group_by(Sample.paper,Annotation.paper,CellType.l2) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(Sample.paper,Annotation.paper) %>%
  mutate(tot = sum(n)) %>%
  ungroup() %>%
  mutate(prop = n/tot) %>%
  ggplot(aes(x = CellType.l2,y = prop,col=Sample.paper))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitterdodge(jitter.width = 0.2,dodge.width = 0.8),shape = 1)+
  theme_bw()+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_y_sqrt()
ggsave("../../out/plot/BoxplotProp_EndoCellTypes.pdf",width = 7,height = 5)

# explore if Procr is expressed in any specific Endo annotation
# build the grouping variable
group_id3 <- paste(scobj_ENDO@meta.data$Sample.paper,
                   scobj_ENDO@meta.data$CellType.l2,
                   sep = "|")
# add it to the meta
scobj_ENDO$group_id3 <- group_id3

Idents(scobj_ENDO) <- "group_id3"
test_dotplot3 <- DotPlot(scobj_ENDO,features = GOI)

# generate the ggplot object
# force the order
df_test2 <- lapply(GOI,function(x){
  test_dotplot3$data %>% 
    filter(features.plot %in% x)
}) %>% 
  bind_rows(.id = "ref")

glimpse(df_test2)

# plot mapping radius
df_test2 %>%
  separate(id,into = c("timepoint","CellType"),sep = "\\|",remove = F) %>%
  # force the order
  mutate(timepoint = factor(timepoint,level = c("Ctrl","DPT_3d","DPT_5d","DTP_7d"))) %>%
  # mutate(phase = factor(phase,level = c("G1","G2M","S"))) %>%
  ggplot(aes(x = CellType,y = timepoint)) +
  geom_point(aes(size = pct.exp, col = avg.exp.scaled))+
  scale_radius(range = c(0, 8)) +
  # facet_grid(~cell_type,scales = "free",space = "free")+
  theme_cowplot()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 90))+
  scale_color_gradient(low = "lightgrey",high = "blue")+
  ggtitle("Procr expression")
ggsave("../../out/plot/DotplotProcr_EndoSubset_Procr.pdf",width = 6,height = 5)

# invert the axis
df_test2 %>%
  separate(id,into = c("timepoint","CellType"),sep = "\\|",remove = F) %>%
  # force the order
  mutate(timepoint = factor(timepoint,level = c("Ctrl","DPT_3d","DPT_5d","DTP_7d"))) %>%
  # mutate(phase = factor(phase,level = c("G1","G2M","S"))) %>%
  ggplot(aes(x = timepoint,y = CellType)) +
  geom_point(aes(size = pct.exp, col = avg.exp.scaled))+
  scale_radius(range = c(0, 8)) +
  # facet_grid(~cell_type,scales = "free",space = "free")+
  theme_cowplot()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 90))+
  scale_color_gradient(low = "lightgrey",high = "blue")+
  ggtitle("Procr expression")
ggsave("../../out/plot/DotplotProcr_EndoSubset_Procr2.pdf",width = 5,height = 5)

# Rinhold suggested to run DE analysis at the single cell level to compare
# Endo vs gCAP trans
# Endo vs gCAP prolif
# Endo vs gCap cytoprotective
# subset only the ENDO cells

# run the DE analysis

# set the ident
Idents(scobj_ENDO) <- "CellType.l2"

# vector of the comparisons
test_id <- c("gCap-Trans.","gCap-Cytoprotective","gCap-Prolif.")

# run the test for the two cell types vs all the other Endo cells
# ident.2: A second identity class for comparison; if NULL, use all other cells for comparison
df_res_treat <- lapply(test_id,function(test){
  ident_1 <- test
  # ident_2 <- "Ctrl"
  res <- SeuratWrappers::RunPresto(object = scobj_ENDO,
                                   ident.1 = ident_1,
                                   # ident.2 = ident_2,
                                   verbose = T,
                                   logfc.threshold = 0,only.pos = F,min.pct = 0.01) %>%
    rownames_to_column("gene")
  return(res)
}) %>%
  setNames(test_id) %>%
  bind_rows(.id = "comparision") %>%
  mutate(comparison = paste0(comparision,"_vs_AllEndo"))

# filter Procr stats
df_res_treat %>%
  filter(gene %in% GOI)
