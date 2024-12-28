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
data.combined <- readRDS("../../out/object/data.combined_harmonySkipIntegration_manualClean.rds")
Idents(data.combined)<-"RNA_snn_res.0.1"

# subset only the presumed ENDO cells
scobj_ENDO <- subset(data.combined,subset = RNA_snn_res.0.1 %in% c(0,9))

GOI <- c("Procr")

# EDA ---------------------------------------------------------------------
# check the global Procr expression
# VlnPlot(object = data.combined,features = GOI)
# ggsave("../../out/plot/VlnProcr_cluster.pdf",width = 6,height = 4)

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

# ENDO only ---------------------------------------------------------------
# split the expresion by cell cycle phase

# average expression
# build the grouping variable
group_id_ENDO_treat <- paste(scobj_ENDO@meta.data$Annotation.paper,
                             scobj_ENDO@meta.data$Sample.paper,
                             scobj_ENDO@meta.data$Phase,
                             sep = "|")
# add it to the meta
scobj_ENDO$group_id_treat <- group_id_ENDO_treat

# set the idents
Idents(scobj_ENDO) <- "group_id_treat"

df_average_ENDO_treat <- AverageExpression(scobj_ENDO,features = GOI) %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(names_to = "group_id",values_to = "avg_exp",-gene)

# build the pattern to extract the metadata
pattern_Annotation.paper_treat <- paste0(unique(scobj_ENDO@meta.data$Annotation.paper),collapse = "|")
pattern_Sample.paper_treat <- paste0(unique(scobj_ENDO@meta.data$Sample.paper),collapse = "|") %>% str_replace_all(pattern = "_",replacement = "\\.")

pattern_Phase_treat <- paste0(unique(scobj_ENDO@meta.data$Phase),collapse = "|") %>% str_replace_all(pattern = "_",replacement = "\\.")

df_average_ENDO_treat %>%
  mutate(Annotation.paper = str_extract_all(group_id,pattern = pattern_Annotation.paper_treat) %>% unlist()) %>%
  mutate(Sample.paper = str_extract_all(group_id,pattern = pattern_Sample.paper_treat) %>% unlist()) %>%
  mutate(Phase = str_extract_all(group_id,pattern = pattern_Phase_treat) %>% unlist()) %>%
  # group_by(RNA_snn_res.0.1) %>%
  # mutate(group_avg = mean(avg_exp)) %>%
  # ungroup() %>%
  # mutate(RNA_snn_res.0.1 = fct_reorder(RNA_snn_res.0.1,group_avg,.desc = T)) %>%
  ggplot(aes(x = Sample.paper,y=avg_exp)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitter(width = 0.2),shape = 1)+
  theme_bw() +
  facet_wrap(Phase~gene)+
  theme(strip.background = element_blank())
ggsave("../../out/plot/AvgExpProcr_CondPhase.pdf",width = 12,height = 4)

# Reinhold suggested to try to color by cell phase rather than splitting
df_average_ENDO_treat %>%
  mutate(Annotation.paper = str_extract_all(group_id,pattern = pattern_Annotation.paper_treat) %>% unlist()) %>%
  mutate(Sample.paper = str_extract_all(group_id,pattern = pattern_Sample.paper_treat) %>% unlist()) %>%
  mutate(Phase = str_extract_all(group_id,pattern = pattern_Phase_treat) %>% unlist()) %>%
  # group_by(RNA_snn_res.0.1) %>%
  # mutate(group_avg = mean(avg_exp)) %>%
  # ungroup() %>%
  # mutate(RNA_snn_res.0.1 = fct_reorder(RNA_snn_res.0.1,group_avg,.desc = T)) %>%
  ggplot(aes(x = Sample.paper,y=avg_exp,col=Phase)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2,dodge.width = 0.8),shape = 1)+
  theme_bw() +
  facet_wrap(~gene)+
  theme(strip.background = element_blank())
ggsave("../../out/plot/AvgExpProcr_CondPhase_02.pdf",width = 5,height = 4)

# test if it is significant
df_average_ENDO_treat %>%
  mutate(Annotation.paper = str_extract_all(group_id,pattern = pattern_Annotation.paper_treat) %>% unlist()) %>%
  mutate(Sample.paper = str_extract_all(group_id,pattern = pattern_Sample.paper_treat) %>% unlist()) %>%
  mutate(Phase = str_extract_all(group_id,pattern = pattern_Phase_treat) %>% unlist()) %>%
  lm(data = .,avg_exp~Sample.paper+Phase) %>%
  summary()

# single cell
# scobj_ENDO$test <- paste0(scobj_ENDO$Sample.paper,"|",scobj_ENDO$Phase)
VlnPlot(scobj_ENDO,features = GOI,group.by = "Sample.paper",split.by = "Phase",split.plot = T)

# set the ident
scobj_ENDO$test <- paste0(scobj_ENDO$Sample.paper,"|",scobj_ENDO$Phase)
Idents(scobj_ENDO) <- "test"
# vector of the comparisons
treat_id <- str_subset(unique(scobj_ENDO$test),pattern = "Ctrl\\|G1",negate = T)

# run the test
df_res_treat <- lapply(treat_id,function(treat){
  ident_1 <- treat
  ident_2 <- "Ctrl|G1"
  res <- SeuratWrappers::RunPresto(object = scobj_ENDO,
                                   ident.1 = ident_1,
                                   ident.2 = ident_2,
                                   verbose = T,
                                   logfc.threshold = 0,only.pos = F,min.pct = 0.01) %>%
    rownames_to_column("gene")
  return(res)
}) %>%
  setNames(treat_id) %>%
  bind_rows(.id = "comparision") %>%
  mutate(comparison = paste0(comparision,"_vs_Ctrl|G1"))

# fitler out Procr
df_res_treat %>%
  filter(gene %in% c("Procr")) %>%
  arrange(p_val_adj)

# manual plotting of the sc data
# extract the expression data
df_exp_ENDO <- FetchData(scobj_ENDO, vars = GOI,layer = "data") |> 
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

head(df_exp_ENDO)

# violin global expression
left_join(scobj_ENDO@meta.data,df_exp_ENDO,by = c("barcode")) %>%
  ggplot(aes(x=Sample.paper,y = count_fix)) +
  geom_violin()+
  geom_point(position = position_jitter(width = 0.2),size=0.1)+
  theme_bw()+scale_y_sqrt()+
  facet_wrap(gene~Phase)+
  theme(strip.background = element_blank())

# alternative plot, tailor the order of the rows
test_dotplot <- DotPlot(scobj_ENDO,features = GOI)

# generate the ggplot object
# force the order
df_test <- lapply(GOI,function(x){
  test_dotplot$data %>% 
    filter(features.plot %in% x)
}) %>% 
  bind_rows(.id = "cell_type")

glimpse(df_test)

# plot mapping radius
df_test %>%
  separate(id,into = c("timepoint","phase"),sep = "\\|",remove = F) %>%
  # force the order
  mutate(timepoint = factor(timepoint,level = c("Ctrl","DPT_3d","DPT_5d","DTP_7d"))) %>%
  mutate(phase = factor(phase,level = c("G1","G2M","S"))) %>%
  ggplot(aes(x = phase,y = timepoint)) +
  geom_point(aes(size = pct.exp, col = avg.exp.scaled))+
  scale_radius(range = c(0, 8)) +
  # facet_grid(~cell_type,scales = "free",space = "free")+
  theme_cowplot()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 90))+
  scale_color_gradient(low = "lightgrey",high = "blue")

ggsave("../../out/plot/DotplotProcr_CondPhase.pdf",width = 4,height = 4)
