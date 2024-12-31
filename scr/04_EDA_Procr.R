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
data.combined <- readRDS("../../out/object/02_data.combined_harmonySkipIntegration_Annotation.rds")
Idents(data.combined)<-"expertAnno.l1"

# subset only the presumed ENDO cells
# scobj_ENDO <- subset(data.combined,subset = expertAnno.l1 %in% c("ENDO"))
# scobj_ENDO <- data.combined

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
  separate(barcodes,into = c("treat","tissue","barcode_id"),sep = "_",remove = F)

head(df_exp)

# # violin global expression
# left_join(data.combined@meta.data,df_exp,by = c("barcode")) %>%
#   ggplot(aes(x=RNA_snn_res.0.1,y = count_fix)) + geom_violin()+geom_point(position = position_jitter(width = 0.2),size=0.1)+theme_bw()+scale_y_sqrt()+facet_wrap(~gene)+
#   theme(strip.background = element_blank())

# # prop of positive cells sample
# LUT_tot_sample <- left_join(data.combined@meta.data,df_exp,by = c("barcode")) %>%
#   group_by(gene,Annotation.paper,Sample.paper) %>%
#   summarise(n_tot = n())
# 
# left_join(data.combined@meta.data,df_exp,by = c("barcode")) %>%
#   group_by(gene,Annotation.paper,Sample.paper,RNA_snn_res.0.1) %>%
#   summarise(n_pos = sum(exp_cat=="pos")) %>%
#   left_join(LUT_tot_sample,by = c("gene","Annotation.paper","Sample.paper")) %>%
#   mutate(prop_pos_sample = n_pos/n_tot) %>%
#   group_by(RNA_snn_res.0.1) %>%
#   mutate(group_avg = mean(prop_pos_sample)) %>%
#   ungroup() %>%
#   mutate(RNA_snn_res.0.1 = fct_reorder(RNA_snn_res.0.1,group_avg,.desc = T)) %>%
#   ggplot(aes(x=RNA_snn_res.0.1,y = prop_pos_sample)) + geom_boxplot(outlier.shape = NA)+geom_point(position = position_jitter(width = 0.2),shape = 1)+theme_bw()+facet_wrap(~gene)+
#   theme(strip.background = element_blank())
# ggsave("../../out/plot/propProcrPos_sample.pdf",width = 6,height = 4)
# 
# 
# # prop pos cells cluster
# left_join(data.combined@meta.data,df_exp,by = c("barcode")) %>%
#   group_by(gene,Annotation.paper,Sample.paper,RNA_snn_res.0.1) %>%
#   summarise(n_tot_cluster = n(),
#             n_pos = sum(exp_cat=="pos")) %>%
#   mutate(prop_pos_cluster = n_pos/n_tot_cluster) %>%
#   group_by(RNA_snn_res.0.1) %>%
#   mutate(group_avg = mean(prop_pos_cluster)) %>%
#   ungroup() %>%
#   mutate(RNA_snn_res.0.1 = fct_reorder(RNA_snn_res.0.1,group_avg,.desc = T)) %>%
#   ggplot(aes(x=RNA_snn_res.0.1,y = prop_pos_cluster)) + geom_boxplot(outlier.shape = NA)+geom_point(position = position_jitter(width = 0.2),shape = 1)+theme_bw()+facet_wrap(~gene)+
#   theme(strip.background = element_blank())
# ggsave("../../out/plot/propProcrPos_cluster.pdf",width = 6,height = 4)

# # average expression global
# # build the grouping variable
# group_id <- paste(data.combined@meta.data$Annotation.paper,
#                   data.combined@meta.data$Sample.paper,
#                   paste0("clu_",data.combined@meta.data$RNA_snn_res.0.1),sep = "|")
# # add it to the meta
# data.combined$group_id <- group_id
# # set the idents
# Idents(data.combined) <- "group_id"
# 
# df_average <- AverageExpression(data.combined,features = GOI) %>%
#   data.frame() %>%
#   rownames_to_column("gene") %>%
#   pivot_longer(names_to = "group_id",values_to = "avg_exp",-gene)
# 
# df_average %>%
#   filter(str_detect(group_id,"clu.0"))
# 
# # build the pettern to extract the metadata
# pattern_Annotation.paper <- paste0(unique(data.combined@meta.data$Annotation.paper),collapse = "|")
# pattern_Sample.paper <- paste0(unique(data.combined@meta.data$Sample.paper),collapse = "|") %>% str_replace_all(pattern = "_",replacement = "\\.")
# 
# df_average %>%
#   mutate(Annotation.paper = str_extract_all(group_id,pattern = pattern_Annotation.paper) %>% unlist()) %>%
#   mutate(Sample.paper = str_extract_all(group_id,pattern = pattern_Sample.paper) %>% unlist()) %>%
#   mutate(RNA_snn_res.0.1 = str_extract_all(group_id,pattern = "clu.\\d+") %>% unlist()) %>%
#   group_by(RNA_snn_res.0.1) %>%
#   mutate(group_avg = mean(avg_exp)) %>%
#   ungroup() %>%
#   mutate(RNA_snn_res.0.1 = fct_reorder(RNA_snn_res.0.1,group_avg,.desc = T)) %>%
#   ggplot(aes(x = RNA_snn_res.0.1,y=avg_exp)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitter(width = 0.2),shape = 1)+
#   theme_bw() +
#   facet_wrap(~gene)+
#   theme(strip.background = element_blank())
# ggsave("../../out/plot/AvgExpProcr_cluster.pdf",width = 6,height = 4)


# ENDO only ---------------------------------------------------------------
# split the expresion by cell cycle phase

# average expression
# build the grouping variable
group_id_ENDO_treat <- paste(data.combined@meta.data$expertAnno.l1,
                             data.combined@meta.data$orig.ident,
                             sep = "|")
# add it to the meta
data.combined$group_id_treat <- group_id_ENDO_treat

# set the idents
Idents(data.combined) <- "group_id_treat"

df_average_ENDO_treat <- AverageExpression(data.combined,features = GOI) %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(names_to = "group_id",values_to = "avg_exp",-gene)

# build the pattern to extract the metadata
pattern_Annotation.paper_treat <- paste0(unique(data.combined@meta.data$expertAnno.l1),collapse = "|")
pattern_Sample.paper_treat <- paste0(unique(data.combined@meta.data$orig.ident),collapse = "|") %>% str_replace_all(pattern = "_",replacement = "\\.")

df_average_ENDO_treat %>%
  mutate(Annotation.paper = str_extract_all(group_id,pattern = pattern_Annotation.paper_treat) %>% unlist()) %>%
  mutate(Sample.paper = str_extract_all(group_id,pattern = pattern_Sample.paper_treat) %>% unlist()) %>%
  # group_by(RNA_snn_res.0.1) %>%
  # mutate(group_avg = mean(avg_exp)) %>%
  # ungroup() %>%
  # mutate(RNA_snn_res.0.1 = fct_reorder(RNA_snn_res.0.1,group_avg,.desc = T)) %>%
  ggplot(aes(x = Sample.paper,y=avg_exp)) +
  geom_boxplot(outlier.shape = NA) +
  # geom_point(position = position_jitter(width = 0.2),shape = 1)+
  geom_point(aes(col=Annotation.paper))+
  geom_line(aes(group = Annotation.paper,col=Annotation.paper))+
  theme_bw() +
  facet_wrap(~gene)+
  theme(strip.background = element_blank())
ggsave("../../out/plot/04_AvgExpProcr_condition.pdf",width = 6,height = 4)

# test if it is significant
df_average_ENDO_treat %>%
  mutate(Annotation.paper = str_extract_all(group_id,pattern = pattern_Annotation.paper_treat) %>% unlist()) %>%
  mutate(Sample.paper = str_extract_all(group_id,pattern = pattern_Sample.paper_treat) %>% unlist()) %>%
  lm(data = .,avg_exp~Sample.paper) %>%
  summary()

# single cell
VlnPlot(data.combined,features = GOI,group.by = "expertAnno.l1")

# alternative plot, tailor the order of the rows
test_dotplot <- DotPlot(data.combined,features = GOI)

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
  separate(id,into = c("expertAnno.l1","treat"),sep = "\\|",remove = F) %>%
  # force the order
  # mutate(id = factor(id,level = c("Ctrl","DPT_3d","DPT_5d","DTP_7d"))) %>%
  ggplot(aes(x = treat,y = expertAnno.l1)) +
  geom_point(aes(size = pct.exp, col = avg.exp.scaled))+
  scale_radius(range = c(0, 8)) +
  # facet_grid(~cell_type,scales = "free",space = "free")+
  theme_cowplot()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 90))+
  scale_color_gradient(low = "lightgrey",high = "blue")

# ggsave("../../out/plot/DotplotProcr_condition.pdf",width = 4,height = 4)


# -------------------------------------------------------------------------
# test run DE analsys at the single cell level
# run the comparison per cell type, regardless of the treatment
Idents(data.combined) <- "expertAnno.l1"
res_cellid <- SeuratWrappers::RunPrestoAll(data.combined,verbose = T,logfc.threshold = 0.25,only.pos = F,min.pct = 0.25)

# check the resutl
res_cellid |>
  filter(gene %in% GOI) |>
  ggplot(aes(x = avg_log2FC,y = -log(p_val_adj),col=cluster,label=cluster))+geom_point()+
  facet_wrap(~gene,scales = "free")+theme_bw()+theme(strip.background = element_blank())+
  geom_vline(xintercept = c(0),linetype="dashed",col="gray")+
  geom_vline(xintercept = c(-0.5,0.5),linetype="dashed",col="red")+
  geom_hline(yintercept = c(-log(0.05)),linetype="dashed",col="red")+
  ggrepel::geom_text_repel(force = 100,segment.alpha=0.5,max.overlaps = 5)

res_cellid %>%
  filter(gene %in% c("Procr"))

# try to run the comparisoin regardless of the cell type
Idents(data.combined) <- "orig.ident"
res_condition <- SeuratWrappers::RunPrestoAll(data.combined,verbose = T,logfc.threshold = 0.25,only.pos = F,min.pct = 0.25)

# check the resutl
res_condition |>
  filter(gene %in% GOI) |>
  ggplot(aes(x = avg_log2FC,y = -log(p_val_adj),col=cluster,label=cluster))+geom_point()+
  facet_wrap(~gene,scales = "free")+theme_bw()+theme(strip.background = element_blank())+
  geom_vline(xintercept = c(0),linetype="dashed",col="gray")+
  geom_vline(xintercept = c(-0.5,0.5),linetype="dashed",col="red")+
  geom_hline(yintercept = c(-log(0.05)),linetype="dashed",col="red")+
  ggrepel::geom_text_repel(force = 100,segment.alpha=0.5,max.overlaps = 5)

res_cellid %>%
  filter(gene %in% c("Procr"))

# run the comparison by condition and cell type
# check number of cells
data.combined@meta.data %>%
  group_by(expertAnno.l1, orig.ident) %>%
  summarise(n = n(), .groups = "drop") %>%  # Summarize with n()
  complete(expertAnno.l1, orig.ident, fill = list(n = 0)) 

id_cell_type <- data.combined@meta.data$expertAnno.l1 %>%
  unique() %>%
  str_subset(pattern = "ENDO|FIBRO",negate = F)

# skip the cell types where there are too few cells
# x <- "ENDO"
df_res <- lapply(id_cell_type, function(x){
  # track the progress
  print(x)
  # set the ident
  scobj_subset <- subset(data.combined,subset = expertAnno.l1 %in% c(x))
  Idents(scobj_subset) <- "orig.ident"
  # vector of the comparisons
  treat_id <- str_subset(unique(scobj_subset$orig.ident),pattern = "young",negate = T)
  
  # run the test. rerutn all the genes
  df_res_treat <- lapply(treat_id,function(treat){
    ident_1 <- treat
    ident_2 <- "young_aorta"
    res <- SeuratWrappers::RunPresto(object = scobj_subset,
                                     ident.1 = ident_1,
                                     ident.2 = ident_2,
                                     verbose = T,
                                     logfc.threshold = 0,only.pos = F,min.pct = 0.01) %>%
      rownames_to_column("gene")
    return(res)
  }) %>%
    setNames(treat_id) %>%
    bind_rows(.id = "test") %>%
    mutate(expertAnno.l1 = x) %>%
    mutate(comparison = paste0(test,"_vs_young_aorta"))
  
  return(df_res_treat)
}) %>%
  bind_rows()


# fitler out Procr
df_res %>%
  filter(gene %in% c("Procr"))

# plot the volcano for the cluster 12 DE
volcano_tot <- df_res %>%
  mutate(DE_cat = case_when(avg_log2FC > 1 & p_val_adj < 0.01~"up",
                            avg_log2FC < (-1) & p_val_adj < 0.01~"down",
                            T~"no"))
# filter(cluster %in% c("IMM"))

ggplot() +
  geom_point(data = volcano_tot%>%filter(DE_cat=="no"),aes(x=avg_log2FC,y=-log10(p_val_adj)),size=1,alpha=0.2)+theme_bw()+
  geom_point(data = volcano_tot%>%filter(DE_cat!="no"),aes(x=avg_log2FC,y=-log10(p_val_adj)),size=1,alpha=0.2,col="red")+
  theme_bw()+
  ggrepel::geom_text_repel(data = volcano_tot%>%filter(DE_cat!="no"),aes(x=avg_log2FC,y=-log10(p_val_adj),label=gene))+
  theme_bw()+
  facet_grid(comparison~expertAnno.l1)+
  theme(strip.background = element_blank())
ggsave("../../out/image/volcano_DE_treatvsCSFctrl_cellIDwise_minpct5.pdf",width = 25,height = 20)

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
  separate(barcodes,into = c("treat","tissue","barcode_id"),sep = "_",remove = F)

head(df_exp_ENDO)

# violin global expression
left_join(scobj_ENDO@meta.data %>% rownames_to_column("barcodes"),df_exp_ENDO,by = c("barcodes")) %>%
  ggplot(aes(x=orig.ident,y = count_fix)) + geom_violin()+geom_point(position = position_jitter(width = 0.2),size=0.1)+theme_bw()+scale_y_sqrt()+facet_wrap(~gene)+
  theme(strip.background = element_blank())

# # alternative plot, tailor the order of the rows
# test_dotplot <- DotPlot(scobj_ENDO,features = GOI)
# 
# # generate the ggplot object
# # force the order
# df_test <- lapply(GOI,function(x){
#   test_dotplot$data %>% 
#     filter(features.plot %in% x)
# }) %>% 
#   bind_rows(.id = "cell_type")
# 
# glimpse(df_test)
# 
# # plot mapping radius
# df_test %>%
#   # force the order
#   # mutate(id = factor(id,level = c("Ctrl","DPT_3d","DPT_5d","DTP_7d"))) %>%
#   ggplot(aes(x = features.plot,y = id)) +
#   geom_point(aes(size = pct.exp, col = avg.exp.scaled))+
#   scale_radius(range = c(5, 6)) +
#   # facet_grid(~cell_type,scales = "free",space = "free")+
#   theme_cowplot()+
#   theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 90))+
#   scale_color_gradient(low = "lightgrey",high = "blue")
# 
# ggsave("../../out/plot/DotplotProcr_condition.pdf",width = 4,height = 4)
