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
VlnPlot(object = data.combined,features = GOI)
ggsave("../../out/plot/VlnProcr_cluster.pdf",width = 6,height = 4)

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

# violin global expression
left_join(data.combined@meta.data,df_exp,by = c("barcode")) %>%
  ggplot(aes(x=RNA_snn_res.0.1,y = count_fix)) + geom_violin()+geom_point(position = position_jitter(width = 0.2),size=0.1)+theme_bw()+scale_y_sqrt()+facet_wrap(~gene)+
  theme(strip.background = element_blank())

# prop of positive cells sample
LUT_tot_sample <- left_join(data.combined@meta.data,df_exp,by = c("barcode")) %>%
  group_by(gene,Annotation.paper,Sample.paper) %>%
  summarise(n_tot = n())

left_join(data.combined@meta.data,df_exp,by = c("barcode")) %>%
  group_by(gene,Annotation.paper,Sample.paper,RNA_snn_res.0.1) %>%
  summarise(n_pos = sum(exp_cat=="pos")) %>%
  left_join(LUT_tot_sample,by = c("gene","Annotation.paper","Sample.paper")) %>%
  mutate(prop_pos_sample = n_pos/n_tot) %>%
  group_by(RNA_snn_res.0.1) %>%
  mutate(group_avg = mean(prop_pos_sample)) %>%
  ungroup() %>%
  mutate(RNA_snn_res.0.1 = fct_reorder(RNA_snn_res.0.1,group_avg,.desc = T)) %>%
  ggplot(aes(x=RNA_snn_res.0.1,y = prop_pos_sample)) + geom_boxplot(outlier.shape = NA)+geom_point(position = position_jitter(width = 0.2),shape = 1)+theme_bw()+facet_wrap(~gene)+
  theme(strip.background = element_blank())
ggsave("../../out/plot/propProcrPos_sample.pdf",width = 6,height = 4)


# prop pos cells cluster
left_join(data.combined@meta.data,df_exp,by = c("barcode")) %>%
  group_by(gene,Annotation.paper,Sample.paper,RNA_snn_res.0.1) %>%
  summarise(n_tot_cluster = n(),
            n_pos = sum(exp_cat=="pos")) %>%
  mutate(prop_pos_cluster = n_pos/n_tot_cluster) %>%
  group_by(RNA_snn_res.0.1) %>%
  mutate(group_avg = mean(prop_pos_cluster)) %>%
  ungroup() %>%
  mutate(RNA_snn_res.0.1 = fct_reorder(RNA_snn_res.0.1,group_avg,.desc = T)) %>%
  ggplot(aes(x=RNA_snn_res.0.1,y = prop_pos_cluster)) + geom_boxplot(outlier.shape = NA)+geom_point(position = position_jitter(width = 0.2),shape = 1)+theme_bw()+facet_wrap(~gene)+
  theme(strip.background = element_blank())
ggsave("../../out/plot/propProcrPos_cluster.pdf",width = 6,height = 4)

# average expression global
# build the grouping variable
group_id <- paste(data.combined@meta.data$Annotation.paper,
                  data.combined@meta.data$Sample.paper,
                  paste0("clu_",data.combined@meta.data$RNA_snn_res.0.1),sep = "|")
# add it to the meta
data.combined$group_id <- group_id
# set the idents
Idents(data.combined) <- "group_id"

df_average <- AverageExpression(data.combined,features = GOI) %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(names_to = "group_id",values_to = "avg_exp",-gene)

df_average %>%
  filter(str_detect(group_id,"clu.0"))

# build the pettern to extract the metadata
pattern_Annotation.paper <- paste0(unique(data.combined@meta.data$Annotation.paper),collapse = "|")
pattern_Sample.paper <- paste0(unique(data.combined@meta.data$Sample.paper),collapse = "|") %>% str_replace_all(pattern = "_",replacement = "\\.")

df_average %>%
  mutate(Annotation.paper = str_extract_all(group_id,pattern = pattern_Annotation.paper) %>% unlist()) %>%
  mutate(Sample.paper = str_extract_all(group_id,pattern = pattern_Sample.paper) %>% unlist()) %>%
  mutate(RNA_snn_res.0.1 = str_extract_all(group_id,pattern = "clu.\\d+") %>% unlist()) %>%
  group_by(RNA_snn_res.0.1) %>%
  mutate(group_avg = mean(avg_exp)) %>%
  ungroup() %>%
  mutate(RNA_snn_res.0.1 = fct_reorder(RNA_snn_res.0.1,group_avg,.desc = T)) %>%
  ggplot(aes(x = RNA_snn_res.0.1,y=avg_exp)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitter(width = 0.2),shape = 1)+
  theme_bw() +
  facet_wrap(~gene)+
  theme(strip.background = element_blank())
ggsave("../../out/plot/AvgExpProcr_cluster.pdf",width = 6,height = 4)


# ENDO only ---------------------------------------------------------------
# split the expresion by cell cycle phase

# average expression
# build the grouping variable
group_id_ENDO_phase <- paste(scobj_ENDO@meta.data$Annotation.paper,
                             scobj_ENDO@meta.data$Sample.paper,
                             scobj_ENDO@meta.data$Phase,
                             sep = "|")
# add it to the meta
scobj_ENDO$group_id_phase <- group_id_ENDO_phase

# set the idents
Idents(scobj_ENDO) <- "group_id_phase"

df_average_ENDO_phase <- AverageExpression(scobj_ENDO,features = GOI) %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(names_to = "group_id",values_to = "avg_exp",-gene)

# build the pattern to extract the metadata
pattern_Annotation.paper_phase <- paste0(unique(scobj_ENDO@meta.data$Annotation.paper),collapse = "|")
pattern_Sample.paper_phase <- paste0(unique(scobj_ENDO@meta.data$Sample.paper),collapse = "|") %>% str_replace_all(pattern = "_",replacement = "\\.")
pattern_phase <- paste0(unique(scobj_ENDO@meta.data$Phase),collapse = "|")

df_average_ENDO_phase %>%
  mutate(Annotation.paper = str_extract_all(group_id,pattern = pattern_Annotation.paper_phase) %>% unlist()) %>%
  mutate(Sample.paper = str_extract_all(group_id,pattern = pattern_Sample.paper_phase) %>% unlist()) %>%
  mutate(phase = str_extract_all(group_id,pattern = pattern_phase) %>% unlist()) %>%
  # group_by(RNA_snn_res.0.1) %>%
  # mutate(group_avg = mean(avg_exp)) %>%
  # ungroup() %>%
  # mutate(RNA_snn_res.0.1 = fct_reorder(RNA_snn_res.0.1,group_avg,.desc = T)) %>%
  ggplot(aes(x = phase,y=avg_exp)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2),shape = 1) +
  theme_bw() +
  facet_wrap(~gene)+
  theme(strip.background = element_blank())
ggsave("../../out/plot/AvgExpProcr_Phase.pdf",width = 6,height = 4)

# reinhold suggested to try to connect the dots from the same timepoints
df_average_ENDO_phase %>%
  mutate(Annotation.paper = str_extract_all(group_id,pattern = pattern_Annotation.paper_phase) %>% unlist()) %>%
  mutate(Sample.paper = str_extract_all(group_id,pattern = pattern_Sample.paper_phase) %>% unlist()) %>%
  mutate(phase = str_extract_all(group_id,pattern = pattern_phase) %>% unlist()) %>%
  # group_by(RNA_snn_res.0.1) %>%
  # mutate(group_avg = mean(avg_exp)) %>%
  # ungroup() %>%
  # mutate(RNA_snn_res.0.1 = fct_reorder(RNA_snn_res.0.1,group_avg,.desc = T)) %>%
  ggplot(aes(x = phase,y=avg_exp)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group = Sample.paper,colour = Sample.paper),shape = 1)+
  geom_line(aes(group = paste(Annotation.paper,Sample.paper),colour = Sample.paper))+
  theme_bw() +
  facet_wrap(~gene)+
  theme(strip.background = element_blank())
ggsave("../../out/plot/AvgExpProcr_Phase_line.pdf",width = 6,height = 4)

# try splitting the timepoints
# calculathe the average per group
df_avg_timepoint <- df_average_ENDO_phase %>%
  mutate(Annotation.paper = str_extract_all(group_id,pattern = pattern_Annotation.paper_phase) %>% unlist()) %>%
  mutate(Sample.paper = str_extract_all(group_id,pattern = pattern_Sample.paper_phase) %>% unlist()) %>%
  mutate(phase = str_extract_all(group_id,pattern = pattern_phase) %>% unlist()) %>%
  group_by(Sample.paper) %>%
  summarise(avg = mean(avg_exp))

df_average_ENDO_phase %>%
  mutate(Annotation.paper = str_extract_all(group_id,pattern = pattern_Annotation.paper_phase) %>% unlist()) %>%
  mutate(Sample.paper = str_extract_all(group_id,pattern = pattern_Sample.paper_phase) %>% unlist()) %>%
  mutate(phase = str_extract_all(group_id,pattern = pattern_phase) %>% unlist()) %>%
  # group_by(RNA_snn_res.0.1) %>%
  # mutate(group_avg = mean(avg_exp)) %>%
  # ungroup() %>%
  # mutate(RNA_snn_res.0.1 = fct_reorder(RNA_snn_res.0.1,group_avg,.desc = T)) %>%
  ggplot(aes(x = phase,y=avg_exp)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2),shape = 1) +
  theme_bw() +
  facet_wrap(gene~Sample.paper,nrow = 1) +
  geom_hline(data = df_avg_timepoint,aes(yintercept = avg),linetype = "dashed",col="red")+
  theme(strip.background = element_blank())
ggsave("../../out/plot/AvgExpProcr_Phase_split.pdf",width = 12,height = 4)

# test if it is significant
df_average_ENDO_phase %>%
  mutate(Annotation.paper = str_extract_all(group_id,pattern = pattern_Annotation.paper_phase) %>% unlist()) %>%
  mutate(Sample.paper = str_extract_all(group_id,pattern = pattern_Sample.paper_phase) %>% unlist()) %>%
  mutate(phase = str_extract_all(group_id,pattern = pattern_phase) %>% unlist()) %>%
  lm(data = .,avg_exp~phase) %>%
  summary()

# single cell
VlnPlot(scobj_ENDO,features = GOI,group.by = "Phase")

# test on single cells
unique(so$cell_id)
unique(so$pathology_class)

# set the ident
Idents(scobj_ENDO) <- "Phase"
# vector of the comparisons
phase_id <- str_subset(unique(scobj_ENDO$Phase),pattern = "G1",negate = T)

# run the test
df_res_phase <- lapply(phase_id,function(phase){
    ident_1 <- phase
    ident_2 <- "G1"
    res <- SeuratWrappers::RunPresto(object = scobj_ENDO,
                                     ident.1 = ident_1,
                                     ident.2 = ident_2,
                                     verbose = T,
                                     logfc.threshold = 0,only.pos = F,min.pct = 0.01) %>%
      rownames_to_column("gene")
    return(res)
}) %>%
  setNames(phase_id) %>%
  bind_rows(.id = "comparision") %>%
  mutate(comparison = paste0(comparision,"_vs_G1"))

# fitler out Procr
df_res_phase %>%
  filter(gene %in% c("Procr"))

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
  ggplot(aes(x=Phase,y = count_fix)) + geom_violin()+geom_point(position = position_jitter(width = 0.2),size=0.1)+theme_bw()+scale_y_sqrt()+facet_wrap(~gene)+
  theme(strip.background = element_blank())

# alternative plot
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
  # force the order
  mutate(id = factor(id,level = c("G1","G2M","S"))) %>%
  ggplot(aes(x = features.plot,y = id)) +
  geom_point(aes(size = pct.exp, col = avg.exp.scaled))+
  scale_radius(range = c(0, 8)) +
  # facet_grid(~cell_type,scales = "free",space = "free")+
  theme_cowplot()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 90))+
  scale_color_gradient(low = "lightgrey",high = "blue")

ggsave("../../out/plot/DotplotProcr_phase.pdf",width = 4,height = 4)
