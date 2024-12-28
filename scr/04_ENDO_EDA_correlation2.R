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
library(GGally)
library(ComplexHeatmap)
library(circlize)
library(GSEABase)
library(GGally)
library(ComplexHeatmap)
library(circlize)

# read in the data --------------------------------------------------------
data.combined <- readRDS("../../out/object/data.combined_harmonySkipIntegration_manualClean.rds")
Idents(data.combined)<-"RNA_snn_res.0.1"

# subset only the presumed ENDO cells
scobj_ENDO <- subset(data.combined,subset = RNA_snn_res.0.1 %in% c(0,9))

# read in the genes from the signature Reinhold recommended
vec_symbol <- getGmt("../../data/GOBP_POSITIVE_REGULATION_OF_VASCULAR_ENDOTHELIAL_CELL_PROLIFERATION.v2024.1.Mm.gmt") %>%
  geneIds() %>%
  unique() %>%
  unlist()

# reinhold also asked to add Mki67 in the dataset
# Ghsr is too lowly expressed
GOI <- c(vec_symbol,"Mki67","Procr") %>%
  str_subset("Ghsr|Igf2|Prok1",negate = T)

# ENDO only ---------------------------------------------------------------
# sc level expression
df_exp <- FetchData(scobj_ENDO, vars = GOI,layer = "data") |> 
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

exp_sc_wide <- df_exp %>%
  dplyr::select(barcodes,gene,count_fix) %>%
  pivot_wider(names_from = gene,values_from = count_fix)

pdf("../../out/plot/ggpair_stemness_sc2.pdf",width = 15,height = 15)
ggpairs(exp_sc_wide[,-1])
dev.off()

# try also to correlate the expression of Procr with the module of the other genes
# score the signatures, both long and short
scobj_ENDO <- Seurat::AddModuleScore(scobj_ENDO,
                                     features = list(GOI),
                                     name = "_stemness")

df_exp %>%
  dplyr::select(barcodes,gene,count_fix) %>%
  dplyr::filter(gene == "Procr") %>%
  pivot_wider(names_from = gene,values_from = count_fix) %>%
  left_join(scobj_ENDO@meta.data %>% rownames_to_column("barcodes")) %>%
  ggplot(aes(x=Procr,y=`_stemness1`)) + geom_point(shape = 1)+geom_smooth(method = "lm")+
  theme_bw()
ggsave("../../out/plot/correlation_stemness_Procr2.pdf",width = 4,height = 4)

test_stemness <- df_exp %>%
  dplyr::select(barcodes,gene,count_fix) %>%
  dplyr::filter(gene == "Procr") %>%
  pivot_wider(names_from = gene,values_from = count_fix) %>%
  left_join(scobj_ENDO@meta.data %>% rownames_to_column("barcodes"))

cor.test(test_stemness$Procr,test_stemness$`_stemness1`)
cor.test(test_stemness$Procr,test_stemness$`_stemness1`,method = "spearman")

# average expression ------------------------------------------------------
# build the grouping variable
# group_id_ENDO <- paste(scobj_ENDO@meta.data$Annotation.paper,
#                        scobj_ENDO@meta.data$Sample.paper,
#                        scobj_ENDO@meta.data$Phase,
#                        sep = "|")
group_id_ENDO <- paste(scobj_ENDO@meta.data$Annotation.paper,
                       scobj_ENDO@meta.data$Sample.paper,
                       sep = "|")
# add it to the meta
scobj_ENDO$group_id <- group_id_ENDO

# set the idents
Idents(scobj_ENDO) <- "group_id"

df_average_ENDO <- AverageExpression(scobj_ENDO,features = GOI) %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(names_to = "group_id",values_to = "avg_exp",-gene)

# build the pattern to extract the metadata
pattern_Annotation.paper <- paste0(unique(scobj_ENDO@meta.data$Annotation.paper),collapse = "|")
pattern_Sample.paper <- paste0(unique(scobj_ENDO@meta.data$Sample.paper),collapse = "|") %>% str_replace_all(pattern = "_",replacement = "\\.")
# pattern_Phase <- paste0(unique(scobj_ENDO@meta.data$Phase),collapse = "|")
# df_exp_avg <- df_average_ENDO %>%
#   mutate(Annotation.paper = str_extract_all(group_id,pattern = pattern_Annotation.paper) %>% unlist()) %>%
#   mutate(Sample.paper = str_extract_all(group_id,pattern = pattern_Sample.paper) %>% unlist()) %>%
#   mutate(Sample.paper = str_extract_all(group_id,pattern = pattern_Phase) %>% unlist())

df_exp_avg <- df_average_ENDO %>%
  mutate(Annotation.paper = str_extract_all(group_id,pattern = pattern_Annotation.paper) %>% unlist()) %>%
  mutate(Sample.paper = str_extract_all(group_id,pattern = pattern_Sample.paper) %>% unlist())

exp_avg_wide <- df_exp_avg %>%
  dplyr::select(group_id,gene,avg_exp) %>%
  pivot_wider(names_from = gene,values_from = avg_exp)

pdf("../../out/plot/ggpair_stemness_avg2.pdf",width = 15,height = 15)
ggpairs(exp_avg_wide[,-1])
dev.off()

# correlation matrix ------------------------------------------------------
test_correlation <- cor(exp_avg_wide %>%
                          column_to_rownames("group_id"))
# plot as heatmap
pdf("../../out/plot/correlation_stemness_pearson_matrix2.pdf",width = 5,height = 4)
Heatmap(test_correlation,
        name = "Correlation \nPearson",
        col = colorRamp2(c(-1, 0, 1), colors = c("blue", "white", "red")),
        row_names_side = "right",
        row_names_gp = gpar(fontsize = 8),
        column_names_side = "bottom",
        column_names_gp = gpar(fontsize = 8),
        row_dend_reorder = FALSE,
        column_dend_reorder = FALSE,
        row_title_gp = gpar(fontsize = 10, fontface = "bold"),
        column_title_gp = gpar(fontsize = 10, fontface = "bold"),
        show_column_names = F,
        show_row_names = T)
dev.off()

test_correlation2 <- cor(exp_avg_wide %>%
                           column_to_rownames("group_id"),
                         method = "spearman")
# plot as heatmap
pdf("../../out/plot/correlation_stemness_spearman_matrix2.pdf",width = 5,height = 4)
Heatmap(test_correlation2,
        name = "Correlation \nSpearman",
        col = colorRamp2(c(-1, 0, 1), colors = c("blue", "white", "red")),
        row_names_side = "right",
        row_names_gp = gpar(fontsize = 8),
        column_names_side = "bottom",
        column_names_gp = gpar(fontsize = 8),
        row_dend_reorder = FALSE,
        column_dend_reorder = FALSE,
        row_title_gp = gpar(fontsize = 10, fontface = "bold"),
        column_title_gp = gpar(fontsize = 10, fontface = "bold"),
        show_column_names = F,
        show_row_names = T)
dev.off()