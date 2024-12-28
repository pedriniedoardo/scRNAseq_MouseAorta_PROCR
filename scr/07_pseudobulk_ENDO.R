# AIM ---------------------------------------------------------------------
# test the pseudobulk analysis on ENDO cells only, comparing different time points

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(GGally)
library(cowplot)
library(ComplexHeatmap)
library(scales)
library(circlize)
library(DESeq2)
library(RNAseqQC)
library(limma)
library(ashr)
library(magick)
library(UpSetR)

# read in the data --------------------------------------------------------
data.combined <- readRDS("../../out/object/data.combined_harmonySkipIntegration_manualClean.rds")
Idents(data.combined)<-"RNA_snn_res.0.1"

# subset only the presumed ENDO cells
scobj <- subset(data.combined,subset = RNA_snn_res.0.1 %in% c(0,9))

DimPlot(scobj,label = T,raster = T,group.by = "RNA_snn_res.0.1")

# wrangling ---------------------------------------------------------------
# explore the full meta
scobj@meta.data %>%
  group_by(Annotation.paper,Sample.paper) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = Annotation.paper,values_from = n)

# aggregate the expression per sample per treatment, donor and cell type
cts_sample_all <- AggregateExpression(object = scobj,
                                      group.by = c("Sample.paper","Annotation.paper"),
                                      assays = 'RNA',
                                      slot = "counts",
                                      return.seurat = FALSE)

# ASTRO processing --------------------------------------------------------
# focus on the ASTRO cluster of cells
# 1. Get counts matrix
counts_ENDO <- cts_sample_all$RNA %>%
  data.frame()
  # rownames_to_column("gene") %>%
  # # pull only the MG cells
  # dplyr::select("gene",contains("ASTRO")) %>%
  # # focus only of the MG from ctrl samples
  # dplyr::select("gene",contains("Ctrl")) %>%
  # # remove the unassigned donors
  # dplyr::select(-contains("unassigned")) %>%
  # column_to_rownames("gene")

# 2. generate sample level metadata
LUT_sample <- scobj@meta.data %>%
  group_by(Annotation.paper,Sample.paper) %>%
  summarise() %>%
  mutate(Sample.paper2 = str_replace_all(Sample.paper,pattern = "_",replace = "\\.")) %>%
  mutate(sample.id = paste(Sample.paper2,Annotation.paper,sep = "_"))

colData_ENDO <- data.frame(sample.id = colnames(counts_ENDO)) %>%
  # separate(samples,into = c(c("orig.ident","origin_01","pathology_class","origin_02")),remove = F,sep = "_") %>%
  left_join(LUT_sample,by = c("sample.id"))

# save matrix and metadata
saveRDS(counts_ENDO,file = "../../out/object/counts_ENDO_filterExp.rds")
saveRDS(colData_ENDO,file = "../../out/object/colData_ENDO_filterExp.rds")

# perform DESeq2 ----------------------------------------------------------
# build the model
treat_full <- colData_ENDO$Sample.paper
design_ENDO <- model.matrix(~ treat_full)
colnames(design_ENDO)[1] <- c("intercept")
# save the disign
saveRDS(design_ENDO,"../../out/object/design_ENDO_filterExp.rds")

# Create DESeq2 object
dds_ENDO <- DESeqDataSetFromMatrix(countData = counts_ENDO,
                                   colData = colData_ENDO,
                                   design = design_ENDO)

# filter
# filter at least 10 read per gene per sample
# keep_ASTRO <- rowSums(counts(dds_ENDO)) >=10
keep_ENDO <- edgeR::filterByExpr(counts(dds_ENDO), group = colData_ENDO$Sample.paper)
# keep_ASTRO <- rowSums(counts(dds_ENDO)) >=10
dds_ENDO_filter <- dds_ENDO[keep_ENDO,]

# plot the raw number of reads per sample
colSums(counts(dds_ENDO_filter)) %>%
  data.frame(tot_reads = .) %>%
  rownames_to_column("sample") %>%
  ggplot(aes(x=sample,y=tot_reads)) + geom_col()+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 45))

# scale the data
vds_ENDO_filter <- vst(dds_ENDO_filter, blind = F)

# clustering samples ------------------------------------------------------
# set seed to control random annotation colors
pdf("../../out/plot/heatmap_SampleCluster_ENDO_filterExp.pdf",width = 10,height = 6)
set.seed(1)
hm_ENDO <- plot_sample_clustering(vds_ENDO_filter,
                                  anno_vars = c("Sample.paper"),
                                  distance = "euclidean")
draw(hm_ENDO,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(2, 2, 2, 30), "mm"))
dev.off()

# PCA ---------------------------------------------------------------------
plot_vsd_ENDO <- plotPCA(vds_ENDO_filter,
                         intgroup = c("Sample.paper","Annotation.paper")) +
  theme_bw()

plot_vsd_ENDO$data %>%
  # mutate(sample = str_extract(name,pattern = "s\\d+")) %>%
  # ggplot(aes(x=PC1,y=PC2,col=condition,label=clone,shape=gender_update)) +
  # ggplot(aes(x=PC1,y=PC2,col=origin_01,label = sample)) +
  ggplot(aes(x=PC1,y=PC2,col=Sample.paper)) +
  geom_point(size =3,alpha=0.6) +
  # ggrepel::geom_text_repel(show.legend = F)+
  scale_x_continuous(expand = expansion(mult = 0.1))+
  theme_bw() + ylab(plot_vsd_ENDO$labels[1]) + xlab(plot_vsd_ENDO$labels[2])
ggsave("../../out/plot/PCA_pseudobulk_ENDO_filterExp.pdf",width = 6,height = 4)

# cell composition --------------------------------------------------------
# # martina suggested to check the relative proportion of different cell types from different sample to possibly justify the reason why this sample seems to be missclassified
# sample_prop_wide <- scobj@meta.data %>%
#   group_by(orig.ident,expertAnno.l1) %>%
#   summarise(n = n()) %>%
#   mutate(tot = sum(n)) %>%
#   mutate(prop = n/tot) %>%
#   # make it long
#   dplyr::select(orig.ident,expertAnno.l1,prop) %>%
#   pivot_wider(names_from = orig.ident,values_from = prop,values_fill = 0) %>%
#   column_to_rownames("expertAnno.l1") %>%
#   # select only the samples in the dataset
#   dplyr::select(colData_ENDO$orig.ident)
# 
# # plot the data as heatmap
# meta_sample_prop <- data.frame(colname = colnames(sample_prop_wide)) %>%
#   left_join(colData_ENDO,by=c("colname"="orig.ident"))
# 
# column_meta_sample_prop <- HeatmapAnnotation(treat = meta_sample_prop$origin_01,
#                                              col = list(treat = c("WM" = "blue",
#                                                                   "CX" = "yellow")))
# 
# ht2_shr_ENDO2 <- Heatmap(sample_prop_wide, show_column_names = T,raster_by_magick = T,show_row_dend = F, use_raster = T,
#                           name = "prop_cell_type",
#                           column_title = "sample",
#                           col = viridis::viridis(option = "turbo",n = 10),
#                           
#                           # row_names_gp = gpar(fontsize = 3),
#                           top_annotation = column_meta_sample_prop,show_row_names = T
#                           # cluster_rows = F,
#                           # right_annotation = row_ha,
#                           # row_split = rep(c(1,2,3,4),c(2,3,4,7))
#                           
# )
# # pdf("../../out/image/heatmap_cellType_pseudobulk_ASTRO_ctrl_refCX.pdf",width = 7,height = 5)
# draw(ht2_shr_ENDO2,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(30,2,2, 2), "mm"))
# # dev.off()
# 
# # check realitve proportion of the subcluster
# scobj_subset <- readRDS("../../out/object/ManualClean/data.combined_WM_CX_ASTRO_harmonySkipIntegration.rds")
# 
# DimPlot(scobj_subset,label = T,raster = T,group.by = "seurat_clusters")
# DimPlot(scobj_subset,label = T,raster = T,group.by = "expertAnno.l1")
# DimPlot(scobj_subset,label = T,raster = T,group.by = "pathology_class")
# 
# # realtive composition of the subclusters
# scobj_subset@meta.data %>%
#   filter(disease == "CTRL") %>%
#   group_by(orig.ident,seurat_clusters) %>%
#   summarise(n = n()) %>%
#   mutate(tot = sum(n)) %>%
#   mutate(prop = n/tot) %>%
#   ggplot(aes(x=orig.ident,y=prop,fill=seurat_clusters)) + geom_col()+theme_cowplot()
# # ggsave("../../out/image/barplot_cellType_pseudobulk_ASTRO_refCX_filterExp.pdf",width = 7,height = 5)
# 
# # render as an heatmap
# sample_prop_wide2 <- scobj_subset@meta.data %>%
#   filter(disease == "CTRL") %>%
#   group_by(orig.ident,seurat_clusters) %>%
#   summarise(n = n()) %>%
#   mutate(tot = sum(n)) %>%
#   mutate(prop = n/tot) %>%
#   # make it long
#   dplyr::select(orig.ident,seurat_clusters,prop) %>%
#   pivot_wider(names_from = orig.ident,values_from = prop,values_fill = 0) %>%
#   column_to_rownames("seurat_clusters")
# 
# # plot the data as heatmap
# meta_sample_prop2 <- data.frame(colname = colnames(sample_prop_wide)) %>%
#   left_join(colData_ENDO,by=c("colname"="orig.ident"))
# 
# column_meta_sample_prop2 <- HeatmapAnnotation(location = meta_sample_prop2$origin_01,
#                                               disease = meta_sample_prop2$disease,
#                                               col = list(location = c("WM" = "blue",
#                                                                       "CX" = "yellow"),
#                                                          disease = c("CTRL" = "green",
#                                                                      "MS" = "red")))
# 
# ht2_shr_MG22 <- Heatmap(sample_prop_wide2,
#                         show_column_names = T,
#                         raster_by_magick = T,
#                         show_row_dend = F, use_raster = T,
#                         name = "prop_cell_type",
#                         column_title = "sample",
#                         col = viridis::viridis(option = "turbo",n = 10),
#                         
#                         # row_names_gp = gpar(fontsize = 3),
#                         top_annotation = column_meta_sample_prop2,show_row_names = T
# )
# # pdf("../../out/image/heatmap_cellType_pseudobulk_ASTRO_refCX_filterExp.pdf",width = 7,height = 5)
# # draw(ht2_shr_MG22,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(30,2,2, 2), "mm"))
# # dev.off()


# run DE ------------------------------------------------------------------
# run DESeq2
ddsHTSeq_ENDO_filter <- DESeq(dds_ENDO_filter)

# Check the coefficients for the comparison
resultsNames(ddsHTSeq_ENDO_filter)

# save the filtered object
saveRDS(ddsHTSeq_ENDO_filter,"../../out/object/ddsHTSeq_pseudobulk_ENDO_filterExp.rds")

# print the contrast
resultsNames(ddsHTSeq_ENDO_filter)

# save the resut of the contrast of interest. notice that is also possible to set the alpha value for the significance (adj-p value)
contrast_ENDO <- makeContrasts(ENDO_Ctrl_vs_D3 = treat_fullDPT_3d,
                               ENDO_Ctrl_vs_D5 = treat_fullDPT_5d,
                               ENDO_Ctrl_vs_D7 = treat_fullDTP_7d,
                               levels = design_ENDO)

# pull the results table
res_ENDO.D3 <- results(ddsHTSeq_ENDO_filter, contrast=contrast_ENDO[,"ENDO_Ctrl_vs_D3"],alpha = 0.05)
res_ENDO.D5 <- results(ddsHTSeq_ENDO_filter, contrast=contrast_ENDO[,"ENDO_Ctrl_vs_D5"],alpha = 0.05)
res_ENDO.D7 <- results(ddsHTSeq_ENDO_filter, contrast=contrast_ENDO[,"ENDO_Ctrl_vs_D7"],alpha = 0.05)
summary(res_ENDO.D3)

# add the gene symbols
df_res <-
  list(res_ENDO.D3 = res_ENDO.D3,
       res_ENDO.D5 = res_ENDO.D5,
       res_ENDO.D7 = res_ENDO.D7) %>%
  lapply(function(x){
    df <- x %>%
      data.frame() %>%
      rownames_to_column(var = "symbol") %>%
      arrange(pvalue)
    # left_join(LUT_gene,by = c("symbol"="gene_name"))
    df
  }) %>%
  bind_rows(.id = "conditionVsCtrl")
# save the table
df_res %>%
  write_tsv("../../out/table/DE_pseudobulk_ENDO_filterExp.tsv")

# shrink ------------------------------------------------------------------
res_ENDO.D3_shr <- lfcShrink(ddsHTSeq_ENDO_filter, res = res_ENDO.D3, type = "ashr")
res_ENDO.D5_shr <- lfcShrink(ddsHTSeq_ENDO_filter, res = res_ENDO.D5, type = "ashr")
res_ENDO.D7_shr <- lfcShrink(ddsHTSeq_ENDO_filter, res = res_ENDO.D7, type = "ashr")

# save the table of DEGs
df_res_shr <-
  list(res_ENDO.D3_shr = res_ENDO.D3_shr,
       res_ENDO.D5_shr = res_ENDO.D5_shr,
       res_ENDO.D7_shr = res_ENDO.D7_shr) %>%
  lapply(function(x){
    df <- x %>%
      data.frame() %>%
      rownames_to_column(var = "symbol") %>%
      arrange(pvalue)
    # left_join(LUT_gene,by = c("symbol"="gene_name"))
    df
  }) %>%
  bind_rows(.id = "conditionVsCtrl")
# save the table
df_res_shr %>%
  write_tsv("../../out/table/DE_pseudobulk_ENDO_shr_filterExp.tsv")

# Another useful diagnostic plot is the histogram of the p values (figure below). This plot is best formed by excluding genes with very small counts, which otherwise generate spikes in the histogram.
df_res %>%
  data.frame()%>%
  dplyr::filter(baseMean>1) %>%
  ggplot(aes(x=pvalue))+geom_histogram(breaks = 0:20/20) +
  facet_wrap(~conditionVsCtrl)+
  theme_bw()+
  theme(strip.background = element_blank())
ggsave("../../out/plot/histogram_pvalue_pseudobulk_ENDO_filterExp.pdf",width = 6,height = 4)

# volcano -----------------------------------------------------------------
# add the info of the genename
plot_volcano <- df_res %>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1&!is.na(symbol)),yes = 1,no = 0)) %>%
  filter(!is.na(col))

plot_volcano %>%
  ggplot(aes(x=log2FoldChange,y=-log(padj)))+
  # geom_point()
  geom_point(data = plot_volcano[plot_volcano$col==0,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.05)+
  geom_point(data = plot_volcano[plot_volcano$col==1,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.5)+
  geom_vline(xintercept = c(-1,1),col="red",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="red",linetype="dashed")+
  scale_color_manual(values = c("black","red"))+theme(legend.position = "none")+
  # ggrepel::geom_text_repel(
  #   data = plot_volcano[plot_volcano$col==1,][1:1000,],
  #   aes(label = symbol),max.overlaps = 1,segment.alpha=0.4,
  #   size = 2,
  #   box.padding = unit(0.35, "lines"),
  #   point.padding = unit(0.3, "lines")) +
  ggrepel::geom_text_repel(
    data = plot_volcano %>% 
      group_by(conditionVsCtrl) %>% 
      arrange(padj) %>% 
      dplyr::slice(1:30) %>% 
      dplyr::filter(abs(log2FoldChange)>1) %>% 
      dplyr::filter(padj<0.05),
    aes(label = symbol),segment.alpha=0.4) +
  facet_wrap(~conditionVsCtrl)+
  theme_bw()+
  theme(strip.background = element_blank())+
  theme(legend.position = "none")
ggsave("../../out/plot/vulcano_plot_pseudobulk_ENDO_filterExp.pdf",width = 10,height = 10)

#
plot_volcano_shr <- df_res_shr %>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1&!is.na(symbol)),yes = 1,no = 0)) %>%
  filter(!is.na(col))

plot_volcano_shr %>%
  ggplot(aes(x=log2FoldChange,y=-log(padj)))+
  # geom_point()
  geom_point(data = plot_volcano_shr[plot_volcano_shr$col==0,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.05)+
  geom_point(data = plot_volcano_shr[plot_volcano_shr$col==1,],aes(x=log2FoldChange,y=-log(padj),col=factor(col)),alpha=0.5)+
  geom_vline(xintercept = c(-1,1),col="red",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="red",linetype="dashed")+
  scale_color_manual(values = c("black","red"))+theme(legend.position = "none")+
  ggrepel::geom_text_repel(
    data = plot_volcano_shr %>% group_by(conditionVsCtrl) %>% arrange(padj) %>% dplyr::slice(1:10)%>% filter(abs(log2FoldChange)>1),
    aes(label = symbol),segment.alpha=0.4) +
  facet_wrap(~conditionVsCtrl)+
  theme_bw()+
  theme(strip.background = element_blank())+
  theme(legend.position = "none")
ggsave("../../out/plot/vulcano_plot_pseudobulk_ENDO_shr_filterExp.pdf",width = 10,height = 10)

# MA plot -----------------------------------------------------------------
df_res %>%
  filter(!is.na(padj)) %>%
  data.frame()%>%
  mutate(color = ifelse(padj<0.05,yes = 1,no = 0))%>%
  mutate(color = factor(ifelse(is.na(color),yes = 0,no = color)))%>%
  arrange(color) %>%
  ggplot(aes(baseMean,y = log2FoldChange,col=color)) + geom_point(alpha=0.2) +
  scale_x_log10() + scale_color_manual(values = c("gray","red")) +
  facet_wrap(~conditionVsCtrl)+
  theme_bw()+
  theme(strip.background = element_blank()) +
  theme(legend.position = "none")
ggsave("../../out/plot/MA_plot_pseudobulk_ENDO_filterExp.pdf",width = 10,height = 10)

# shr
df_res_shr %>%
  filter(!is.na(padj)) %>%
  data.frame()%>%
  mutate(color = ifelse(padj<0.05,yes = 1,no = 0))%>%
  mutate(color = factor(ifelse(is.na(color),yes = 0,no = color)))%>%
  arrange(color) %>%
  ggplot(aes(baseMean,y = log2FoldChange,col=color)) + geom_point(alpha=0.2) +
  scale_x_log10() + scale_color_manual(values = c("gray","red")) +
  facet_wrap(~conditionVsCtrl)+
  theme_bw()+
  theme(strip.background = element_blank()) +
  theme(legend.position = "none")
ggsave("../../out/plot/MA_plot_pseudobulk_ENDO_shr_filterExp.pdf",width = 10,height = 10)

# heatmaps ----------------------------------------------------------------
# pull the scaled values
mat_filter_ASTRO <- assay(vds_ENDO_filter) %>%
  data.frame() %>%
  # dplyr::select(contains(c("_0_","_6_"))) %>%
  as.matrix()

# the DEGs plot
DEG_2 <- df_res_shr %>%
  data.frame()%>%
  # add a clor variable in case significant
  mutate(col=ifelse(((padj<0.05)&abs(log2FoldChange)>1),yes = 1,no = 0)) %>%
  dplyr::filter(col==1) %>%
  pull(symbol) %>%
  unique()

# mat_filter <- assay(vds_filter) %>%
#   data.frame() %>%
#   # dplyr::select(contains(c("_0_","_6_"))) %>%
#   as.matrix()

mat_shr_ENDO <- mat_filter_ASTRO[rownames(vds_ENDO_filter) %in% DEG_2, ]
mat2_shr_ENDO <- (mat_shr_ENDO - rowMeans(mat_shr_ENDO))/rowSds(mat_shr_ENDO,useNames = TRUE)
#
meta_sample_ENDO <- data.frame(colname = colnames(mat2_shr_ENDO)) %>%
  left_join(colData_ENDO,by=c("colname"="sample.id"))

# make the column of the matrix more readable
colnames(mat2_shr_ENDO) <- meta_sample_ENDO$colname

column_ha_shr_ENDO <- HeatmapAnnotation(treat = meta_sample_ENDO$Sample.paper,
                                         # gender = meta_sample_ENDO$sex,
                                         # facility = meta_sample_ENDO$facility,
                                         col = list(treat = c("Ctrl" = "green",
                                                              "DPT_3d" = "gray80",
                                                              "DPT_5d" = "gray30",
                                                              "DTP_7d" = "black")
                                                    # gender = c("M" = "cyan",
                                                    #            "F" = "pink"),
                                                    # facility = c("NIH" = "black",
                                                    #              "HSR" = "gray")
                                                    ))

ht2_shr_ENDO <- Heatmap(mat2_shr_ENDO, show_column_names = T,raster_by_magick = T,show_row_dend = F, use_raster = T,
                        name = "exp",
                        column_title = "ENDO_shr",
                        # row_names_gp = gpar(fontsize = 3),
                        top_annotation = column_ha_shr_ENDO,show_row_names = F,
                        # cluster_rows = F,
                        # right_annotation = row_ha,
                        # row_split = rep(c(1,2,3,4),c(2,3,4,7))
                         
)
pdf("../../out/plot/heatmap_DEG_plot_pseudobulk_ENDO_shr_filterExp.pdf",width = 7,height = 14)
draw(ht2_shr_ENDO,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(30,2,2, 2), "mm"))
dev.off()

# # upset plot --------------------------------------------------------------
# # read in the table of DEGs
# df_res_shr <- read_tsv("../../out/table/DE_pseudobulk_ASTRO_ctrl_refCX_shr.tsv")
# 
# # build a list of common elements belonging to each set of fegs
# list_DE_up <- df_res_shr %>%
#   split(f = .$conditionVsCX) %>%
#   lapply(function(x){
#     x %>%
#       filter(padj < 0.05,log2FoldChange>1) %>%
#       pull(symbol) %>%
#       unique()
#   })
# 
# glimpse(list_DE_up)
# 
# list_DE_down <- df_res_shr %>%
#   split(f = .$conditionVsCX) %>%
#   lapply(function(x){
#     x %>%
#       filter(padj < 0.05,log2FoldChange<(-1)) %>%
#       pull(symbol) %>%
#       unique()
#   })
# glimpse(list_DE_down)
# 
# # try the upset plot version
# # library(UpSetR)
# pdf("../../out/image/upset_DEG_UP_plot_pseudobulk_ASTRO_refBASELINE_shr.pdf",width = 14,height = 7)
# upset(fromList(list_DE_up), order.by = "freq",nsets = 7)
# dev.off()
# 
# pdf("../../out/image/upset_DEG_DOWN_plot_pseudobulk_ASTRO_refBASELINE_shr.pdf",width = 14,height = 7)
# upset(fromList(list_DE_down), order.by = "freq",nsets = 7)
# dev.off()
# 
# # pull the intersections
# df1_UP <- lapply(list_DE_up,function(x){
#   data.frame(gene = x)
# }) %>%
#   bind_rows(.id = "path")
# 
# df1_DOWN <- lapply(list_DE_down,function(x){
#   data.frame(gene = x)
# }) %>%
#   bind_rows(.id = "path")
# 
# head(df1_UP)
# head(df1_DOWN)
# 
# df2_UP <- data.frame(gene=unique(unlist(list_DE_up)))
# df2_DOWN <- data.frame(gene=unique(unlist(list_DE_down)))
# 
# head(df2_UP)
# head(df2_DOWN)
# 
# df_int_UP <- lapply(df2_UP$gene,function(x){
#   # pull the name of the intersections
#   intersection <- df1_UP %>%
#     dplyr::filter(gene==x) %>%
#     arrange(path) %>%
#     pull("path") %>%
#     paste0(collapse = "|")
# 
#   # build the dataframe
#   data.frame(gene = x,int = intersection)
# }) %>%
#   bind_rows()
# 
# df_int_DOWN <- lapply(df2_DOWN$gene,function(x){
#   # pull the name of the intersections
#   intersection <- df1_DOWN %>%
#     dplyr::filter(gene==x) %>%
#     arrange(path) %>%
#     pull("path") %>%
#     paste0(collapse = "|")
# 
#   # build the dataframe
#   data.frame(gene = x,int = intersection)
# }) %>%
#   bind_rows()
# 
# df_int_UP %>%
#   write_tsv("../../out/table/upset_DEG_UP_plot_pseudobulk_ASTRO_refBASELINE_shr.tsv")
# 
# df_int_DOWN %>%
#   write_tsv("../../out/table/upset_DEG_DOWN_plot_pseudobulk_ASTRO_refBASELINE_shr.tsv")
# 
# head(df_int_UP,n=20)
# head(df_int_DOWN,n=20)
# 
# df_int_UP %>%
#   group_by(int) %>%
#   summarise(n=n()) %>%
#   arrange(desc(n))
# 
# df_int_DOWN %>%
#   group_by(int) %>%
#   summarise(n=n()) %>%
#   arrange(desc(n))

# PLOT DISPERSION ---------------------------------------------------------
pdf("../../out/plot/dispersion_plot_pseudobulk_ENDO_filterExp.pdf",width = 5,height = 5)
plotDispEsts(ddsHTSeq_ENDO_filter)
dev.off()