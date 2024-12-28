# AIM ---------------------------------------------------------------------
# try to run harmony by merging the matrices from the individual objects. this will allow the skipping of the regular integration. the regular integration is needed as to run harmony the matrices should have the same number/order of the genes.
# I will be skippint he QC steps as I already have access the the final whitelist of the barcodes retained in the final object

# libraries ---------------------------------------------------------------
library("Seurat")
library("hdf5r")
library("tidyverse")
library("scater")
library("robustbase")
library("harmony")
library("ComplexHeatmap")
library("homologene")

# read in the data --------------------------------------------------------
# data have been downloaded from the geo repo:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116718
file_id <- c("GSM3259313_young_aorta.txt",
             "GSM3259314_old_aorta.txt")

# read them all
list_mat <- map(file_id,function(x){
  test <- read.table(paste0("../../data/",x),row.names = 1)
  return(test)
}) %>%
  setNames(str_extract(file_id,pattern = "old_aorta|young_aorta"))

# check the dimension of the matrices
lapply(list_mat,function(x){
  dim(x)
})

# check a portion of the matrices
lapply(list_mat,function(x){
  x[1:10,1:10]
})

# we do not have access to the metadata from the authors
# LUT <- read_csv("../../data/GSE211335_Global_seurat_metadata.csv")

# create a LUT for coverting the project name
LUT_sample <- data.frame(name = c("old_aorta","young_aorta"),
                         id = c("old_aorta","young_aorta"))

# read in the LUT to convert the gene names form mouse to human
homologeneData2_240528 <- read_tsv("../../data/Homologene_240528") %>%
  as.data.frame()

# wrangling ---------------------------------------------------------------
# notice that in this specifc case, I already have access to the list of filtered barcodes from each matrix, therefore I can skip the QC step
# create a seurat object from each raw matrix
# mat <- list_mat[[1]]
# mat_name <- names(list_mat)[1]
data.list <- pmap(list(list_mat,names(list_mat)),function(mat,mat_name){
  
  # define the final mat name
  mat_name_update <- LUT_sample %>%
    filter(name == mat_name) %>%
    pull(id)
  
  # for now retain all the cells. we are going to filter them out afterwards
  datasc <- CreateSeuratObject(counts = mat, project = mat_name_update, min.cells = 0, min.features = 0)
  
  # datasc <- NormalizeData(datasc) %>%
  #   FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  #   ScaleData() %>%
  #   RunPCA() %>% 
  #   FindNeighbors(dims = 1:30) %>%
  #   FindClusters() %>%
  #   RunUMAP(dims = 1:30) 
  
  # return the object
  return(datasc)
})

# merge the individual objetcs to create a single count matrix
data.list.id <- LUT_sample %>%
  dplyr::slice(match(names(list_mat),LUT_sample$name)) %>%
  pull(id)

# make a full merge of the tables
# NOTICE: make sure to use the list of the objects and not the list of matrices!
data.combined.all <- merge(data.list[[1]], y = data.list[-1], add.cell.ids = data.list.id, project = "MouseAortaProcr")

# check the size of the dataset
data.combined.all

# confirm the total number of cells
lapply(data.list,function(x){
  dim(x@assays$RNA$counts)[2]
}) %>%
  unlist() %>%
  sum()

# join the layers to be able to create a single object
# this step is needed for the new assay 5 structure
data.combined.all[["RNA"]] <- JoinLayers(data.combined.all[["RNA"]])

# notice it is critacal that all matrices have the same dimension
# I need to create a single object to add the cell cycle scoring and other metadata. I decided to keep all the cells and all the genes for now.
sobj_total <- CreateSeuratObject(counts = data.combined.all@assays$RNA$counts,
                                 project = "AortaProcr",
                                 meta.data = data.combined.all@meta.data,
                                 min.cells = 0, min.features = 0) %>%
  # this is needed as the cell cycle scoring is done on the data slot, which would be empty
  Seurat::NormalizeData(verbose = T)

# after creating the object I do not need the list and the merged matrix, free up some space
remove(data.list,data.combined.all)
gc()

# add the cell cycle analysis
DefaultAssay(sobj_total) <- "RNA"

# once updated pick the annotation of interest
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

s.genes.mouse <- homologene(s.genes, inTax = 9606, outTax = 10090,db = homologeneData2_240528) %>%
  pull(`10090`) %>%
  unique()
g2m.genes.mouse <- homologene(g2m.genes, inTax = 9606, outTax = 10090,db = homologeneData2_240528) %>%
  pull(`10090`) %>%
  unique()

sobj_total <- CellCycleScoring(sobj_total, s.features = s.genes.mouse, g2m.features = g2m.genes.mouse)
sobj_total$percent.mt <- PercentageFeatureSet(sobj_total, pattern = "^mt-")
sobj_total$percent.ribo <- PercentageFeatureSet(sobj_total, pattern = "^Rp[sl][[:digit:]]|^Rplp[[:digit:]]|^Rpsa")
# add also the percentage of globin. in this dataset it is not meaningful as there is no blood
sobj_total$percent.globin <- Seurat::PercentageFeatureSet(sobj_total,pattern = "^Hb[^(p)]")

# filter the whitelisted barcodes -----------------------------------------
# according to the authors the dataset should be filtered using the following paramters:
# Specifically, cells were first filtered to have at least 1000 UMIs, at least 100 genes and at most 10% mitochondrial gene expression



# read in the full draft dataset
# meta_whitelist <- LUT %>% 
#   # dplyr::rename("rowname" = "...1") %>%
#   select(barcode = "...1",
#          orig.ident.paper =  orig.ident,
#          nCount_RNA.paper = nCount_RNA,
#          nFeature_RNA.paper = nFeature_RNA,
#          percent.mito.paper = percent.mito,
#          Annotation.paper = Annotation,
#          Sample.paper = Sample,
#          RNA_snn_res.0.25.paper = RNA_snn_res.0.25,
#          seurat_clusters.paper = seurat_clusters)

# build the whitelist of the barcodes. in particular filter out all the barcodes not present in the final metadata shared by the authors

# # add all the original metadata
# meta_new <- sobj_total@meta.data %>%
#   rownames_to_column("rowname") %>%
#   separate(rowname,into = c("treat","tissue","barcode_id"),sep = "_",remove = F)

# # build the final meta
# meta_new_total <- meta_new %>%
#   left_join(meta_whitelist,by = c("barcode"),suffix=c(".new",".whitelist")) %>%
#   column_to_rownames("rowname") %>%
#   # build the whitelist condition
#   mutate(whitelist = case_when(is.na(orig.ident.paper) ~ 0,
#                                T ~ 1))
# 
# # check the dimension
# dim(meta_new_total)
# 
# # update the meta
# sobj_total@meta.data <- meta_new_total

# # verify the dimansions are correct
# dim(LUT)
# dim(meta_new_total)
# table(sobj_total@meta.data$whitelist)

# apply the filtering on the whitelist
# sobj_total <- subset(sobj_total,subset = whitelist == 1)

# confirm the dimenstions
# dim(sobj_total)

# integration processing --------------------------------------------------
# check the scale matrix
sobj_total@assays$RNA
# pull all the genes to scale
# all.genes <- rownames(sobj_total)

# rescale the data for regressing out the sources of variation do not scale all the genes. if needed I can scale them before the heatmap call. for speeding up the computation I will keep 
sobj_total <- sobj_total %>%
  # skip the normalizatio that has been already performed at the beginning
  # Seurat::NormalizeData(verbose = T) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 4000) %>%
  # I can scale the missing features afterwards now focus on the highly variable one for speed purposes
  ScaleData(vars.to.regress = c("percent.mt","nCount_RNA","S.Score","G2M.Score"), verbose = T) %>% 
  # run this if you want to scale all the variables
  # ScaleData(vars.to.regress = c("percent.mt.harmony","nCount_RNA.harmony","S.Score.harmony","G2M.Score.harmony"), verbose = T,features = all.genes) %>% 
  RunPCA(npcs = 30, verbose = T) %>% 
  RunUMAP(reduction = "pca", dims = 1:30,return.model = TRUE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30) %>%
  FindClusters(resolution = seq(0.1, 1, by = 0.1)) %>%
  identity()

# check the status of dataset preintegration
DimPlot(sobj_total,group.by = "orig.ident",raster = T)

# Run Harmony -------------------------------------------------------------
# The simplest way to run Harmony is to pass the Seurat object and specify which variable(s) to integrate out. RunHarmony returns a Seurat object, updated with the corrected Harmony coordinates. Let's set plot_convergence to TRUE, so we can make sure that the Harmony objective function gets better with each round.
sobj_total_h <- sobj_total %>%
  RunHarmony("orig.ident", plot_convergence = TRUE)
# Harmony with two or more covariates
# Do the same with your Seurat object:
# seuratObject <- RunHarmony(seuratObject, c("dataset", "donor", "batch_id"))
# To directly access the new Harmony embeddings, use the Embeddings command.
harmony_embeddings <- Embeddings(sobj_total_h, 'harmony')
harmony_embeddings[1:5, 1:5]
# Let's make sure that the datasets are well integrated in the first 2 dimensions after Harmony.
# DimPlot(object = sobj_total_h, reduction = "harmony", pt.size = .1, group.by = "sample_id")

# Downstream analysis -----------------------------------------------------
# Many downstream analyses are performed on low dimensional embeddings, not gene expression. To use the corrected Harmony embeddings rather than PCs, set reduction = 'harmony'. For example, let's perform the UMAP and Nearest Neighbor analyses using the Harmony embeddings.
sobj_total_h <- sobj_total_h %>%
  RunUMAP(reduction = "harmony", dims = 1:30,return.model = TRUE) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = seq(0.1, 1, by = 0.1)) %>%
  # FindClusters(resolution = 0.5) %>%
  identity()

# verify that all the relevant slots are filled
sobj_total_h@assays$RNA$counts[1:20,1:10]
sobj_total_h@assays$RNA$data[1:20,1:10]
sobj_total_h@assays$RNA$scale.data[1:20,1:10]

dim(sobj_total_h@assays$RNA$counts)
dim(sobj_total_h@assays$RNA$data)
dim(sobj_total_h@assays$RNA$scale.data)

DimPlot(sobj_total_h,group.by = "orig.ident",raster = F)
DimPlot(sobj_total_h,split.by = "orig.ident",raster = F,ncol = 4)

# -------------------------------------------------------------------------
# # add a costum annotation
# df_meta <- sobj_total_h@meta.data %>% 
#   data.frame() %>% 
#   mutate(orig_alt = case_when(orig.ident %in% c("s31","s27","s26","s25")~"wm_new",
#                               T ~ origin))
# 
# sobj_total_h$origin_alt <- df_meta$orig_alt
# -------------------------------------------------------------------------

# sobj_total_h <- readRDS("../../out/object/data.combined_harmonySkipIntegration_AllSoupX_01000_6000_15.rds")
saveRDS(sobj_total_h,"../../out/object/data.combined_harmonySkipIntegration.rds")
