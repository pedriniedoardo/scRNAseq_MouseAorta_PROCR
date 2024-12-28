data.combined <- readRDS("../../out/object/data.combined_harmonySkipIntegration_manualClean.rds")
test <- read.table("../../data/raf.txt",sep = "\t")

# load the final metadata provided by the authors
LUT <- read_csv("../../data/GSE211335_Global_seurat_metadata.csv")

dim(data.combined)
dim(LUT)
dim(test)

table(test$CellType)
table(test$cluster_high_res)
table(test$cluster_high_res)
table(test$RNA_snn_res.0.25)

# explore what cells are missing
original_meta <- data.combined@meta.data

df_test <- original_meta %>%
  left_join(test %>%
              rownames_to_column("barcode") %>%
              select(barcode,cluster_high_res,CellType),by = "barcode")

df_test %>%
  filter(is.na(CellType)) %>%
  group_by(Annotation.paper,Sample.paper,RNA_snn_res.0.1) %>%
  summarise(n = n())

df_test %>%
  filter(is.na(CellType)) %>%
  group_by(RNA_snn_res.0.1) %>%
  summarise(n = n())

df_test %>%
  filter(is.na(CellType)) %>%
  group_by(Annotation.paper,Sample.paper) %>%
  summarise(n = n())
