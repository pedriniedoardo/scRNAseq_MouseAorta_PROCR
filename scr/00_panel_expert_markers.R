#Mouse ENDOTHELIAL SUBSET
c("Cd36", "Itga1", "Kit","Spcs2", "Calm3", "Smad6", "Tns1","Nrp1", "Bmpr2", "Zbtb20","Klf2","Rgs6","Casp3","Cdh5")

#Mouse Stromal Cell subset 

## Acta2+ Pericytes
FeaturePlot(seurat, features = c("Mustn1", "Des", "Acta2", "Myh11", "Tagln", "Tpm2"), order = T)
## Pericytes

FeaturePlot(seurat, features = c("Gucy1a1", "Pdzd2", "Gucy1b1", "Ndufa4l2", "Vtn", "Cspg4"), order = T)
FeaturePlot(seurat, features = c("Cspg4","Des","Rgs5","Pdgfrb","Mcam"), order = T, label = T)
VlnPlot(seurat, features = c("Cspg4","Des","Rgs5","Pdgfrb","Mcam"), ncol =5)
FeaturePlot(seurat, features = c("Pdgfra", "Pdgfrb"), order = T)
VlnPlot(seurat, features = c("Pdgfra", "Pdgfrb"))

## Fibroblasts

FeaturePlot(seurat, features = c("Col13a1", "Col14a1", "Col1a1"), order = T)
FeaturePlot(seurat, features = c("Mt1", "Mt2", "Hp", "Timp1", "Tagln"), order = T)

## Mesothelial

FeaturePlot(seurat, features = "Msln", order = T)

# Cluster Plots
FeaturePlot(seurat, features = c("Cspg4", "Pdzd2", "Acta2", "Des", "Tagln", "Pdgfrb", "Pdgfra", "Col13a1", "Col14a1",  "Msln"), order = T, ncol =5)

#Mouse Myeloid Subset

## Neutrophils
FeaturePlot(seurat, features = c("S100a9","S100a8","Stfa2","Retnlg"))

## Macrophages
FeaturePlot(seurat, features = "Cd68", order = T, label = T)

### Monocyte subtypes
FeaturePlot(seurat, features = c("Cd14", "Fcgr3", "Dab2", "Plac8"), order = T, label = T)

#Ly6c2 and Ccr2 for classical
#Fcgr4 for non-classical
FeaturePlot(seurat, features = c("Ly6c2", "Ccr2", "Fcgr4"), order = T)

#Interstitial Macs
FeaturePlot(seurat, feature = c( "C1qa", "Ms4a7", "Cd74", "H2-Aa"), label = T, order = T)

### Alveolar macrophages
FeaturePlot(seurat, feature = c("Bhlhe41", "Lpl", "Plet1", "Mrc1"), label = TRUE, order = TRUE)
FeaturePlot(seurat, feature = "Olr1")
FeaturePlot(seurat, feature = c("Itgax","Ca4","Marco","Csf2rb"), order = T)

## Dendritic Cells
FeaturePlot(seurat, feature = c("Siglech", "Flt3", "Zbtb46","Mreg"), label = F, order = TRUE)
FeaturePlot(seurat, features = c("Xcr1", "Mreg", "H2-Aa", "H2-Eb1"), order = T)
FeaturePlot(seurat, features = "Cd209a", order = T)
FeaturePlot(seurat, features = c("Ccl17", "Cyb561a3"), order = T)
FeaturePlot(seurat, features = c("Itgam","Itgae","Fcgr1","Itgax"), order = T)
FeaturePlot(seurat, features=c("Batf3","Clec10a", "Cst3"),label=T, order=T)
FeaturePlot(seurat, features = c("Epsti1", "Fscn1","Ccr7","Traf1"), order = T)

## Mast or Basophils
FeaturePlot(seurat, feature = c("Cpa3", "Mcpt4", "Mcpt8l3", "Cma1"), label = TRUE, order = TRUE)
FeaturePlot(seurat, features = c("Ccl3", "Ccl4", "Il6", "Ifitm1"), order = T)

# FeaturePlots for myeloid discrimination
FeaturePlot(seurat, features = c("Marco", "S100a9", "Fcgr3", "Cd14", "Ly6c2", "Fcgr4","H2-Aa", "C1qa", "Mreg", "Cd209a", "Itgae", "Cpa3"), ncol = 6, order=T)

#Mouse Lymphoid Subset

FeaturePlot(seurat, features = c("Cd3d", "Nkg7", "Cd79b"), label =T, order=T, ncol=3)

FeaturePlot(seurat, features = c("Ptprc", "Cldn5"), order = T)

## Cluster 9 - Stromal (fibroblast/ SMC) - ECM secreting cells
FeaturePlot(seurat, features = c("Mgp", "Bgn", "Ogn", "Fn1"), order = T)

## Cluster 10 - platelets
FeaturePlot(seurat, features = c("Ppbp", "Pf4"))

# B cell markers
FeaturePlot(seurat, features = c("Cd79b","Cd79a","Ms4a1","Cd74"), label = T)
FeaturePlot(seurat, features = c("Plac8", "Vpreb3", "Iglc1"), order = T)


# NK Cells
FeaturePlot(seurat, features = c("Nkg7","Prf1","Gzma","Gzmm"), label = T)
FeaturePlot(seurat, features = "Xcl1")
FeaturePlot(seurat, features = c("Klrd1", "Trgc2"))

# T Cells
FeaturePlot(seurat, features = c("Cd3d", "Cd8a", "Cd4"), order = T, label = T)
VlnPlot(seurat, features = c("Cd8a", "Cd4"), split.by = "Sample")
FeaturePlot(seurat, features = c("S100a4","Cxcr6"))
FeaturePlot(seurat, features = c("Cd4","Tcf7","Lef1","Dusp10"), order = T, label = T)

FeaturePlot(seurat, features = c("Foxp3","Il2ra", "Itgb1", "Icos"), order = T)
FeaturePlot(seurat, features = c("Ctla4", "S100a6"))
FeaturePlot(seurat, feature = "Dusp2", order =T)
FeaturePlot(seurat, features = c("Rag1", "Dntt"), order = T)


FeaturePlot(seurat, features = c("Nkg7", "Gzma", "Cd3d", "Cd4", "Cd8a","Cxcr6", "Il2ra", "Dntt", "Cd79a", "Cd74"), order=T, ncol = 5)
