#载入包
if(!require(multtest))BiocManager::install("multtest")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(patchwork))install.packages("patchwork")
if(!require(R.utils))install.packages("R.utils")
if(!require(ggplot2))BiocManager::install("ggplot2")
if(!require(clustree))install.packages("clustree")
if(!require(cowplot))install.packages("cowplot")
if(!require(DropletUtils))install.packages("DropletUtils")
if(!require(DoubletFinder))install.packages("DoubletFinder")
if(!require(tidyverse))install.packages("tidyverse")
if(!require(Rcpp))install.packages("Rcpp")
if(!require(pheatmap))install.packages("pheatmap")
if(!require(ggsignif)) install.packages("ggsignif")
if(!require(gridExtra)) install.packages("gridExtra")
if(!require(presto)) install.packages("presto")
if(!require(rvcheck)) install.packages("rvcheck", repos = "http://cran.us.r-project.org")
suppressWarnings(suppressMessages(if(!require(clusterProfiler))BiocManager::install("clusterProfiler")))
if(!require(clusterProfiler)) install.packages("clusterProfiler")
if(!require(AnnotationDbi)) install.packages("AnnotationDbi")
suppressWarnings(suppressMessages(if(!require(org.Mm.eg.db))BiocManager::install("org.Mm.eg.db")))
suppressWarnings(suppressMessages(if(!require(org.Hs.eg.db))BiocManager::install("org.Hs.eg.db")))
suppressWarnings(suppressMessages(library(dplyr)))
library(dplyr)
library(clusterProfiler)
library(gridExtra)
library(ggsignif)
library(argparse)
library(SeuratWrappers)
library(clustree)
library(patchwork)
library(AnnoProbe)
library(infercnv)
library(grid)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggrepel)
rm(list = ls());gc()
setwd("C:/Users/wangyanhong/Desktop/LNLC")
scRNAseq.sct.int <- readRDS("new/LNLC_FAM.rds")
B_cell <- readRDS('B_cell.rds')
Endothelial_cell <- readRDS('Endothelial_cell.rds')
Epithelial_cell <- readRDS('Epithelial_cell.rds')
Fibroblast <- readRDS('Fibroblast.rds')
Mast_cell <- readRDS('Mast_cell.rds')
Myeloid_cell <- readRDS('Myeloid_cell.rds')
NK_cell <- readRDS('NK_cell.rds')
T_cell <- readRDS('T_cell.rds')
setwd("C:/Users/wangyanhong/Desktop/LNLC/figure")
set.seed(123)
#定义细胞亚群
complete_celltype_list <- data.frame(
  ClusterID = 0:15,  # 假设 ClusterID 从 0 到 12
  celltype = c("T_cell","T_cell",  "Myeloid_cell", "Epithelial_cell","NK_cell", "Mast_cell","B_cell","Myeloid_cell","Fibroblast",
                     "Endothelial_cell", "T_cell","B_cell","Epithelial_cell","Fibroblast", "Myeloid_cell","Myeloid_cell"))
# 更新 celltype 数据框架
celltype <- complete_celltype_list
head(celltype)
celltype
table(celltype$celltype)
scRNAseq.sct.int@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNAseq.sct.int@meta.data[which(scRNAseq.sct.int@meta.data$RNA_snn_res.0.1 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(scRNAseq.sct.int@meta.data$celltype)
selected_genes=  c("CD3D", "CD3G", "TRAC",
                   "EPCAM", "MUC1", "CDH1", "EPCAM",
                   "LYZ", "C1QA", "C1QB",
                   "CD79A","IGHA1","IGKC",
                   "KIT", "MS4A2",
                   "GNLY", "KLRD1",
                   "ACTA2", "COL1A1", "COL3A1",
                   "VWF", "CDH5","CLDN5")
DotPlot(scRNAseq.sct.int, features = unique(selected_genes)) + RotatedAxis()
ggsave("DotPlot_selected_genes_20240516.pdf", width = 10, height = 5)
# # 定义每个亚群的分辨率和降维维度
params <- list(
  Myeloid_cell = list(resolution = 0.10, dims = 1:40),
  Endothelial_cell = list(resolution = 0.40, dims = 1:40),
  Mast_cell = list(resolution = 0.20, dims = 1:40),
  Fibroblast = list(resolution = 0.30, dims = 1:30),
  B_cell = list(resolution = 0.30, dims = 1:30),
  T_cell = list(resolution = 0.20, dims = 1:40),
  Epithelial_cell = list(resolution = 0.30, dims = 1:30),
  NK_cell = list(resolution = 0.10, dims = 1:40)
)
# 如果提取后直接从这部运行可以保留大群分群# 
Endothelial_cell<- FindNeighbors(Endothelial_cell, reduction = "integrated.cca",dims = 1:40) 
Endothelial_cell<- FindClusters(Endothelial_cell, resolution = 0.40, cluster.name = "cca_clusters")
# UMAP的结果
Endothelial_cell <- RunUMAP(Endothelial_cell, reduction = "integrated.cca",  dims = 1:40, 
                            reduction.name = "umap.cca")
DimPlot(Endothelial_cell, group.by = "cca_clusters", label = T)
ggsave("Endothelial_cell/umap_Endothelial_cell_Tumor_after_integrated.cca.pdf", width = 10, height = 8)
table(Endothelial_cell@meta.data$seurat_clusters)
table(Endothelial_cell@meta.data$orig.ident)
# 定义内皮细胞亚群
complete_celltype_list <- data.frame(
  ClusterID = 0:7,  # 假设 ClusterID 从 0 到 7
  celltype = c("SELE_vascular_ECs_C0", "SLC6A4_vascular_ECs_C1", "INSR_vascular_ECs_C2", 
               "MT1X_ECs_C3", "ING1_ECs_C4", "PTGIS_vascular_ECs_C5", 
               "MGAT4C_lymphatic_ECs_C6", "PRR16_capillary-arterial_ECs_C7")
)
# 更新 celltype 列为 "NA"
Endothelial_cell@meta.data$celltype <- "NA"
# 确保 seurat_clusters 是字符型
Endothelial_cell@meta.data$seurat_clusters <- as.character(Endothelial_cell@meta.data$seurat_clusters)
# 更新 celltype 列
for (i in 1:nrow(complete_celltype_list)) {
  cluster_id <- as.character(complete_celltype_list$ClusterID[i])  # 确保 cluster_id 为字符型
  cell_type <- complete_celltype_list$celltype[i]
  
  matching_cells <- which(Endothelial_cell@meta.data$seurat_clusters == cluster_id)
  if (length(matching_cells) > 0) {
    Endothelial_cell@meta.data[matching_cells, 'celltype'] <- cell_type
  }
}
# 检查更新结果
print(table(Endothelial_cell@meta.data$celltype))
selected_gene=   c("PECAM1", "CDH5", 
                   "PROX1", "PDPN",
                   "TGFB2", "GLUL") 
# 使用 DotPlot 函数绘制图形，按照更新后的细胞类型进行分组
dot_plot <- DotPlot(Endothelial_cell, features = unique(selected_gene), group.by = "celltype") + RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))  # 调整 x 轴标签的显示
# 保存 DotPlot 图
ggsave("Endothelial_cell/DotPlot_Endothelial_cell0519.pdf", plot = dot_plot, width = 10, height = 8)
# 提取内皮细胞的注释信息
endothelial_metadata <- Endothelial_cell@meta.data[, c("celltype", "seurat_clusters")]
# 确保 scRNAseq.sct.int 和 Endothelial_cell 的细胞名匹配
common_cells <- intersect(rownames(scRNAseq.sct.int@meta.data), rownames(endothelial_metadata))
length(common_cells)  # 检查共有细胞数
# 检查内皮细胞的总细胞数
total_cells <- sum(table(Endothelial_cell@meta.data$seurat_clusters))
cat("内皮细胞总数：", total_cells, "\n")

# 验证 `scRNAseq.sct.int` 中是否包含所有的内皮细胞
common_cells <- intersect(rownames(scRNAseq.sct.int@meta.data), rownames(Endothelial_cell@meta.data))
num_common_cells <- length(common_cells)
cat("共有细胞数：", num_common_cells, "\n")

# 检查共有细胞是否匹配
if (total_cells != num_common_cells) {
  cat("警告：内皮细胞总数与共有细胞数不匹配！\n")
} else {
  cat("内皮细胞总数与共有细胞数匹配。\n")
}
# 创建一个新的列以存储内皮细胞的注释信息
scRNAseq.sct.int@meta.data$sub_celltype <- "NA"
# 将注释信息添加到原始数据对象中
scRNAseq.sct.int@meta.data[common_cells, "sub_celltype"] <- endothelial_metadata[common_cells, "celltype"]
# 检查结果
print(table(scRNAseq.sct.int@meta.data$sub_celltype))
####髓系亚群####
# 如果提取后直接从这部运行可以保留大群分群# 
Myeloid_cell <- FindNeighbors(Myeloid_cell, reduction = "integrated.cca", dims = 1:40) 
Myeloid_cell <- FindClusters(Myeloid_cell, resolution = 0.10, cluster.name = "cca_clusters")
# UMAP的结果
Myeloid_cell <- RunUMAP(Myeloid_cell, reduction = "integrated.cca", dims = 1:40, 
                        reduction.name = "umap.cca")
DimPlot(Myeloid_cell, group.by = "cca_clusters", label = TRUE)
ggsave("Myeloid_cell/umap_Myeloid_cell_Tumor_after_integrated.cca.pdf", width = 10, height = 8)
print(table(Myeloid_cell@meta.data$seurat_clusters))
print(table(Myeloid_cell@meta.data$orig.ident))
# 定义髓系细胞亚群
complete_celltype_list <- data.frame(
  ClusterID = 0:7,  # 假设 ClusterID 从 0 到 7
  celltype = c("FABP4_Macrophages_C0", "SELENOP_Macrophages_C1", "CD1E_cDC2s_C2", "FCGR3B_Monocytes_C3", "SPC25_Macrophages_C4", 
               "CHRDL1_Mast_cells_C5", "SCT_pDCs_C6", "SFTPD_Mast_cells_C7")
)
# 更新 celltype 列为 "NA"
Myeloid_cell@meta.data$celltype <- "NA"
# 确保 seurat_clusters 是字符型
Myeloid_cell@meta.data$seurat_clusters <- as.character(Myeloid_cell@meta.data$seurat_clusters)
# 更新 celltype 列
for (i in 1:nrow(complete_celltype_list)) {
  cluster_id <- as.character(complete_celltype_list$ClusterID[i])  # 确保 cluster_id 为字符型
  cell_type <- complete_celltype_list$celltype[i]
  
  matching_cells <- which(Myeloid_cell@meta.data$seurat_clusters == cluster_id)
  if (length(matching_cells) > 0) {
    Myeloid_cell@meta.data[matching_cells, 'celltype'] <- cell_type
  }
}
# 绘制带有注释的 UMAP 图
umap_plot <- DimPlot(Myeloid_cell, reduction = "umap", group.by = "celltype", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP Plot of Myeloid Cells") +
  theme(plot.title = element_text(hjust = 0.5))  # 标题居中
# 保存 UMAP 图
ggsave("Myeloid_cell/UMAP_Myeloid_cell_annotated_0519.pdf", plot = umap_plot, width = 10, height = 8)

# 提取髓系细胞的注释信息
myeloid_metadata <- Myeloid_cell@meta.data[, c("celltype", "seurat_clusters")]
# 确保 scRNAseq.sct.int 和 Myeloid_cell 的细胞名匹配
common_cells <- intersect(rownames(scRNAseq.sct.int@meta.data), rownames(myeloid_metadata))
length(common_cells)  # 检查共有细胞数
# 将注释信息添加到原始数据对象中
scRNAseq.sct.int@meta.data[common_cells, "sub_celltype"] <- myeloid_metadata[common_cells, "celltype"]
# 检查结果
print(table(scRNAseq.sct.int@meta.data$sub_celltype))

###### Mast_cell 亚群分析 Mast_cell = list(resolution = 0.20, dims = 1:40
Mast_cell <- FindNeighbors(Mast_cell, reduction = "integrated.cca", dims = 1:40)
Mast_cell <- FindClusters(Mast_cell, resolution = 0.20, cluster.name = "cca_clusters")
Mast_cell <- RunUMAP(Mast_cell, reduction = "integrated.cca", dims = 1:40, reduction.name = "umap.cca")
DimPlot(Mast_cell, group.by = "cca_clusters", label = TRUE)
ggsave("C:/Users/wangyanhong/Desktop/LNLC/figure/umap_Mast_cell_Tumor_after_integrated.cca.pdf", width = 10, height = 8)
print(table(Mast_cell@meta.data$seurat_clusters))
print(table(Mast_cell@meta.data$orig.ident))
# 定义肥大细胞亚群
complete_celltype_list <- data.frame(
  ClusterID = 0:2,  # 假设 ClusterID 从 0 到 2
  celltype = c("NTRK1_Mast_cells_C0", "CTSG_Mast_cells_C1", "GNLY_Mast_cells_C2")
)
# 更新 celltype 列为 "NA"
Mast_cell@meta.data$celltype <- "NA"
# 确保 seurat_clusters 是字符型
Mast_cell@meta.data$seurat_clusters <- as.character(Mast_cell@meta.data$seurat_clusters)
# 更新 celltype 列
for (i in 1:nrow(complete_celltype_list)) {
  cluster_id <- as.character(complete_celltype_list$ClusterID[i])  # 确保 cluster_id 为字符型
  cell_type <- complete_celltype_list$celltype[i]
  
  matching_cells <- which(Mast_cell@meta.data$seurat_clusters == cluster_id)
  if (length(matching_cells) > 0) {
    Mast_cell@meta.data[matching_cells, 'celltype'] <- cell_type
  }
}
# 肥大细胞亚群的 marker gene 列表
marker_genes <- c("CPA3", "MALAT1", "GNLY")
# 绘制 DotPlot 图
dot_plot <- DotPlot(Mast_cell, features = unique(marker_genes), group.by = "celltype") + RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))  # 调整 x 轴标签的显示
# 保存 DotPlot 图
ggsave("Mast_cell/DotPlot_Mast_cell.pdf", plot = dot_plot, width = 10, height = 8)
# 绘制带有注释的 UMAP 图
umap_plot <- DimPlot(Mast_cell, reduction = "umap", group.by = "celltype", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP Plot of Mast Cells") +
  theme(plot.title = element_text(hjust = 0.5))  # 标题居中
# 保存 UMAP 图
ggsave("Mast_cell/UMAP_Mast_cell_annotated.pdf", plot = umap_plot, width = 10, height = 8)
# 提取肥大细胞的注释信息
Mast_cell_metadata <- Mast_cell@meta.data[, c("celltype", "seurat_clusters")]
# 确保 scRNAseq.sct.int 和 Mast_cell 的细胞名匹配
common_cells <- intersect(rownames(scRNAseq.sct.int@meta.data), rownames(Mast_cell_metadata))
length(common_cells)  # 检查共有细胞数
# 将注释信息添加到原始数据对象中
scRNAseq.sct.int@meta.data[common_cells, "sub_celltype"] <- Mast_cell_metadata[common_cells, "celltype"]
# 检查结果
print(table(scRNAseq.sct.int@meta.data$sub_celltype))

# Fibroblast 亚群分析 Fibroblast = list(resolution = 0.30, dims = 1:30
Fibroblast <- FindNeighbors(Fibroblast, reduction = "integrated.cca", dims = 1:30)
Fibroblast <- FindClusters(Fibroblast, resolution = 0.30, cluster.name = "cca_clusters")
Fibroblast <- RunUMAP(Fibroblast, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
DimPlot(Fibroblast, group.by = "cca_clusters", label = TRUE)
ggsave("C:/Users/wangyanhong/Desktop/LNLC/figure/umap_Fibroblast_Tumor_after_integrated.cca.pdf", width = 10, height = 8)
print(table(Fibroblast@meta.data$seurat_clusters))
print(table(Fibroblast@meta.data$orig.ident))
# 定义成纤维细胞亚群
complete_celltype_list <- data.frame(
  ClusterID = 0:7,  # 假设 ClusterID 从 0 到 4
  celltype = c("PI16_Inflammatory_CAFs_C0","LAMC3_Myofibrotic_CAFs_C1", "MMP11_Inflammatory_CAFs_C2", "CD74_Inflammatory_CAFs_C3", "FGFR4_Inflammatory_CAFs_C4", 
               "DES_Myofibrotic_CAFs_C5", "CPA3_Antigen-presenting_CAFs_C6", "SILC1_EMT-like_CAFs_C7")
)
# 更新 celltype 列为 "NA"
Fibroblast@meta.data$celltype <- "NA"
# 确保 seurat_clusters 是字符型
Fibroblast@meta.data$seurat_clusters <- as.character(Fibroblast@meta.data$seurat_clusters)
# 更新 celltype 列
for (i in 1:nrow(complete_celltype_list)) {
  cluster_id <- as.character(complete_celltype_list$ClusterID[i])  # 确保 cluster_id 为字符型
  cell_type <- complete_celltype_list$celltype[i]
  
  matching_cells <- which(Fibroblast@meta.data$seurat_clusters == cluster_id)
  if (length(matching_cells) > 0) {
    Fibroblast@meta.data[matching_cells, 'celltype'] <- cell_type
  }
}
# 髓系亚群的 marker gene 列表
marker_genes <- c("CXCL14", "PDGFRA", "MEF2C", "COL4A1", "C1QA", "C1QB", "SILC1")
# 绘制 DotPlot 图
dot_plot <- DotPlot(Fibroblast, features = unique(marker_genes), group.by = "cca_clusters") + RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))  # 调整 x 轴标签的显示
# 保存 DotPlot 图
ggsave("Fibroblast/DotPlot_Fibroblast_cell.pdf", plot = dot_plot, width = 10, height = 8)
# 绘制带有注释的 UMAP 图
umap_plot <- DimPlot(Fibroblast, reduction = "umap", group.by = "celltype", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP Plot of Fibroblast Cells") +
  theme(plot.title = element_text(hjust = 0.5))  # 标题居中
# 保存 UMAP 图
ggsave("Fibroblast/UMAP_Fibroblast_cell_annotated.pdf", plot = umap_plot, width = 10, height = 8)
# 提取成纤维细胞的注释信息
Fibroblast_metadata <- Fibroblast@meta.data[, c("celltype", "seurat_clusters")]
# 确保 scRNAseq.sct.int 和 Fibroblast 的细胞名匹配
common_cells <- intersect(rownames(scRNAseq.sct.int@meta.data), rownames(Fibroblast_metadata))
length(common_cells)  # 检查共有细胞数
# 将注释信息添加到原始数据对象中
scRNAseq.sct.int@meta.data[common_cells, "sub_celltype"] <- Fibroblast_metadata[common_cells, "celltype"]
# 检查结果
print(table(scRNAseq.sct.int@meta.data$sub_celltype))

# B 细胞亚群分析 B_cell = list(resolution = 0.30, dims = 1:30)
B_cell <- FindNeighbors(B_cell, reduction = "integrated.cca", dims = 1:30)
B_cell <- FindClusters(B_cell, resolution = 0.30, cluster.name = "cca_clusters")
B_cell <- RunUMAP(B_cell, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
DimPlot(B_cell, group.by = "cca_clusters", label = TRUE)
ggsave("C:/Users/wangyanhong/Desktop/LNLC/figure/umap_B_cell_Tumor_after_integrated.cca.pdf", width = 10, height = 8)
print(table(B_cell@meta.data$seurat_clusters))
print(table(B_cell@meta.data$orig.ident))
# 定义 B 细胞亚群
complete_celltype_list <- data.frame(
  ClusterID = 0:5,  # 假设 ClusterID 从 0 到 5
  celltype = c("NIBAN3_naive_B_cells_C0","MT2A_naive_B_cells_C1", "DUSP4_naive_B_cells_C2", "CHAC1_Plasma_B_cells_C3", "IGLV6-57_Plasma_B_cells_C4", 
               "IGLC2_Plasma_B_cells_C5")
)
# 更新 celltype 列为 "NA"
B_cell@meta.data$celltype <- "NA"
# 确保 seurat_clusters 是字符型
B_cell@meta.data$seurat_clusters <- as.character(B_cell@meta.data$seurat_clusters)
# 更新 celltype 列
for (i in 1:nrow(complete_celltype_list)) {
  cluster_id <- as.character(complete_celltype_list$ClusterID[i])  # 确保 cluster_id 为字符型
  cell_type <- complete_celltype_list$celltype[i]
  
  matching_cells <- which(B_cell@meta.data$seurat_clusters == cluster_id)
  if (length(matching_cells) > 0) {
    B_cell@meta.data[matching_cells, 'celltype'] <- cell_type
  }
}
# B 细胞亚群的 marker gene 列表
marker_genes <- c("CD69", "VPREB3", "CD52", 'CD38', 'XBP1', 'JCHAIN')
# 绘制 DotPlot 图
dot_plot <- DotPlot(B_cell, features = unique(marker_genes), group.by = "celltype") + RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))  # 调整 x 轴标签的显示
# 保存 DotPlot 图
ggsave("B_cell/DotPlot_B_cell.pdf", plot = dot_plot, width = 10, height = 8)
# 绘制带有注释的 UMAP 图
umap_plot <- DimPlot(B_cell, reduction = "umap", group.by = "celltype", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP Plot of B Cells") +
  theme(plot.title = element_text(hjust = 0.5))  # 标题居中
# 保存 UMAP 图
ggsave("B_cell/UMAP_B_cell_annotated.pdf", plot = umap_plot, width = 10, height = 8)
# 提取 B 细胞的注释信息
B_cell_metadata <- B_cell@meta.data[, c("celltype", "seurat_clusters")]
# 确保 scRNAseq.sct.int 和 B_cell 的细胞名匹配
common_cells <- intersect(rownames(scRNAseq.sct.int@meta.data), rownames(B_cell_metadata))
length(common_cells)  # 检查共有细胞数
# 将注释信息添加到原始数据对象中
scRNAseq.sct.int@meta.data[common_cells, "sub_celltype"] <- B_cell_metadata[common_cells, "celltype"]
# 检查结果
print(table(scRNAseq.sct.int@meta.data$sub_celltype))

# T 细胞亚群分析T_cell = list(resolution = 0.20, dims = 1:40)
T_cell <- FindNeighbors(T_cell, reduction = "integrated.cca", dims = 1:40)
T_cell <- FindClusters(T_cell, resolution = 0.20, cluster.name = "cca_clusters")
T_cell <- RunUMAP(T_cell, reduction = "integrated.cca", dims = 1:40, reduction.name = "umap.cca")
DimPlot(T_cell, group.by = "cca_clusters", label = TRUE)
ggsave("C:/Users/wangyanhong/Desktop/LNLC/figure/umap_T_cell_Tumor_after_integrated.cca.pdf", width = 10, height = 8)
print(table(T_cell@meta.data$seurat_clusters))
print(table(T_cell@meta.data$orig.ident))
# 定义 T 细胞亚群
complete_celltype_list <- data.frame(
  ClusterID = 0:4,  # 假设 ClusterID 从 0 到 7
  celltype = c("CD40LG_Naive_T_Cells_C0", "ENC1_CD8_Exhausted_T_cells_C1", "KLRC1_CD8_Follicular_Helper_T_Cells_C2", 
               "FOXP3_T_Regulatory_Cells_C3", "NEAT1_MTgene_T_Cells_C4")
)
# 更新 celltype 列为 "NA"
T_cell@meta.data$celltype <- "NA"
# 确保 seurat_clusters 是字符型
T_cell@meta.data$seurat_clusters <- as.character(T_cell@meta.data$seurat_clusters)
# 更新 celltype 列
for (i in 1:nrow(complete_celltype_list)) {
  cluster_id <- as.character(complete_celltype_list$ClusterID[i])  # 确保 cluster_id 为字符型
  cell_type <- complete_celltype_list$celltype[i]
  
  matching_cells <- which(T_cell@meta.data$seurat_clusters == cluster_id)
  if (length(matching_cells) > 0) {
    T_cell@meta.data[matching_cells, 'celltype'] <- cell_type
  }
}
# T 细胞的 marker gene 列表
marker_genes_t <- c("CD4","CD8A","CD8B", "TCF7", "CCR7","PDCD1","TOX", "HOPX", "ZNF683","IL2RA", "FOXP3", "IKZF2")
# 绘制 T 细胞的 DotPlot 图
dot_plot_t <- DotPlot(T_cell, features = unique(marker_genes_t), group.by = "celltype") + RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))  # 调整 x 轴标签的显示
# 保存 T 细胞 DotPlot 图
ggsave("T_cell/DotPlot_T_cell23.pdf", plot = dot_plot_t, width = 10, height = 8)
# 绘制带有注释的 UMAP 图
umap_plot <- DimPlot(T_cell, reduction = "umap", group.by = "celltype", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP Plot of T Cells") +
  theme(plot.title = element_text(hjust = 0.5))  # 标题居中

# 保存 UMAP 图
ggsave("T_cell/UMAP_T_cell_annotated.pdf", plot = umap_plot, width = 10, height = 8)
# 提取 T 细胞的注释信息
T_cell_metadata <- T_cell@meta.data[, c("celltype", "seurat_clusters")]
# 确保 scRNAseq.sct.int 和 T_cell 的细胞名匹配
common_cells <- intersect(rownames(scRNAseq.sct.int@meta.data), rownames(T_cell_metadata))
length(common_cells)  # 检查共有细胞数
# 将注释信息添加到原始数据对象中
scRNAseq.sct.int@meta.data[common_cells, "sub_celltype"] <- T_cell_metadata[common_cells, "celltype"]
# 检查结果
print(table(scRNAseq.sct.int@meta.data$sub_celltype))            
######上皮细胞亚群分析######Epithelial_cell = list(resolution = 0.30, dims = 1:30)
Epithelial_cell <- FindNeighbors(Epithelial_cell, reduction = "integrated.cca", dims = 1:30)
Epithelial_cell <- FindClusters(Epithelial_cell, resolution = 0.30, cluster.name = "cca_clusters")
Epithelial_cell <- RunUMAP(Epithelial_cell, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
DimPlot(Epithelial_cell, group.by = "cca_clusters", label = TRUE)
ggsave("C:/Users/wangyanhong/Desktop/LNLC/figure/umap_Epithelial_cell_Tumor_after_integrated.cca.pdf", width = 10, height = 8)
print(table(Epithelial_cell@meta.data$seurat_clusters))
print(table(Epithelial_cell@meta.data$orig.ident))
# 定义 Epithelial 细胞亚群
complete_celltype_list <- data.frame(
  ClusterID = 0:8,  # 假设 ClusterID 从 0 到 4
  celltype = c("PTN_malignant_Epithelial_Cells_C0", "CTSE_malignant_Epithelial_Cells_C1","HHIP_Epithelial_Cells_C2", "MRPS33_malignant_Epithelial_Cells_C3","ARMH2_malignant_Epithelial_Cells_C4", 
               "NABP1_malignant_Epithelial_Cells_C5", "EFCAB10_malignant_Epithelial_Cells_C6","CPA3_Epithelial_Cells_C7",  "RTKN2_malignant_Epithelial_Cells_C8")
)
# 更新 celltype 列为 "NA"
Epithelial_cell@meta.data$celltype <- "NA"
# 确保 seurat_clusters 是字符型
Epithelial_cell@meta.data$seurat_clusters <- as.character(Epithelial_cell@meta.data$seurat_clusters)
# 更新 celltype 列
for (i in 1:nrow(complete_celltype_list)) {
  cluster_id <- as.character(complete_celltype_list$ClusterID[i])  # 确保 cluster_id 为字符型
  cell_type <- complete_celltype_list$celltype[i]
  
  matching_cells <- which(Epithelial_cell@meta.data$seurat_clusters == cluster_id)
  if (length(matching_cells) > 0) {
    Epithelial_cell@meta.data[matching_cells, 'celltype'] <- cell_type
  }
}
# 提取 Epithelial 细胞的注释信息
Epithelial_cell_metadata <- Epithelial_cell@meta.data[, c("celltype", "seurat_clusters")]
# 确保 scRNAseq.sct.int 和 Epithelial_cell 的细胞名匹配
common_cells <- intersect(rownames(scRNAseq.sct.int@meta.data), rownames(Epithelial_cell_metadata))
length(common_cells)  # 检查共有细胞数
# 将注释信息添加到原始数据对象中
scRNAseq.sct.int@meta.data[common_cells, "sub_celltype"] <- Epithelial_cell_metadata[common_cells, "celltype"]
# 检查结果
print(table(scRNAseq.sct.int@meta.data$sub_celltype))
##### NK 细胞亚群分析NK_cell = list(resolution = 0.10, dims = 1:40)####
NK_cell <- FindNeighbors(NK_cell, reduction = "integrated.cca", dims = 1:40)
NK_cell <- FindClusters(NK_cell, resolution = 0.10, cluster.name = "cca_clusters")
NK_cell <- RunUMAP(NK_cell, reduction = "integrated.cca", dims = 1:40, reduction.name = "umap.cca")
DimPlot(NK_cell, group.by = "cca_clusters", label = TRUE)
ggsave("C:/Users/wangyanhong/Desktop/LNLC/figure/umap_NK_cell_Tumor_after_integrated.cca.pdf", width = 10, height = 8)
print(table(NK_cell@meta.data$seurat_clusters))
print(table(NK_cell@meta.data$orig.ident))
# 定义 NK 细胞亚群
complete_celltype_list <- data.frame(
  ClusterID = 0:2,  # 假设 ClusterID 从 0 到 2
  celltype = c("MYOM2_mature_NK_Cells_C0", "SPRY1_immature_NK_Cells_C1", "CD5_NK_T_Cells_C2" ))
# 更新 celltype 列为 "NA"
NK_cell@meta.data$celltype <- "NA"
# 确保 seurat_clusters 是字符型
NK_cell@meta.data$seurat_clusters <- as.character(NK_cell@meta.data$seurat_clusters)
# 更新 celltype 列
for (i in 1:nrow(complete_celltype_list)) {
  cluster_id <- as.character(complete_celltype_list$ClusterID[i])  # 确保 cluster_id 为字符型
  cell_type <- complete_celltype_list$celltype[i]
  
  matching_cells <- which(NK_cell@meta.data$seurat_clusters == cluster_id)
  if (length(matching_cells) > 0) {
    NK_cell@meta.data[matching_cells, 'celltype'] <- cell_type
  }
}
# NK 细胞亚群的 marker gene 列表
marker_genes <- c("FCGR3A","PRF1","KLRC1","CD3C","CD3G")
# 绘制 DotPlot 图
dot_plot <- DotPlot(NK_cell, features = unique(marker_genes), group.by = "celltype") + RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))  # 调整 x 轴标签的显示
# 保存 DotPlot 图
ggsave("NK_cell/DotPlot_NK_cell.pdf", plot = dot_plot, width = 10, height = 8)
# 绘制带有注释的 UMAP 图
umap_plot <- DimPlot(NK_cell, reduction = "umap", group.by = "celltype", label = TRUE, pt.size = 0.5) +
  ggtitle("UMAP Plot of NK Cells") +
  theme(plot.title = element_text(hjust = 0.5))  # 标题居中
# 保存 UMAP 图
ggsave("NK_cell/UMAP_NK_cell_annotated.pdf", plot = umap_plot, width = 10, height = 8)

# 提取 NK 细胞的注释信息
NK_cell_metadata <- NK_cell@meta.data[, c("celltype", "seurat_clusters")]
# 确保 scRNAseq.sct.int 和 NK_cell 的细胞名匹配
common_cells <- intersect(rownames(scRNAseq.sct.int@meta.data), rownames(NK_cell_metadata))
length(common_cells)  # 检查共有细胞数
# 将注释信息添加到原始数据对象中
scRNAseq.sct.int@meta.data[common_cells, "sub_celltype"] <- NK_cell_metadata[common_cells, "celltype"]
# 检查结果
print(table(scRNAseq.sct.int@meta.data$sub_celltype))
# 绘制 亚群细分的UMAP 图，显示整合后的细胞亚群注释
umap_plot <- DimPlot(scRNAseq.sct.int, reduction = "umap", group.by = "sub_celltype",  repel = TRUE, pt.size = 0.2) +
  theme(legend.position = "right", legend.text = element_text(size = 8), legend.key.size = unit(0.3, "cm")) +
  labs(title = "UMAP with Sub-celltype Annotations") +
  guides(colour = guide_legend(ncol = 1))  # 调整图注列数
# 显示并保存 UMAP 图
print(umap_plot)
ggsave("C:/Users/wangyanhong/Desktop/LNLC/figure/UMAP_with_Sub_celltype_Annotations01.pdf", plot = umap_plot, width = 15, height = 10)
saveRDS(scRNAseq.sct.int, "scRNAseq.sct.int_cell_FAM.rds")
scRNAseq.sct.int <- readRDS("scRNAseq.sct.int_cell_FAM.rds")

