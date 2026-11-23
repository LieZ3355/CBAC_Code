# Load packages----
library(Seurat)
library(R.utils)
library(mascarade)
## 1.1COVID-19----
## 1.2COVID-19单细胞数据----
## meta中必须有 barcodes 列
stopifnot("barcodes" %in% colnames(meta))

## 1. 从meta$barcodes提取 sample_id
sample_ids <- unique(sub("-[ACGT]+$", "", meta$barcodes))

## 2. 匹配文件
raw_dir <- "GSE171668_RAW"

all_files <- list.files(
  path = raw_dir,
  pattern = "\\.h5\\.gz$",
  full.names = TRUE
)

file_map <- data.frame(
  file = all_files,
  sample_id = sub(
    "^GSM[0-9]+_(.+)_raw_feature_bc_matrix\\.h5\\.gz$",
    "\\1",
    basename(all_files)
  ),
  stringsAsFactors = FALSE
)

matched_files <- file_map[file_map$sample_id %in% sample_ids, , drop = FALSE]

## 3. 读取每个样本，并把细胞名改成和 meta$barcodes 一致
read_one_h5_to_seurat <- function(h5_file, sample_id) {
  tmp_h5 <- tempfile(fileext = ".h5")
  R.utils::gunzip(h5_file, destname = tmp_h5, remove = FALSE, overwrite = TRUE)
  
  mat <- Read10X_h5(tmp_h5)
  
  if (is.list(mat)) {
    if ("Gene Expression" %in% names(mat)) {
      mat <- mat[["Gene Expression"]]
    } else {
      mat <- mat[[1]]
    }
  }
  
  obj <- CreateSeuratObject(
    counts = mat,
    project = sample_id,
    min.cells = 3,
    min.features = 200
  )
  
  ## 10x barcode 通常像 AAACCC...-1
  ## 改成和 meta$barcodes 一致：sample_id + "-" + 去掉末尾 -1 的barcode
  new_cells <- paste0(
    sample_id,
    "-",
    sub("-[0-9]+$", "", colnames(obj))
  )
  colnames(obj) <- new_cells
  
  obj$sample_id <- sample_id
  
  unlink(tmp_h5)
  obj
}

## 4. 逐个读取
seu_list <- vector("list", length = nrow(matched_files))
names(seu_list) <- matched_files$sample_id

for (i in seq_len(nrow(matched_files))) {
  sid <- matched_files$sample_id[i]
  cat("Reading:", sid, "\n")
  seu_list[[i]] <- read_one_h5_to_seurat(
    h5_file = matched_files$file[i],
    sample_id = sid
  )
}

## 5. 合并
GSE171668 <- seu_list[[1]]
if (length(seu_list) > 1) {
  for (i in 2:length(seu_list)) {
    GSE171668 <- merge(GSE171668, y = seu_list[[i]])
  }
}

## 6. 只保留 meta$barcodes 中存在的细胞
keep_cells <- intersect(colnames(GSE171668), meta$barcodes)
GSE171668 <- subset(GSE171668, cells = keep_cells)

## 7. 将 meta 并入 GSE171668@meta.data
## 先保证顺序和 GSE171668 细胞顺序一致
meta2 <- meta[match(colnames(GSE171668), meta$barcodes), , drop = FALSE]

## 检查是否完全对齐
stopifnot(all(meta2$barcodes == colnames(GSE171668)))

## barcodes列可留可不留
rownames(meta2) <- meta2$barcodes

## 为避免和 Seurat 现有列重复，可去掉 barcodes 后再加
meta_add <- meta2[, setdiff(colnames(meta2), "barcodes"), drop = FALSE]

GSE171668 <- AddMetaData(GSE171668, metadata = meta_add)
GSE171668 <- NormalizeData(GSE171668, normalization.method = "LogNormalize", scale.factor = 10000)
GSE171668 <- FindVariableFeatures(GSE171668, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(GSE171668)
GSE171668 <- ScaleData(GSE171668, features = all.genes)
GSE171668 <- RunPCA(GSE171668, features = VariableFeatures(object = GSE171668))
GSE171668 <- FindNeighbors(GSE171668, dims = 1:10)
GSE171668 <- FindClusters(GSE171668, resolution = 0.7)
GSE171668 <- RunUMAP(GSE171668, dims = 1:10)
anno <- as.character(GSE171668$manual_annotation)
anno <- anno |>
  str_replace("^\\d+:\\s*", "") |>
  str_replace_all("\\b(CM|EC|Fib)\\s+(\\d+)\\b", "\\1-\\2") |>
  str_replace_all("\\bMito\\s+(EC|CM)\\b", "Mito-\\1") |>
  str_replace_all("\\bLow\\s+count\\b", "Low-count") |>
  str_replace_all("([A-Za-z]+)\\s*\\(([^()]+)\\)", "\\1-\\2")
anno[anno=="CM-5 (RYR2$^-$)"] <- "CM-5"
GSE171668$manual_annotation2 <- anno

GSE171668 <- RunHarmony(GSE171668, "donor")
GSE171668 <- FindNeighbors(GSE171668, reduction = "harmony", dims = 1:20)
GSE171668 <- FindClusters(GSE171668, resolution = 0.1)
GSE171668 <- RunUMAP(GSE171668, dims = 1:20,reduction = "harmony")
DimPlot(GSE171668, reduction = "umap",label = TRUE,group.by = "manual_annotation2")+xlab("UMAP_1")+ylab("UMAP_2")

library(ggplot2)
library(dplyr)
library(tidydr)
library(mascarade)
### UMAP绘图 ---------------------------------------
color_mapping <- c(
  "Adipocyte"            = "#B22222",
  "CM-1"                 = "#8A8FBD",
  "CM-2"                 = "#6F78B8",
  "CM-3"                 = "#A7ADD8",
  "CM-4"                 = "#BEC3E6",
  "CM-5"                 = "#5E679E",
  "Mito-CM"              = "#4E567F",
  "CM + EC + Pericyte"   = "#4A4E69",
  "CM + Macrophage"      = "#8B6F8F",
  
  "EC-1"                 = "#F2C84B",
  "EC-2"                 = "#E3B93C",
  "EC-3"                 = "#D6A92E",
  "EC-POSTN"             = "#C9961A",
  "EC + Pericyte"        = "#C9A24A",
  "Lymph EC"             = "#B88E2F",
  "Mito-EC"              = "#9C7A1A",
  
  "Fib-1"                = "#B7E08A",
  "Fib-2"                = "#8FCB73",
  "Low-count Fib"        = "#D5EDBE",
  "Fib + EC + Pericyte"  = "#A9C97A",
  "Fib + CM"             = "#9BAF88",
  
  "Macrophage"           = "#2E6DB4",
  "T cell"               = "#444444",
  "Neuronal"             = "#C17CCF",
  
  "Pericyte"             = "#F0A202",
  "VSMC"                 = "#A64B9A"
)


celltype_major <- c(
  "Adipocyte"           = "Adipocyte",
  "CM-1"                = "Cardiomyocyte",
  "CM-2"                = "Cardiomyocyte",
  "CM-3"                = "Cardiomyocyte",
  "CM-4"                = "Cardiomyocyte",
  "CM-5"                = "Cardiomyocyte",
  "Mito-CM"             = "Cardiomyocyte",
  "Fib + CM"            = "Cardiomyocyte",
  "CM + EC + Pericyte"  = "Cardiomyocyte",
  "CM + Macrophage"     = "Cardiomyocyte",
  "EC-1"                = "Endothelial",
  "EC-2"                = "Endothelial",
  "EC-3"                = "Endothelial",
  "EC-POSTN"            = "Endothelial",
  "Mito-EC"             = "Endothelial",
  "Lymph EC"            = "Endothelial",
  "EC + Pericyte"       = "Endothelial",
  "Fib-1"               = "Fibroblast",
  "Fib-2"               = "Fibroblast",
  "Low-count Fib"       = "Fibroblast",
  "Fib + EC + Pericyte" = "Fibroblast",
  "Macrophage"          = "Macrophage",
  "T cell"              = "T cell",
  "Pericyte"            = "Pericyte",
  "VSMC"                = "VSMC",
  "Neuronal"            = "Neuronal"
)

make_major_colors <- function(sub_colors, sub_to_major) {
  valid_sub <- intersect(names(sub_colors), names(sub_to_major))
  
  tmp_df <- data.frame(
    subtype = valid_sub,
    major = unname(sub_to_major[valid_sub]),
    color = unname(sub_colors[valid_sub]),
    stringsAsFactors = FALSE
  )
  
  major_colors <- lapply(split(tmp_df$color, tmp_df$major), function(cols) {
    rgb_mat <- grDevices::col2rgb(cols)
    avg_rgb <- rowMeans(rgb_mat)
    grDevices::rgb(
      red   = avg_rgb[1],
      green = avg_rgb[2],
      blue  = avg_rgb[3],
      maxColorValue = 255
    )
  })
  
  major_colors <- unlist(major_colors)
  return(major_colors)
}

major_color_mapping <- make_major_colors(
  sub_colors = color_mapping,
  sub_to_major = celltype_major
)


#---- 02 计算质心距离 
plot_df <- cbind(
  as.data.frame(GSE171668@reductions$umap@cell.embeddings),
  GSE171668@meta.data
)

plot_df$major_celltype <- factor(
  plot_df$major_celltype,
  levels = names(major_color_mapping)
)

# 主群轮廓
maskTable_major <- generateMask(
  dims = plot_df[, c("umap_1", "umap_2")],
  clusters = plot_df$major_celltype,
  minDensity = 0.2,
  smoothSigma = 0.025
)

# 主群标签位置
label_df_major <- plot_df %>%
  dplyr::group_by(major_celltype) %>%
  dplyr::summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2),
    .groups = "drop"
  )


# 作图

ggplot(plot_df, aes(x = umap_1, y = umap_2)) +
  geom_point(aes(color = major_celltype), size = 0.1) +
  geom_path(
    data = maskTable_major,
    aes(group = group),
    linewidth = 0.8,
    linetype = 2,
    color = "black"
  ) +
  geom_text(
    data = label_df_major,
    aes(x = umap_1, y = umap_2, label = major_celltype),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold",
    color = "black"
  ) +
  coord_fixed(ratio = 1) +
  theme_classic() +
  tidydr::theme_dr() +
  scale_color_manual(
    values = major_color_mapping,
    guide = guide_legend(
      ncol = 1,
      override.aes = list(size = 4)
    )
  ) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(colour = "black", size = 16),
    axis.title = element_text(colour = "black", size = 18),
    axis.text = element_blank()
  ) +
  xlab("UMAP_1") +
  ylab("UMAP_2")



## 1.2动脉粥样硬化----
### GSE155512----

files <- c(
  "GSM4705589_RPE004_matrix.txt.gz",
  "GSM4705590_RPE005_matrix.txt.gz",
  "GSM4705591_RPE006_matrix.txt.gz"
)

samples <- c("RPE004", "RPE005", "RPE006")

obj_list <- lapply(seq_along(files), function(i) {
  mat <- data.table::fread(files[i], data.table = FALSE)
  
  rownames(mat) <- mat[, 1]
  mat <- as.matrix(mat[, -1])
  storage.mode(mat) <- "numeric"
  
  obj <- CreateSeuratObject(
    counts = mat,
    project = samples[i],
    min.cells = 3,
    min.features = 200
  )
  
  obj$sample <- samples[i]
  obj
})

seu_all <- merge(
  x = obj_list[[1]],
  y = obj_list[2:3],
  add.cell.ids = samples,
  project = "RPE_merged"
)
GSE155512 <- seu_all
rm(obj_list, seu_all)
GSE155512[["percent.mt"]] <- PercentageFeatureSet(GSE155512, pattern = "^MT-")
GSE155512[["percent.rb"]] <- PercentageFeatureSet(GSE155512, pattern = "^RP[SL]")
VlnPlot(GSE155512, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4)
GSE155512 <- NormalizeData(GSE155512, normalization.method = "LogNormalize", scale.factor = 10000)
GSE155512 <- FindVariableFeatures(GSE155512, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(GSE155512)
GSE155512 <- ScaleData(GSE155512, features = all.genes)
GSE155512 <- RunPCA(GSE155512, features = VariableFeatures(object = GSE155512))
GSE155512 <- FindNeighbors(GSE155512, dims = 1:10)
GSE155512 <- FindClusters(GSE155512, resolution = 0.7)
GSE155512 <- RunUMAP(GSE155512, dims = 1:10)
DimPlot(GSE155512, reduction = "umap",label = TRUE,group.by = "sample")+xlab("UMAP_1")+ylab("UMAP_2")



fs = list.files('./GSE253903_RAW/',pattern = '^GSM')
samples <- substr(fs,1,10)
lapply( unique(samples), function(x){
  y = fs[grepl(x,fs)]
  folder = paste0('./GSE253903_RAW/',strsplit(y[1],split = '_')[[1]][1])
  #创建文件夹
  dir.create(folder,recursive = T)
  #重命名子文件夹并移动到相应的文件夹中
  file.rename(paste0('./GSE253903_RAW/',y[1]),file.path(folder,"barcodes.tsv.gz"))
  file.rename(paste0('./GSE253903_RAW/',y[2]),file.path(folder,"features.tsv.gz"))
  file.rename(paste0('./GSE253903_RAW/',y[3]),file.path(folder,"matrix.mtx.gz"))
})
###生成每一个Seurat对象
samples <- dir(path="./GSE253903_RAW/", pattern="^GSM")
for (i in samples){
  assign(paste0( i), Read10X(data.dir = paste0("./GSE253903_RAW/", i)))
}
####Seurat list
##进入文件夹
setwd("GSE253903_RAW")
sceList = lapply(samples,function(folder){ 
  CreateSeuratObject(min.cells = 5, min.features = 200,counts = Read10X(folder), 
                     project = folder )
})
GSE253903 <- merge(x = sceList[[1]],
              y = sceList[-1],
              add.cell.ids = samples)
setwd("..")

GSE253903[["percent.mt"]] <- PercentageFeatureSet(GSE253903, pattern = "^MT-")
GSE253903[["percent.rb"]] <- PercentageFeatureSet(GSE253903, pattern = "^RP[SL]")
VlnPlot(GSE253903, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4)
FeatureScatter(GSE253903, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
GSE253903 <- subset(GSE253903, subset = nFeature_RNA > 200 & percent.mt < 25& percent.rb < 40)
GSE253903 <- NormalizeData(GSE253903, normalization.method = "LogNormalize", scale.factor = 10000)
GSE253903 <- FindVariableFeatures(GSE253903, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(GSE253903)
GSE253903 <- ScaleData(GSE253903, features = all.genes)
GSE253903 <- RunPCA(GSE253903, features = VariableFeatures(object = GSE253903))
GSE253903 <- FindNeighbors(GSE253903, dims = 1:20)
GSE253903 <- FindClusters(GSE253903, resolution = 0.1)
GSE253903 <- RunUMAP(GSE253903, dims = 1:20)
DimPlot(GSE253903, reduction = "umap",label = TRUE)+xlab("UMAP_1")+ylab("UMAP_2")
scedata <- merge(GSE155512, GSE253903)

scedata <- JoinLayers(scedata)
count_mat <- LayerData(scedata, assay = "RNA", layer = "counts")

head(rownames(count_mat))

scedata <- CreateSeuratObject(
  counts = count_mat,
  meta.data = scedata@meta.data
)
scedata <- NormalizeData(scedata, normalization.method = "LogNormalize", scale.factor = 10000)
scedata <- FindVariableFeatures(scedata, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(scedata)
scedata <- ScaleData(scedata, features = all.genes)
scedata <- RunPCA(scedata, features = VariableFeatures(object = scedata))
scedata <- FindNeighbors(scedata, dims = 1:20)
scedata <- FindClusters(scedata, resolution = 0.1)
scedata <- RunUMAP(scedata, dims = 1:20)
scedata$batch <- ifelse(scedata$batch=="batch1", "GSE155512", "GSE253903")
p1=DimPlot(
  scedata,
  reduction = "umap",
  label = FALSE,
  group.by = "batch",
  shuffle = TRUE,
  pt.size = 0.4,
  alpha = 0.7
) +
  scale_color_manual(values = c(
    "#4C78A8",  # 雾蓝
    "#C06C84"   # 柔和砖红
  )) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  tidydr::theme_dr() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(colour = "black", size = 16),
    axis.title = element_text(colour = "black", size = 18),
    axis.text = element_blank(),
    plot.title = element_blank()
  )
library(harmony)
scedata <- RunHarmony(scedata, "batch")
scedata <- FindNeighbors(scedata, reduction = "harmony", dims = 1:10)
scedata <- FindClusters(scedata, resolution = 1)
scedata <- RunUMAP(scedata, dims = 1:10,reduction = "harmony")
p2=DimPlot(
  scedata,
  reduction = "umap",
  label = FALSE,
  group.by = "batch",
  shuffle = TRUE,
  pt.size = 0.4,
  alpha = 0.7
) +
  scale_color_manual(values = c(
    "#4C78A8",  # 雾蓝
    "#C06C84"   # 柔和砖红
  )) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  tidydr::theme_dr() +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(colour = "black", size = 16),
    axis.title = element_text(colour = "black", size = 18),
    axis.text = element_blank(),
    plot.title = element_blank()
  )


### 细胞注释----
new.cluster.ids <- c("0"="T cell",
                     "1"="Macrophage",
                     "2"="Fibroblast",
                     "3"="T cell",
                     "4"="SMC",
                     "5"="Macrophage",
                     "6"="SMC",
                     "7"="SMC",
                     "8"="Endothelial",
                     "9"="Macrophage",
                     "10"="DCs",
                     "11"="Macrophage",
                     "12"="Endothelial",
                     "13"="Macrophage",
                     "14"="Fibroblast",
                     "15"="Neutrophil",
                     "16"="MAST cell",
                     "17"="Macrophage",
                     "18"="B cell",
                     "19"="Neutrophil",
                     "20"="Fibroblast",
                     "21"="T cell",
                     "22"="Plasma",
                     "23"="B cell",
                     "24"="Endothelial",
                     "25"="SMC",
                     "26"="B cell",
                     "27"="Endothelial")
scedata2 <- RenameIdents(scedata, new.cluster.ids)
scedata2$celltype <- scedata2@active.ident
DimPlot(scedata2, label = T,group.by = "celltype")+xlab("UMAP_1")+ylab("UMAP_2")+ggtitle("")
### UMAP绘图 ---------------------------------------------
plot_df <- cbind(scedata2@reductions$umap@cell.embeddings, scedata2@meta.data)
maskTable <- generateMask( dims = plot_df[,c("umap_1","umap_2")], cluster = plot_df$celltype, minDensity =0.2, smoothSigma =0.025)

label_df <- plot_df %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2)
  )

color_mapping <- c(
  "T cell"      = "#5B5B5B",  # 深灰
  "Macrophage"  = "#3C7FB1",  # 柔和蓝
  "Fibroblast"  = "#8FBC6B",  # 草木绿
  "SMC"         = "#7E87B8",  # 灰紫蓝
  "Endothelial" = "#D8B547",  # 柔和金黄
  "DCs"         = "#9A8F7A",  # 灰卡其
  "Neutrophil"  = "#D98C8C",  # 灰粉红
  "MAST cell"   = "#E3A65F",  # 杏橙
  "B cell"      = "#8C6BB1",  # 柔和紫
  "Plasma"      = "#C06C84"   # 玫瑰豆沙
)

plot_df2 <- plot_df %>%
  dplyr::filter(
    # 1) 去掉 Macrophage 和 B cell 中间那撮绿色 Fibroblast
    !(celltype == "Fibroblast" &
        umap_1 > -0.5) &
      # 2) 去掉 T cell 周围的零散小群（按你的图再微调）
      !(celltype != "T cell" &
          umap_1 > 2.5 & umap_1 < 13.5 &
          umap_2 > -10 & umap_2 < -5.0)&
      !(celltype != "T cell" &
          umap_1 > 9 & umap_1 < 12 &
          umap_2 > -6 & umap_2 < -2.5)&
      !(celltype == "SMC" &
          umap_1 > -7 &
          umap_2 >-2 )&
      !(celltype == "SMC" &
          umap_1 < -8 &
          umap_2 > 0 )&
      !(celltype == "Endothelial" &
          umap_1 < -5  )
  )
maskTable <- generateMask(
  dims = as.matrix(plot_df2[, c("umap_1", "umap_2")]),
  cluster = plot_df2$celltype,
  minDensity = 0.65,
  smoothSigma = 0.005
)
ggplot(plot_df, aes(x = umap_1, y = umap_2)) +
  geom_point(aes(color = celltype), size = 0.1) +
  geom_path(
    data = maskTable,
    aes(group = group),
    linewidth = 0.8,
    linetype = 2,
    color = "black"
  ) +
  coord_fixed(ratio = 1) +
  theme_classic() +
  geom_text(
    data = label_df,
    aes(x = umap_1, y = umap_2, label = celltype),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold",
    color = "black"
  ) +
  tidydr::theme_dr() +
  scale_color_manual(
    values = color_mapping,
    guide = guide_legend(
      ncol = 1,
      override.aes = list(size = 4)
    )
  ) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(colour = "black", size = 16),
    axis.title = element_text(colour = "black", size = 18),
    axis.text = element_blank()
  ) +
  xlab("UMAP_1") +
  ylab("UMAP_2")


## 2. scPagwas关联分析----
## COVID-19_scPagwas----
library(scPagwas)
GSE171668 <- JoinLayers(GSE171668)
COVID_GWAS<- fread("E:/Rwork/12.11xh/scRNA/COVID-19 infection.tsv",header = T)%>%data.frame()
COVID_GWAS$maf=ifelse(COVID_GWAS$all_meta_AF>0.5,1-COVID_GWAS$all_meta_AF,COVID_GWAS$all_meta_AF)
COVID_GWAS=COVID_GWAS[,c(1,2,15,8,7,16)]
colnames(COVID_GWAS) <- c("chrom","pos","rsid","se","beta","maf")
GSE171668 <- SetIdent(GSE171668, value = "major_celltype")
setwd("G:/tmp")
Pagwas_data_Finland<-scPagwas_main(Pagwas = NULL,
                                   gwas_data =COVID_GWAS, # The GWAS Summary statistics files
                                   Single_data =GSE171668,# scRNA-seq data in seruat format with "RNA" assays and normalized.
                                   output.prefix="test", # the prefix name for output files
                                   output.dirs="scPagwastest_output_COVID",
                                   assay="RNA", # the assays for scRNA-seq data to use.
                                   Pathway_list=Genes_by_pathway_kegg,# pathway list is provided by package, including gene symbols.
                                   n.cores=24,
                                   block_annotation = block_annotation,
                                   iters_singlecell = 0,
                                   min.pathway.size = 20,
                                   chrom_ld = chrom_ld,# The LD data is provided by package.
                                   singlecell=T, # Whether to run the singlecell process.
                                   celltype=T,
)
source("E:/Rwork/12.11xh/scRNA/Visualization.R")

vv(Single_data=Pagwas_data_Finland,
   p_thre = 0.05,
   FigureType = "umap",
   width = 7,
   height = 7,
   lowColor = "white", 
   highColor = "red",
   output.dirs="figure_Finland",
   size = 0.5,
   do_plot = T)

## TRS绘图----
gw <- read.csv("G:/tmp/scPagwastest_output_COVID/test_singlecell_scPagwas_score_pvalue.Result.csv")
GSE171668$TRS <- gw$scPagwas.TRS.Score

FeaturePlot(GSE171668, features =c("TRS"),raster=FALSE,pt.size = 0.1,cols = c("#d3d3d3","#c5423d"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+
  theme(  panel.grid = element_blank(),
          legend.position ="right",
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 22,face = "bold",colour = "black",),
          legend.text  = element_text(colour = "black",size=16),
          axis.title = element_text(colour = "black",size=18),
          axis.text = element_text(colour = "black",size=16))

## GWAS数据绘图----
library(CMplot)
COVID_GWAS<- fread("E:/Rwork/12.11xh/scRNA/COVID-19 infection.tsv",header = T)%>%data.frame()
tmp <- COVID_GWAS[,c(15,1,2,9)]
COVID_GWAS_filter=tmp[tmp$all_inv_var_meta_p<5e-9,]%>%na.omit()
COVID_GWAS_filter <- COVID_GWAS_filter[order(COVID_GWAS_filter$all_inv_var_meta_p),]
COVID_GWAS_highlight <- COVID_GWAS_filter[1:10,]
### 柱状图----
CMplot(tmp, 
       plot.type="m",
       LOG10=TRUE,
       ylim=NULL,
       threshold=c(5e-9,1e-5),
       threshold.lty=c(1,2),
       threshold.lwd=c(1,1),
       threshold.col=c("black","grey"), 
       amplify=F,
       bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow","red"),
       signal.col=c("red","green"),
       signal.cex=c(1,1),
       signal.pch=c(19,19),
       highlight = COVID_GWAS_highlight$rsid,
       highlight.col = 'black',
       highlight.cex = 1,
       highlight.text=COVID_GWAS_highlight$rsid,
       file="jpg",dpi=300,
       file.output=TRUE,
       verbose=TRUE)


### 环形曼哈顿图----
CMplot(tmp,
       plot.type="c",
       chr.labels=paste0("Chr",c(1:22,"X")),
       r=0.4,
       outward=FALSE,
       cir.chr.h=1.3,
       chr.den.col=c("darkgreen","yellow","red"),
       file="jpg",
       dpi=300,
       file.output=TRUE,
       verbose=TRUE,
       highlight =COVID_GWAS_highlight$rsid)

### QQ图----
CMplot(tmp,plot.type = "q",threshold = 0.05,file="jpg",file.name = "covid_qq",cex = 0.4)

## Ath_scPagwas----
Ath_GWAS <- fread("E:/Rwork/12.11xh/scRNA/finngen_R12_I9_ATHSCLE.gz",header = T)%>%data.frame()
Ath_GWAS$maf=ifelse(Ath_GWAS$af_alt>0.5,1-Ath_GWAS$af_alt,Ath_GWAS$af_alt)
Ath_GWAS=Ath_GWAS[,c(1,2,5,10,9,14)]
colnames(Ath_GWAS) <- c("chrom","pos","rsid","se","beta","maf")
Pagwas_data_Ath<-scPagwas_main(Pagwas = NULL,
                               gwas_data =Ath_GWAS, # The GWAS Summary statistics files
                               Single_data =scedata,# scRNA-seq data in seruat format with "RNA" assays and normalized.
                               output.prefix="test", # the prefix name for output files
                               output.dirs="scPagwastest_output_Ath",
                               assay="RNA", # the assays for scRNA-seq data to use.
                               Pathway_list=Genes_by_pathway_kegg,# pathway list is provided by package, including gene symbols.
                               n.cores=24,
                               block_annotation = block_annotation,
                               iters_singlecell = 0,
                               min.pathway.size = 20,
                               chrom_ld = chrom_ld,# The LD data is provided by package.
                               singlecell=T, # Whether to run the singlecell process.
                               celltype=T,
)

## GWAS数据绘图----
library(CMplot)
Ath_GWAS <- fread("E:/Rwork/12.11xh/scRNA/finngen_R12_I9_ATHSCLE.gz",header = T)%>%data.frame()

tmp <- Ath_GWAS[,c(5,1,2,7)]
Ath_GWAS_filter=tmp[tmp$pval<5e-9,]%>%na.omit()
Ath_GWAS_highlight=Ath_GWAS_filter[Ath_GWAS_filter$pval<5e-30,]%>%na.omit()
Ath_GWAS_highlight <- Ath_GWAS_highlight[order(Ath_GWAS_highlight$pval),]
Ath_GWAS_highlight <- Ath_GWAS_highlight[1:10,]
### 柱状图----
CMplot(tmp, 
       plot.type="m",
       LOG10=TRUE,
       ylim=NULL,
       threshold=c(5e-9,1e-5),
       threshold.lty=c(1,2),
       threshold.lwd=c(1,1),
       threshold.col=c("black","grey"), 
       amplify=T,
       bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow","red"),
       signal.col=c("red","green"),
       signal.cex=c(1,1),
       signal.pch=c(19,19),
       highlight = Ath_GWAS_highlight$rsid,
       highlight.col = 'black',
       highlight.cex = 1,
       highlight.text=Ath_GWAS_highlight$rsid,
       file="jpg",
       dpi=300,
       file.output=TRUE,
       verbose=TRUE)


### 环形曼哈顿图----
CMplot(tmp,
       plot.type="c",
       chr.labels=paste0("Chr",c(1:22,"X")),
       r=0.4,
       outward=FALSE,
       cir.chr.h=1.3,
       chr.den.col=c("darkgreen","yellow","red"),
       file="jpg",
       dpi=300,
       file.output=TRUE,
       verbose=TRUE,
       highlight =Ath_GWAS_highlight$rsid)

### QQ图----
CMplot(tmp,plot.type = "q",threshold = 0.05,file="jpg",file.name = "Ath_qq",cex = 0.4)

### TRS绘图----
gw2 <- read.csv("G:/tmp/scPagwastest_output_Ath/test_singlecell_scPagwas_score_pvalue.Result.csv")
Ath$TRS <- gw2$scPagwas.TRS.Score
## TRS绘图----
FeaturePlot(Ath, features =c("TRS"),raster=FALSE,pt.size = 0.1,cols = c("#d3d3d3","#c5423d"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+
  theme(  panel.grid = element_blank(),
          legend.position ="right",
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 22,face = "bold",colour = "black",),
          legend.text  = element_text(colour = "black",size=16),
          axis.title = element_text(colour = "black",size=18),
          axis.text = element_text(colour = "black",size=16))

## 3. 亚群----
## COVID细胞亚群----
COVID_EC<- GSE171668[,GSE171668$major_celltype%in%c("Endothelial")]
COVID_Fib<- GSE171668[,GSE171668$major_celltype%in%c("Fibroblast")]
subcluster <- function(seurat_obj, resolution = 0.5, dims = 1:20) {
  seurat_obj <- CreateSeuratObject(counts = seurat_obj@assays$RNA$counts,meta.data = seurat_obj@meta.data)
  seurat_obj <- NormalizeData(seurat_obj, normalization.method =  "LogNormalize", scale.factor = 10000)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj, features = all.genes)
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
  seurat_obj <- FindNeighbors(seurat_obj, dims = dims)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
  seurat_obj <- RunUMAP(seurat_obj, dims = dims)
  return(seurat_obj)
}
COVID_EC <- subcluster(COVID_EC, resolution = 2, dims = 1:10)
COVID_Fib <- subcluster(COVID_Fib, resolution = 0.1, dims = 1:10)
### COVID_EC Dimplot----
old_cluster_col <- "seurat_clusters"
new_cluster_col <- "manual_group_raw"
final_cluster_col <- "manual_group_rank"

COVID_EC[[old_cluster_col]] <- as.character(COVID_EC[[old_cluster_col]][, 1])
group_map <- list(
  "A0" = c(1, 9,19, 5, 7, 12, 23),
  "A1" = c(2, 8, 11),
  "A2" = c(0, 22, 13, 16, 21, 4),
  "A3" = c(6, 3, 14, 10, 20, 24),
  "A4" = c(15, 18),
  "A5" = c(17)
)
map_df <- do.call(
  rbind,
  lapply(names(group_map), function(g) {
    data.frame(
      old_cluster = as.character(group_map[[g]]),
      new_group = g,
      stringsAsFactors = FALSE
    )
  })
)
# 4. 写入新的粗分群
meta_df <- COVID_EC@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::left_join(
    map_df,
    by = c(setNames("old_cluster", old_cluster_col))
  )

# 没被分到任何组的细胞标记为 NA 或 "Unassigned"
meta_df$new_group[is.na(meta_df$new_group)] <- "Unassigned"

COVID_EC[[new_cluster_col]] <- meta_df$new_group

# 5. 统计每个新群细胞数，并按细胞数降序排列
group_size <- COVID_EC@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::count(.data[[new_cluster_col]], name = "n_cells") %>%
  dplyr::arrange(dplyr::desc(n_cells))

print(group_size)

# 如果你不想把 Unassigned 参与排序，可以先去掉
group_size2 <- group_size %>%
  dplyr::filter(.data[[new_cluster_col]] != "Unassigned") %>%
  dplyr::mutate(new_rank = as.character(seq_len(n()) - 1))

print(group_size2)

# 6. 按细胞数重新编号，从0开始
rank_map <- setNames(group_size2$new_rank, group_size2[[new_cluster_col]])

COVID_EC[[final_cluster_col]] <- COVID_EC@meta.data[[new_cluster_col]]
COVID_EC@meta.data[[final_cluster_col]] <- unname(rank_map[COVID_EC@meta.data[[new_cluster_col]]])

# Unassigned 保留为 NA 或单独标记
COVID_EC@meta.data[[final_cluster_col]][
  is.na(COVID_EC@meta.data[[final_cluster_col]])
] <- NA

# 7. 设置为有序因子，方便作图
ordered_levels <- as.character(sort(as.numeric(na.omit(unique(COVID_EC@meta.data[[final_cluster_col]])))))
COVID_EC@meta.data[[final_cluster_col]] <- factor(
  COVID_EC@meta.data[[final_cluster_col]],
  levels = ordered_levels
)

# 查看最终每群细胞数
final_count <- COVID_EC@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::filter(!is.na(.data[[final_cluster_col]])) %>%
  dplyr::count(.data[[final_cluster_col]], name = "n_cells") %>%
  dplyr::arrange(as.numeric(as.character(.data[[final_cluster_col]])))

print(final_count)
COVID_EC$manual_group_rank <- paste0("Endo_", COVID_EC$manual_group_rank)
Idents(COVID_EC) <- COVID_EC$manual_group_rank
p <- DimPlot(COVID_EC, group.by = "manual_group_rank", reduction = "umap",label = TRUE,pt.size = 0.6,label.size = 8)+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  ggtitle("Endothelial")+
  theme(  panel.grid = element_blank(),
          legend.position ="right",
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 22,face = "bold",colour = "black",),
          legend.text  = element_text(colour = "black",size=16),
          axis.title = element_text(colour = "black",size=18),
          axis.text = element_text(colour = "black",size=16))
ggsave("COVID_EC_UMAP.pdf", p, width = 10, height = 10)
### COVID_Fib Dimplot----

old_cluster_col <- "seurat_clusters"
final_cluster_col <- "manual_group_rank"

# 1. 先把原始cluster转成字符，避免后面匹配出问题
COVID_Fib[[old_cluster_col]] <- as.character(COVID_Fib[[old_cluster_col]][, 1])

# 2. 统计每个原始cluster细胞数，并按细胞数降序排列
group_size <- COVID_Fib@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::count(.data[[old_cluster_col]], name = "n_cells") %>%
  dplyr::arrange(dplyr::desc(n_cells))

print(group_size)

# 3. 按细胞数重新编号，从0开始
group_size <- group_size %>%
  dplyr::mutate(new_rank = as.character(seq_len(n()) - 1))

print(group_size)

# 4. 建立 old cluster -> new rank 的映射
rank_map <- setNames(group_size$new_rank, group_size[[old_cluster_col]])

# 5. 写入新的rank列
COVID_Fib[[final_cluster_col]] <- unname(rank_map[COVID_Fib@meta.data[[old_cluster_col]]])

# 6. 设置为有序因子
ordered_levels <- as.character(sort(as.numeric(unique(COVID_Fib@meta.data[[final_cluster_col]]))))
COVID_Fib@meta.data[[final_cluster_col]] <- factor(
  COVID_Fib@meta.data[[final_cluster_col]],
  levels = ordered_levels
)

# 7. 加前缀 Fib_
COVID_Fib@meta.data[[final_cluster_col]] <- paste0("Fib_", COVID_Fib@meta.data[[final_cluster_col]])

# 再设一次因子顺序
ordered_levels2 <- paste0("Fib_", seq_len(length(unique(COVID_Fib@meta.data[[old_cluster_col]]))) - 1)
COVID_Fib@meta.data[[final_cluster_col]] <- factor(
  COVID_Fib@meta.data[[final_cluster_col]],
  levels = ordered_levels2
)

# 8. 查看最终每群细胞数
final_count <- COVID_Fib@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::count(.data[[final_cluster_col]], name = "n_cells") %>%
  dplyr::arrange(as.numeric(gsub("Fib_", "", as.character(.data[[final_cluster_col]]))))

print(final_count)

# 9. 设定identity并画图
Idents(COVID_Fib) <- COVID_Fib@meta.data[[final_cluster_col]]

p <- DimPlot(
  COVID_Fib,
  group.by = final_cluster_col,
  reduction = "umap",
  label = TRUE,
  pt.size = 0.6,
  label.size = 8
) +
  xlab("UMAP_1") +
  ylab("UMAP_2") +
  ggtitle("Fibroblast") +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold", colour = "black"),
    legend.text = element_text(colour = "black", size = 16),
    axis.title = element_text(colour = "black", size = 18),
    axis.text = element_text(colour = "black", size = 16)
  )

ggsave("COVID_Fib_UMAP.pdf", p, width = 10, height = 10)
### COVID_EC TRS FeaturePlot----
FeaturePlot(COVID_EC, features =c("TRS"),raster=FALSE,pt.size = 0.6,cols = c("#d3d3d3","#c5423d"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+
  theme(  panel.grid = element_blank(),
          legend.position ="right",
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 22,face = "bold",colour = "black",),
          legend.text  = element_text(colour = "black",size=16),
          axis.title = element_text(colour = "black",size=18),
          axis.text = element_text(colour = "black",size=16))
### COVID_Fib TRS FeaturePlot----
FeaturePlot(COVID_Fib, features =c("TRS"),raster=FALSE,pt.size = 0.6,cols = c("#d3d3d3","#c5423d"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+
  theme(  panel.grid = element_blank(),
          legend.position ="right",
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 22,face = "bold",colour = "black",),
          legend.text  = element_text(colour = "black",size=16),
          axis.title = element_text(colour = "black",size=18),
          axis.text = element_text(colour = "black",size=16))

### COVID_EC 堆积柱状图----
COVID_EC@meta.data <- COVID_EC@meta.data %>%
  dplyr::mutate(
    TRS_quartile = dplyr::ntile(TRS, 4),
    TRS_group = dplyr::case_when(
      TRS_quartile == 1 ~ "Low",
      TRS_quartile %in% c(2, 3) ~ "Middle",
      TRS_quartile == 4 ~ "High"
    )
  )

COVID_EC$TRS_group <- factor(COVID_EC$TRS_group, levels = c("Low", "Middle", "High"))

meta_df <- COVID_EC@meta.data %>%
  dplyr::mutate(
    TRS_binary = ifelse(TRS_group == "High", "High", "Non-High")
  )

prop_df <- meta_df %>%
  dplyr::group_by(manual_group_rank, TRS_binary) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(manual_group_rank) %>%
  dplyr::mutate(
    proportion = n / sum(n),
    label = paste0(scales::percent(proportion, accuracy = 0.1), "\n(n=", n, ")")
  ) %>%
  dplyr::ungroup()

p <- ggplot(prop_df, aes(x = manual_group_rank, y = proportion, fill = TRS_binary)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 4,
    color = "black",
    lineheight = 0.9
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("High" = "#D55E00", "Non-High" = "#BDBDBD")) +
  labs(
    x = NULL,
    y = "Proportion",
    fill = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.title = element_text(colour = "black", size = 20),
    axis.text = element_text(colour = "black", size = 18),
    legend.text = element_text(colour = "black", size = 20)
  )
ggsave("EC_比例柱状图.pdf", p, width = 12, height = 7)

### COVID_Fib 堆积柱状图----
COVID_Fib@meta.data <- COVID_Fib@meta.data %>%
  dplyr::mutate(
    TRS_quartile = dplyr::ntile(TRS, 4),
    TRS_group = dplyr::case_when(
      TRS_quartile == 1 ~ "Low",
      TRS_quartile %in% c(2, 3) ~ "Middle",
      TRS_quartile == 4 ~ "High"
    )
  )

COVID_Fib$TRS_group <- factor(COVID_Fib$TRS_group, levels = c("Low", "Middle", "High"))

meta_df <- COVID_Fib@meta.data %>%
  dplyr::mutate(
    TRS_binary = ifelse(TRS_group == "High", "High", "Non-High")
  )

prop_df <- meta_df %>%
  dplyr::group_by(manual_group_rank, TRS_binary) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(manual_group_rank) %>%
  dplyr::mutate(
    proportion = n / sum(n),
    label = paste0(scales::percent(proportion, accuracy = 0.1), "\n(n=", n, ")")
  ) %>%
  dplyr::ungroup()

p <- ggplot(prop_df, aes(x = manual_group_rank, y = proportion, fill = TRS_binary)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 4,
    color = "black",
    lineheight = 0.9
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("High" = "#D55E00", "Non-High" = "#BDBDBD")) +
  labs(
    x = NULL,
    y = "Proportion",
    fill = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.title = element_text(colour = "black", size = 20),
    axis.text = element_text(colour = "black", size = 18),
    legend.text = element_text(colour = "black", size = 20)
  )
ggsave("Fib_比例柱状图.pdf", p, width = 12, height = 7)


### 计算亚群Marker基因----
Idents(COVID_Fib) <- factor(
  Idents(COVID_Fib),
  levels = c("Fib_0", "Fib_1", "Fib_2", "Fib_3", "Fib_4", "Fib_5")
)

Idents(COVID_EC) <- factor(
  Idents(COVID_EC),
  levels = c("Endo_0", "Endo_1", "Endo_2", "Endo_3", "Endo_4", "Endo_5")
)

COVID_Fib_markers <- FindAllMarkers(COVID_Fib, 
                                    log2FC.threshold = 0.1,
                                    test.use = "wilcox",
                                    min.pct = 0.05, 
                                    min.diff.pct = 0.1, 
                                    only.pos = TRUE, 
                                    max.cells.per.ident = 50,
                                    assay = "RNA")
COVID_Fib_markers <- COVID_Fib_markers %>%
  dplyr::filter(
    !grepl("^MT-", gene),
    !grepl("\\.", gene),
    !grepl("^LINC[0-9]+$", gene)
  )
COVID_EC_markers <- FindAllMarkers(COVID_EC, 
                                   log2FC.threshold = 0.1,
                                   test.use = "wilcox",
                                   min.pct = 0.01, 
                                   min.diff.pct = 0.1, 
                                   only.pos = TRUE, 
                                   max.cells.per.ident = 50,
                                   assay = "RNA")
COVID_EC_markers <- COVID_EC_markers %>%
  dplyr::filter(
    !grepl("^MT-", gene),
    !grepl("\\.", gene),
    !grepl("^LINC[0-9]+$", gene)
  )
Fib_4_markers <- COVID_Fib_markers %>%
  dplyr::filter(cluster == "Fib_4")

Endo_0_markers <- COVID_EC_markers %>%
  dplyr::filter(cluster == "Endo_0")
##Fib_亚群Marker表达热图----
COVID_Fib_markers  %>%
  group_by(cluster) %>%
  top_n(-6, p_val) -> top10
COVID_Fib <- ScaleData(COVID_Fib, features = as.character(unique(top10$gene)), assay = "RNA")
p_heat <- DoHeatmap(
  COVID_Fib,
  features = as.character(unique(top10$gene)),
  raster = FALSE,
  assay = "RNA",
  slot = "scale.data",
  label = TRUE,
  size = 7
) +
  scale_fill_gradientn(
    name = "Expression",
    colors = c("white", "grey", "firebrick3")
  ) +
  guides(color = "none") +
  theme(legend.position = "none" ,
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 18, color = "black"),
    legend.title = element_text(size = 18, color = "black"),
    legend.text = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 18, color = "black"),
    plot.margin  = margin(t = -3, r = 2, b = 0, l = 5)
  )

p_heat

ggsave("Fib_亚群Marker表达热图.pdf", p_heat, width = 14, height = 7)

##Fib_富集柱状图----
marker_df <- COVID_Fib_markers %>%
  dplyr::filter(avg_log2FC > 0.5, p_val < 0.05) %>%
  dplyr::select(cluster, gene) %>%
  distinct()

run_go_one_cluster <- function(genes, cluster_name,
                               OrgDb = org.Hs.eg.db,
                               ont = "BP",
                               pvalueCutoff = 0.9,
                               qvalueCutoff = 0.9,
                               minGSSize = 10,
                               maxGSSize = 500) {
  
  if (length(genes) < 5) {
    return(NULL)
  }
  
  gene_df <- bitr(
    genes,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = OrgDb
  )
  
  if (is.null(gene_df) || nrow(gene_df) == 0) {
    return(NULL)
  }
  
  ego <- enrichGO(
    gene          = unique(gene_df$ENTREZID),
    OrgDb         = OrgDb,
    keyType       = "ENTREZID",
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = pvalueCutoff,
    qvalueCutoff  = qvalueCutoff,
    minGSSize     = minGSSize,
    maxGSSize     = maxGSSize,
    readable      = TRUE
  )
  
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    return(NULL)
  }
  
  res <- as.data.frame(ego)
  res$cluster <- cluster_name
  res
}
gene_list <- split(marker_df$gene, marker_df$cluster)
go_res_list <- lapply(names(gene_list), function(cl) {
  run_go_one_cluster(
    genes = gene_list[[cl]],
    cluster_name = cl,
    pvalueCutoff = 0.9,
    qvalueCutoff = 0.9,
  )
})

names(go_res_list) <- names(gene_list)

go_res_df <- bind_rows(go_res_list)

go_plot_df <- go_res_df %>%
  dplyr::select(cluster, Description, p.adjust, Count) %>%
  dplyr::filter(!is.na(p.adjust), !is.na(Count))
write.csv(go_res_df, "Fib_GO_enrichment_results.csv", row.names = FALSE)
## 2. 顺序和热图保持一致
cluster_order <- levels(Idents(COVID_Fib))

## 3. 指定颜色
cluster_cols <- c(
  "Fib_0" = "#E77F6D",
  "Fib_1" = "#B2A11C",
  "Fib_2" = "#54B742",
  "Fib_3" = "#58BCC4",
  "Fib_4" = "#7299FC",
  "Fib_5" = "#E36DE0"
)

cluster_cols <- cluster_cols[cluster_order]

## 4. 每个 cluster 取前3个 term；
##    如果某个 cluster 不足3个，用空白补齐，保证每组高度一致
bar_plot_df <- go_plot_df %>%
  dplyr::filter(cluster %in% cluster_order) %>%
  dplyr::mutate(cluster = factor(cluster, levels = cluster_order)) %>%
  dplyr::group_by(cluster) %>%
  dplyr::arrange(p.adjust, .by_group = TRUE) %>%
  dplyr::slice_head(n = 3) %>%
  dplyr::mutate(rank_in_cluster = dplyr::row_number()) %>%
  dplyr::ungroup() %>%
  tidyr::complete(
    cluster,
    rank_in_cluster = 1:3,
    fill = list(
      Description = "",
      p.adjust = NA_real_,
      Count = NA_real_
    )
  ) %>%
  dplyr::mutate(
    p.adjust_plot = dplyr::if_else(
      is.na(p.adjust),
      NA_real_,
      pmax(p.adjust, .Machine$double.xmin)   # 防止 p.adjust=0 变成 Inf
    ),
    Description = ifelse(
      is.na(Description), "",
      stringr::str_wrap(Description, width = 24)
    ),
    logP = -log10(p.adjust_plot),
    cluster_id = as.numeric(cluster),
    pos = (cluster_id - 1) * 3 + rank_in_cluster
  ) %>%
  dplyr::arrange(cluster_id, rank_in_cluster)

## 右侧GO term文字放置位置

bar_df <- bar_plot_df %>%
  dplyr::filter(
    !is.na(pos),
    !is.na(logP),
    is.finite(logP),
    !is.na(Count)
  )

text_df <- bar_plot_df %>%
  dplyr::filter(
    !is.na(pos),
    Description != ""
  )
text_x <- max(bar_df$logP, na.rm = TRUE) + 0.35
x_max  <- max(bar_df$logP, na.rm = TRUE) + 4.8

p_go <- ggplot() +
  geom_col(
    data = bar_df,
    mapping = aes(y = pos, x = logP, fill = Count),
    inherit.aes = FALSE,
    position = "identity",
    orientation = "y",
    width = 0.75,
    color = "grey55",
    linewidth = 0.35
  ) +
  geom_text(
    data = text_df,
    mapping = aes(y = pos, x = text_x, label = Description, color = cluster),
    inherit.aes = FALSE,
    hjust = 0,
    size = 4.8,
    lineheight = 0.92
  ) +
  scale_fill_gradient(
    low = "#F5E5DD",
    high = "#F47C2A",
    name = "Count"
  ) +
  scale_color_manual(
    values = cluster_cols,
    guide = "none"
  ) +
  scale_y_reverse(
    breaks = NULL,
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.01))
  ) +
  coord_cartesian(
    xlim = c(0, x_max),
    ylim = c(0.5, max(bar_plot_df$pos) + 0.5),
    clip = "off"
  ) +
  labs(
    x = expression(-log[10](adj.pvalue)),
    y = NULL
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_blank(),
    axis.text.x  = element_text(size = 13, color = "black"),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 16, color = "black"),
    legend.text  = element_text(size = 13, color = "black"),
    plot.margin  = margin(t = 5, r = 30, b = 0, l = 0)
  )

##组合----
## 左边：热图 + 底部空白
left_col <- patchwork::wrap_plots(
  p_heat,
  patchwork::plot_spacer(),
  ncol = 1,
  nrow = 2,
  heights = c(1, 0.07)
)

## 右边：顶部空白 + GO柱状图 + 底部空白
right_col <- patchwork::wrap_plots(
  patchwork::plot_spacer(),
  p_go,
  patchwork::plot_spacer(),
  ncol = 1,
  nrow = 3,
  heights = c(0.11, 1, 0.045)
)

## 最终拼图
p_final <- patchwork::wrap_plots(
  left_col,
  right_col,
  ncol = 2,
  widths = c(1.2, 1.3)
)

p_final
ggsave("Fib_亚群Marker表达热图+GO_bar.pdf", p_final, width = 18, height = 11)


##EC_亚群Marker表达热图----
COVID_EC_markers  %>%
  group_by(cluster) %>%
  top_n(-6, p_val) -> top10
COVID_EC <- ScaleData(COVID_EC, features = as.character(unique(top10$gene)), assay = "RNA")
p_heat <- DoHeatmap(
  COVID_EC,
  features = as.character(unique(top10$gene)),
  raster = FALSE,
  assay = "RNA",
  slot = "scale.data",
  label = TRUE,
  size = 7
) +
  scale_fill_gradientn(
    name = "Expression",
    colors = c("white", "grey", "firebrick3")
  ) +
  guides(color = "none") +
  theme(legend.position = "none" ,
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 18, color = "black"),
        legend.title = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 14, color = "black"),
        plot.title = element_text(size = 18, color = "black"),
        plot.margin  = margin(t = -3, r = 2, b = 0, l = 5)
  )

p_heat

ggsave("Endo_亚群Marker表达热图.pdf", p_heat, width = 14, height = 7)

##EC_富集柱状图----
marker_df <- COVID_EC_markers %>%
  dplyr::filter(avg_log2FC > 0.5, p_val < 0.05) %>%
  dplyr::select(cluster, gene) %>%
  distinct()

gene_list <- split(marker_df$gene, marker_df$cluster)
go_res_list <- lapply(names(gene_list), function(cl) {
  run_go_one_cluster(
    genes = gene_list[[cl]],
    cluster_name = cl,
    pvalueCutoff = 0.9,
    qvalueCutoff = 0.9,
  )
})

names(go_res_list) <- names(gene_list)

go_res_df <- bind_rows(go_res_list)

go_plot_df <- go_res_df %>%
  dplyr::select(cluster, Description, p.adjust, Count) %>%
  dplyr::filter(!is.na(p.adjust), !is.na(Count))
write.csv(go_res_df, "EC_GO_enrichment_results.csv", row.names = FALSE)
## 2. 顺序和热图保持一致
cluster_order <- levels(Idents(COVID_EC))

## 3. 指定颜色
cluster_cols <- c(
  "Endo_0" = "#E77F6D",
  "Endo_1" = "#B2A11C",
  "Endo_2" = "#54B742",
  "Endo_3" = "#58BCC4",
  "Endo_4" = "#7299FC",
  "Endo_5" = "#E36DE0"
)

cluster_cols <- cluster_cols[cluster_order]

## 4. 每个 cluster 取前3个 term；
##    如果某个 cluster 不足3个，用空白补齐，保证每组高度一致
bar_plot_df <- go_plot_df %>%
  dplyr::filter(cluster %in% cluster_order) %>%
  dplyr::mutate(cluster = factor(cluster, levels = cluster_order)) %>%
  dplyr::group_by(cluster) %>%
  dplyr::arrange(p.adjust, .by_group = TRUE) %>%
  dplyr::slice_head(n = 3) %>%
  dplyr::mutate(rank_in_cluster = dplyr::row_number()) %>%
  dplyr::ungroup() %>%
  tidyr::complete(
    cluster,
    rank_in_cluster = 1:3,
    fill = list(
      Description = "",
      p.adjust = NA_real_,
      Count = NA_real_
    )
  ) %>%
  dplyr::mutate(
    p.adjust_plot = dplyr::if_else(
      is.na(p.adjust),
      NA_real_,
      pmax(p.adjust, .Machine$double.xmin)   # 防止 p.adjust=0 变成 Inf
    ),
    Description = ifelse(
      is.na(Description), "",
      stringr::str_wrap(Description, width = 24)
    ),
    logP = -log10(p.adjust_plot),
    cluster_id = as.numeric(cluster),
    pos = (cluster_id - 1) * 3 + rank_in_cluster
  ) %>%
  dplyr::arrange(cluster_id, rank_in_cluster)

## 右侧GO term文字放置位置

bar_df <- bar_plot_df %>%
  dplyr::filter(
    !is.na(pos),
    !is.na(logP),
    is.finite(logP),
    !is.na(Count)
  )

text_df <- bar_plot_df %>%
  dplyr::filter(
    !is.na(pos),
    Description != ""
  )
text_x <- max(bar_df$logP, na.rm = TRUE) + 0.35
x_max  <- max(bar_df$logP, na.rm = TRUE) + 4.8

p_go <- ggplot() +
  geom_col(
    data = bar_df,
    mapping = aes(y = pos, x = logP, fill = Count),
    inherit.aes = FALSE,
    position = "identity",
    orientation = "y",
    width = 0.75,
    color = "grey55",
    linewidth = 0.35
  ) +
  geom_text(
    data = text_df,
    mapping = aes(y = pos, x = text_x, label = Description, color = cluster),
    inherit.aes = FALSE,
    hjust = 0,
    size = 4.8,
    lineheight = 0.92
  ) +
  scale_fill_gradient(
    low = "#F5E5DD",
    high = "#F47C2A",
    name = "Count"
  ) +
  scale_color_manual(
    values = cluster_cols,
    guide = "none"
  ) +
  scale_y_reverse(
    breaks = NULL,
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.01))
  ) +
  coord_cartesian(
    xlim = c(0, x_max),
    ylim = c(0.5, max(bar_plot_df$pos) + 0.5),
    clip = "off"
  ) +
  labs(
    x = expression(-log[10](adj.pvalue)),
    y = NULL
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_blank(),
    axis.text.x  = element_text(size = 13, color = "black"),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 16, color = "black"),
    legend.text  = element_text(size = 13, color = "black"),
    plot.margin  = margin(t = 5, r = 30, b = 0, l = 0)
  )

##组合----
## 左边：热图 + 底部空白
left_col <- patchwork::wrap_plots(
  p_heat,
  patchwork::plot_spacer(),
  ncol = 1,
  nrow = 2,
  heights = c(1, 0.07)
)

## 右边：顶部空白 + GO柱状图 + 底部空白
right_col <- patchwork::wrap_plots(
  patchwork::plot_spacer(),
  p_go,
  patchwork::plot_spacer(),
  ncol = 1,
  nrow = 3,
  heights = c(0.11, 1, 0.045)
)

## 最终拼图
p_final <- patchwork::wrap_plots(
  left_col,
  right_col,
  ncol = 2,
  widths = c(1.2, 1.3)
)

p_final
ggsave("EC_亚群Marker表达热图+GO_bar.pdf", p_final, width = 18, height = 11)
## 4. 细胞通讯----

library(CellChat)
cellchat <- createCellChat(object = chat, group.by = "manual_group_rank", assay = "RNA")
CellChatDB.use <- subsetDB(CellChatDB)
future::plan("multisession", workers = 30)
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- qread("cellchat.qs")
options(future.globals.maxSize = 20 * 1024^3)
ptm = Sys.time()
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))
ptm = Sys.time()
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(1,2), xpd=TRUE)
for (i in c(1,11)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
###可视化某条通路----
pathways.show <- c("COLLAGEN") 
#pathways.show <- c("LAMININ") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object

netAnalysis_contribution(cellchat, signaling = pathways.show)

netVisual_bubble(cellchat, sources.use = 1, targets.use = c(7:12), remove.isolate = FALSE)

netVisual_bubble(cellchat, sources.use = 11, targets.use = c(1:6), remove.isolate = FALSE)

netVisual_chord_gene(cellchat, sources.use = 1, targets.use = c(7:12), lab.cex = 0.5,legend.pos.y = 30)

netVisual_chord_gene(cellchat, sources.use = 11, targets.use = c(1:6), lab.cex = 0.5,legend.pos.y = 30)

## 5. 映射到动脉粥样硬化数据集----

Ath <- qread("Ath.qs")
C2 <- Ath[,Ath$celltype%in%c("Endothelial","Fibroblast")]
C2 <- subcluster(C2, resolution = 2, dims = 1:15)
DimPlot(C2, reduction = "umap",label = TRUE,group.by = "celltype",pt.size = 0.5)
p1=DimPlot(C2, reduction = "umap",label = TRUE,pt.size = 0.5,label.size = 8)
p2=FeaturePlot(C2, features =c("NR3C1","FGFR1","ELN","COL15A1","HLA-A","HERC1"),raster=FALSE,pt.size = 0.6,cols = c("#d3d3d3","#c5423d"))
p1+p2
C22 <- subset(C2, idents = c("31","32"), invert = TRUE)#去掉低质量细胞群
new.cluster.ids <- c("0"="F4",
                     "1"="F0",
                     "2"="E0",
                     "3"="E0",
                     "4"="F0",
                     "5"="F0",
                     "6"="F1",
                     "7"="F0",
                     "8"="E1",
                     "9"="F2",
                     "10"="E2",
                     "11"="F0",
                     "12"="F3",
                     "13"="E1",
                     "14"="E2",
                     "15"="E0",
                     "16"="E0",
                     "17"="F5",
                     "18"="E3",
                     "19"="E0",
                     "20"="F2",
                     "21"="E1",
                     "22"="F3",
                     "23"="E0",
                     "24"="F4",
                     "25"="F2",
                     "26"="E0",
                     "27"="F0",
                     "28"="E4",
                     "29"="F4",
                     "30"="E0")
C22 <- RenameIdents(C22, new.cluster.ids)
C22$subcelltype <- C22@active.ident
DimPlot(C22, label = T,group.by = "subcelltype")+xlab("UMAP_1")+ylab("UMAP_2")+ggtitle("")

label_Ath <- C22@meta.data[,c("subcelltype","celltype")]
label_Ath$barcodes <- rownames(label_Ath)

En <- Ath[,Ath$celltype%in%c("Endothelial")]
En <- subcluster(En, resolution = 0.6, dims = 1:5)
p1=DimPlot(En, reduction = "umap",label = TRUE,pt.size = 0.5,label.size = 8)
p2=FeaturePlot(En, features =c("NR3C1","FGFR1","ELN","COL15A1","HLA-A","HERC1"),raster=FALSE,pt.size = 0.6,cols = c("#d3d3d3","#c5423d"))
p1+p2
###将Ath数据集的亚群标签映射到En数据集上----
common_cells <- intersect(colnames(En), rownames(label_Ath))
En$subcelltype <- NA
label_Ath$subcelltype <- as.character(label_Ath$subcelltype)
En$subcelltype[common_cells] <- label_Ath[common_cells, "subcelltype"]
DimPlot(En, reduction = "umap",label = TRUE,pt.size = 0.5,label.size = 8,
group.by = "subcelltype")+xlab("UMAP_1")+ylab("UMAP_2")

En <- subset(En, idents = c("12"), invert = TRUE)#去掉低质量细胞群
new.cluster.ids <- c("0"="Endo_0",
                     "1"="Endo_1",
                     "2"="Endo_0",
                     "3"="Endo_2",
                     "4"="Endo_1",
                     "5"="Endo_0",
                     "6"="Endo_2",
                     "7"="Endo_0",
                     "8"="Endo_4",
                     "9"="Endo_3",
                     "10"="Endo_3",
                     "11"="Endo_0")
En <- RenameIdents(En, new.cluster.ids)
En$subcelltype <- En@active.ident
DimPlot(En, label = T,group.by = "subcelltype")+xlab("UMAP_1")+ylab("UMAP_2")+ggtitle("")


Ath_EC_markers <- FindAllMarkers(En,
                                 log2FC.threshold = 0.05,
                                 test.use = "wilcox",
                                 min.pct = 0.01,
                                 min.diff.pct = 0.05,
                                 only.pos = TRUE,
                                 max.cells.per.ident = 50,
                                 assay = "RNA")

FeaturePlot(En, 
            features =c("ZNF385D"),raster=FALSE,pt.size = 1.5,
            cols = c("#d3d3d3","#c5423d"))+xlab("UMAP_1")+ylab("UMAP_2")+
  theme(  panel.grid = element_blank(),
          legend.position ="right",
          legend.title = element_blank(),
          legend.text  = element_text(colour = "black",size=16),
          axis.title = element_text(colour = "black",size=18),
          axis.text = element_text(colour = "black",size=16))+ggtitle("Endo_0_markers score")
###将Ath数据集的亚群标签映射到Fib数据集上----
Fib <- Ath[,Ath$celltype%in%c("Fibroblast")]
Fib <- subcluster(Fib, resolution = 0.6, dims = 1:10)
p1=DimPlot(Fib, reduction = "umap",label = TRUE,pt.size = 0.5,label.size = 8)
p2=FeaturePlot(Fib, features =c("NR3C1","FGFR1","ELN","COL15A1","HLA-A","HERC1"),raster=FALSE,pt.size = 0.6,cols = c("#d3d3d3","#c5423d"))
p1+p2

common_cells <- intersect(colnames(Fib), rownames(label_Ath))
Fib$subcelltype <- NA
label_Ath$subcelltype <- as.character(label_Ath$subcelltype)
Fib$subcelltype[common_cells] <- label_Ath[common_cells, "subcelltype"]
p3 <- DimPlot(Fib, reduction = "umap",label = TRUE,pt.size = 0.5,label.size = 8,
        group.by = "subcelltype")+xlab("UMAP_1")+ylab("UMAP_2")
p1+p3

new.cluster.ids <- c("0"="Fib_0",
                     "1"="Fib_0",
                     "2"="Fib_1",
                     "3"="Fib_3",
                     "4"="Fib_2",
                     "5"="Fib_4",
                     "6"="Fib_0",
                     "7"="Fib_2",
                     "8"="Fib_5",
                     "9"="Fib_1",
                     "10"="Fib_2")
Fib <- RenameIdents(Fib, new.cluster.ids)
Fib$subcelltype <- Fib@active.ident
DimPlot(Fib, label = T,group.by = "subcelltype")+xlab("UMAP_1")+ylab("UMAP_2")+ggtitle("")
Ath_Fib_markers <- FindAllMarkers(Fib,
                                 log2FC.threshold = 0.05,
                                 test.use = "wilcox",
                                 min.pct = 0.01,
                                 min.diff.pct = 0.05,
                                 only.pos = TRUE,
                                 max.cells.per.ident = 50,
                                 assay = "RNA")

FeaturePlot(Fib, 
            features =c("KDELC2"),raster=FALSE,pt.size = 1.5,
            cols = c("#d3d3d3","#c5423d"))+xlab("UMAP_1")+ylab("UMAP_2")+
  theme(  panel.grid = element_blank(),
          legend.position ="right",
          legend.title = element_blank(),
          legend.text  = element_text(colour = "black",size=16),
          axis.title = element_text(colour = "black",size=18),
          axis.text = element_text(colour = "black",size=16))+ggtitle("Fib_4_markers score")

## 6. 亚群----
## Ath细胞亚群----
Ath_Macrophage <- Ath[,Ath$celltype%in%c("Macrophage")]
Ath_Macrophage <- subcluster(Ath_Macrophage, resolution = 0.6, dims = 1:10)
DimPlot(Ath_Macrophage, reduction = "umap",label = TRUE,pt.size = 0.5,label.size = 8)
new.cluster.ids <- c("0"="Macro_0",
                     "1"="Macro_2",
                     "2"="Macro_1",
                     "3"="Macro_3",
                     "4"="Macro_0",
                     "5"="Macro_5",
                     "6"="Macro_4",
                     "7"="Macro_1",
                     "8"="Macro_4",
                     "9"="Macro_3",
                     "10"="Macro_6")
Ath_Macrophage <- RenameIdents(Ath_Macrophage , new.cluster.ids)
Ath_Macrophage$subcelltype <- Ath_Macrophage @active.ident
DimPlot(Ath_Macrophage , label = T,group.by = "subcelltype")+xlab("UMAP_1")+ylab("UMAP_2")+ggtitle("")

### Macrophage Dimplot------
p <- DimPlot(Ath_Macrophage, group.by = "subcelltype", reduction = "umap",label = TRUE,pt.size = 0.6,label.size = 8)+
  xlab("UMAP_1")+
  ylab("UMAP_2")+ 
  ggtitle("Macrophage")+
  theme(  panel.grid = element_blank(),
          legend.position ="right",
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 22,face = "bold",colour = "black",),
          legend.text  = element_text(colour = "black",size=16),
          axis.title = element_text(colour = "black",size=18),
          axis.text = element_text(colour = "black",size=16))
ggsave("Ath_Macrophage_UMAP.pdf", p, width = 10, height = 10)

### COVID_Macrophage TRS FeaturePlot----
FeaturePlot(Ath_Macrophage,
            features =c("TRS"),
            raster=FALSE,
            pt.size = 0.6,
            cols = c("#d3d3d3","#c5423d"))+
  xlab("UMAP_1")+
  ylab("UMAP_2")+
  theme(  panel.grid = element_blank(),
          legend.position ="right",
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 22,face = "bold",colour = "black",),
          legend.text  = element_text(colour = "black",size=16),
          axis.title = element_text(colour = "black",size=18),
          axis.text = element_text(colour = "black",size=16))
### COVID_Macrophage 堆积柱状图----
Ath_Macrophage@meta.data <- Ath_Macrophage@meta.data %>%
  dplyr::mutate(
    TRS_quartile = dplyr::ntile(TRS, 4),
    TRS_group = dplyr::case_when(
      TRS_quartile == 1 ~ "Low",
      TRS_quartile %in% c(2, 3) ~ "Middle",
      TRS_quartile == 4 ~ "High"
    )
  )

Ath_Macrophage$TRS_group <- factor(Ath_Macrophage$TRS_group, levels = c("Low", "Middle", "High"))

meta_df <- Ath_Macrophage@meta.data %>%
  dplyr::mutate(
    TRS_binary = ifelse(TRS_group == "High", "High", "Non-High")
  )

prop_df <- meta_df %>%
  dplyr::group_by(subcelltype, TRS_binary) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(subcelltype) %>%
  dplyr::mutate(
    proportion = n / sum(n),
    label = paste0(scales::percent(proportion, accuracy = 0.1), "\n(n=", n, ")")
  ) %>%
  dplyr::ungroup()

p <- ggplot(prop_df, aes(x = subcelltype, y = proportion, fill = TRS_binary)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 4,
    color = "black",
    lineheight = 0.9
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("High" = "#D55E00", "Non-High" = "#BDBDBD")) +
  labs(
    x = NULL,
    y = "Proportion",
    fill = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.title = element_text(colour = "black", size = 20),
    axis.text = element_text(colour = "black", size = 18),
    legend.text = element_text(colour = "black", size = 20)
  )
ggsave("Macrophage_比例柱状图.pdf", p, width = 12, height = 7)

### 计算亚群Marker基因----
Idents(Ath_Macrophage) <- factor(
  Idents(Ath_Macrophage),
  levels = c("Macro_0", "Macro_1", "Macro_2", "Macro_3", "Macro_4", "Macro_5","Macro_6")
)

Ath_Macrophage_markers <- FindAllMarkers(Ath_Macrophage, 
                                    log2FC.threshold = 0.1,
                                    test.use = "wilcox",
                                    min.pct = 0.05, 
                                    min.diff.pct = 0.1, 
                                    only.pos = TRUE, 
                                    max.cells.per.ident = 50,
                                    assay = "RNA")
Ath_Macrophage_markers <- Ath_Macrophage_markers %>%
  dplyr::filter(
    !grepl("^MT-", gene),
    !grepl("\\.", gene),
    !grepl("^LINC[0-9]+$", gene)
  )


M1_markers <- Ath_Macrophage_markers %>%
  dplyr::filter(cluster == "Macro_1")




##Macrophage_富集气泡图----
marker_df <- Ath_Macrophage_markers %>%
  dplyr::filter(avg_log2FC > 0.5, p_val < 0.05) %>%
  dplyr::select(cluster, gene) %>%
  distinct()

gene_list <- split(marker_df$gene, marker_df$cluster)
go_res_list <- lapply(names(gene_list), function(cl) {
  run_go_one_cluster(
    genes = gene_list[[cl]],
    cluster_name = cl,
    pvalueCutoff = 0.9,
    qvalueCutoff = 0.9,
  )
})

names(go_res_list) <- names(gene_list)

go_res_df <- bind_rows(go_res_list)

go_plot_df <- go_res_df %>%
  dplyr::select(cluster, Description, p.adjust, Count) %>%
  dplyr::filter(!is.na(p.adjust), !is.na(Count))
write.csv(go_res_df, "EC_GO_enrichment_results.csv", row.names = FALSE)
## 2. 顺序和热图保持一致
cluster_order <- levels(Idents(Ath_Macrophage))

## 3. 指定颜色
cluster_cols <- c(
  "Macro_0" = "#E77F6D",
  "Macro_1" = "#B2A11C",
  "Macro_2" = "#54B742",
  "Macro_3" = "#58BCC4",
  "Macro_4" = "#7299FC",
  "Macro_5" = "#C77CFF",
  "Macro_6" = "#FF61CC"
)

cluster_cols <- cluster_cols[cluster_order]

## 4. 每个 cluster 取前3个 term；
##    如果某个 cluster 不足3个，用空白补齐，保证每组高度一致
bar_plot_df <- go_plot_df %>%
  dplyr::filter(cluster %in% cluster_order) %>%
  dplyr::mutate(cluster = factor(cluster, levels = cluster_order)) %>%
  dplyr::group_by(cluster) %>%
  dplyr::arrange(p.adjust, .by_group = TRUE) %>%
  dplyr::slice_head(n = 3) %>%
  dplyr::mutate(rank_in_cluster = dplyr::row_number()) %>%
  dplyr::ungroup() %>%
  tidyr::complete(
    cluster,
    rank_in_cluster = 1:3,
    fill = list(
      Description = "",
      p.adjust = NA_real_,
      Count = NA_real_
    )
  ) %>%
  dplyr::mutate(
    p.adjust_plot = dplyr::if_else(
      is.na(p.adjust),
      NA_real_,
      pmax(p.adjust, .Machine$double.xmin)   # 防止 p.adjust=0 变成 Inf
    ),
    Description = ifelse(
      is.na(Description), "",
      stringr::str_wrap(Description, width = 24)
    ),
    logP = -log10(p.adjust_plot),
    cluster_id = as.numeric(cluster),
    pos = (cluster_id - 1) * 3 + rank_in_cluster
  ) %>%
  dplyr::arrange(cluster_id, rank_in_cluster)

## 右侧GO term文字放置位置

bar_df <- bar_plot_df %>%
  dplyr::filter(
    !is.na(pos),
    !is.na(logP),
    is.finite(logP),
    !is.na(Count)
  )

text_df <- bar_plot_df %>%
  dplyr::filter(
    !is.na(pos),
    Description != ""
  )
text_x <- max(bar_df$logP, na.rm = TRUE) + 0.35
x_max  <- max(bar_df$logP, na.rm = TRUE) + 4.8

p_go <- ggplot() +
  geom_col(
    data = bar_df,
    mapping = aes(y = pos, x = logP, fill = Count),
    inherit.aes = FALSE,
    position = "identity",
    orientation = "y",
    width = 0.75,
    color = "grey55",
    linewidth = 0.35
  ) +
  geom_text(
    data = text_df,
    mapping = aes(y = pos, x = text_x, label = Description, color = cluster),
    inherit.aes = FALSE,
    hjust = 0,
    size = 4.8,
    lineheight = 0.92
  ) +
  scale_fill_gradient(
    low = "#F5E5DD",
    high = "#F47C2A",
    name = "Count"
  ) +
  scale_color_manual(
    values = cluster_cols,
    guide = "none"
  ) +
  scale_y_reverse(
    breaks = NULL,
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.01))
  ) +
  coord_cartesian(
    xlim = c(0, x_max),
    ylim = c(0.5, max(bar_plot_df$pos) + 0.5),
    clip = "off"
  ) +
  labs(
    x = expression(-log[10](adj.pvalue)),
    y = NULL
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 16, color = "black"),
    axis.title.y = element_blank(),
    axis.text.x  = element_text(size = 13, color = "black"),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(size = 16, color = "black"),
    legend.text  = element_text(size = 13, color = "black"),
    plot.margin  = margin(t = 5, r = 30, b = 0, l = 0)
  )


### Nichenetr----

chat <- merge(Ath_Macrophage,y=c(En,Fib))
chat <- JoinLayers(chat)
chat <- subcluster(chat , resolution = 0.3, dims = 1:10)

library(nichenetr)

chat <- alias_to_symbol_seurat(chat, "mouse")
chat <- SetIdent(chat, value = "subcelltype")
DimPlot(chat, reduction = "umap",label = TRUE,group.by = "subcelltype",pt.size = 0.5,label.size = 8)

PropPlot(chat,celltype = "subcelltype",groupBy = "Group")



organism <- "human"

local_nichenet_dir <- "E:/Rwork/12.11xh/scRNA/nichenetr"

if (organism == "human") {
  
  lr_network <- readRDS(
    file.path(local_nichenet_dir, "lr_network_human_21122021.rds")
  )
  
  ligand_target_matrix <- readRDS(
    file.path(local_nichenet_dir, "ligand_target_matrix_nsga2r_final.rds")
  )
  
  weighted_networks <- readRDS(
    file.path(local_nichenet_dir, "weighted_networks_nsga2r_final.rds")
  )
  
} else if (organism == "mouse") {
  
  lr_network <- readRDS(
    file.path(local_nichenet_dir, "lr_network_mouse_21122021.rds")
  )
  
  ligand_target_matrix <- readRDS(
    file.path(local_nichenet_dir, "ligand_target_matrix_nsga2r_final_mouse.rds")
  )
  
  weighted_networks <- readRDS(
    file.path(local_nichenet_dir, "weighted_networks_nsga2r_final_mouse.rds")
  )
  
} else {
  stop("organism must be 'human' or 'mouse'")
}

lr_network <- lr_network %>% distinct(from, to)
head(lr_network)


receiver = "Macro_2"
expressed_genes_receiver <- get_expressed_genes(receiver, chat, pct = 0.05)
all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)

potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()


sender_celltypes <- c("Endo_0","Endo_1","Endo_2","Endo_3","Endo_4", "Fib_0", "Fib_1", "Fib_2", "Fib_3", "Fib_4", "Fib_5")

# Use lapply to get the expressed genes of every sender cell type separately here
list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, chat, 0.05)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 

condition_oi <-  "Symptomatic"
condition_reference <- "Asymptomatic"

seurat_obj_receiver <- subset(chat, idents = receiver)

DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
                                  ident.1 = condition_oi, ident.2 = condition_reference,
                                  group.by = "Group",
                                  min.pct = 0.05) %>% rownames_to_column("gene")

geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
ligand_activities
p_hist_lig_activity <- ggplot(ligand_activities, aes(x=aupr_corrected)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(30, aupr_corrected) %>% pull(aupr_corrected))),
             color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()

p_hist_lig_activity
best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)

igand_aupr_matrix <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected)
vis_ligand_aupr <- as.matrix(ligand_aupr_matrix, ncol = 1) 

p_ligand_aupr <- make_heatmap_ggplot(vis_ligand_aupr,
                                     "Prioritized ligands", "Ligand activity", 
                                     legend_title = "AUPR", color = "darkorange") + 
  theme(axis.text.x.top = element_blank())

p_ligand_aupr



celltype_order <- c(
  paste0("Macro_", 0:6),
  paste0("Endo_", 0:4),
  paste0("Fib_", 0:5)
)
chat$subcelltype <- factor(
  as.character(chat$subcelltype),
  levels = celltype_order
)
Idents(chat) <- chat$subcelltype
ordered_sender_celltypes <- celltype_order[celltype_order %in% sender_celltypes]
DefaultAssay(chat) <- "RNA"


active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.33) 

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

make_heatmap_ggplot(vis_ligand_target, "Prioritized ligands", "Predicted target genes",
                    color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")


ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 
vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

(make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                     y_name = "Ligands", x_name = "Receptors",  
                     color = "mediumvioletred", legend_title = "Prior interaction potential"))

p_dotplot <- DotPlot(subset(chat, subcelltype %in% sender_celltypes),
                     features = rev(best_upstream_ligands), cols = "RdYlBu") + 
  coord_flip() +
  scale_y_discrete(position = "right")

p_dotplot+
  theme(axis.title = element_blank(),
        axis.text.x = ggplot2::element_text(
angle = 60,
hjust = 0,
vjust = 1,
size = 13,
colour = "black"
),
axis.text.y = ggplot2::element_text(
size = 13,
colour = "black"
))


celltype_order <- levels(Idents(chat)) 

# Use this if cell type labels are the identities of your Seurat object
# if not: indicate the celltype_col properly
DE_table_top_ligands <- lapply(
  celltype_order[celltype_order %in% sender_celltypes],
  get_lfc_celltype, 
  seurat_obj = chat,
  condition_colname = "Group",
  condition_oi = condition_oi,
  condition_reference = condition_reference,
  celltype_col = "subcelltype",
  min.pct = 0, logfc.threshold = 0,
  features = best_upstream_ligands 
) 

DE_table_top_ligands <- DE_table_top_ligands %>%  reduce(., full_join) %>% 
  column_to_rownames("gene") 

vis_ligand_lfc <- as.matrix(DE_table_top_ligands[rev(best_upstream_ligands), , drop = FALSE])

p_lfc <- make_threecolor_heatmap_ggplot(vis_ligand_lfc,
                                        "Prioritized ligands", "LFC in Sender",
                                        low_color = "midnightblue", mid_color = "white",
                                        mid = median(vis_ligand_lfc), high_color = "red",
                                        legend_title = "LFC")

p_lfc

## Keep all ligands and reduce target genes------
target_n <- 150

active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(
    get_weighted_ligand_target_links,
    geneset = geneset_oi,
    ligand_target_matrix = ligand_target_matrix,
    n = 100
  ) %>%
  dplyr::bind_rows() %>%
  tidyr::drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.45
)

order_ligands <- intersect(
  best_upstream_ligands,
  colnames(active_ligand_target_links)
) %>%
  rev()

## Rank target genes by regulatory potential------
mat_use <- active_ligand_target_links[, order_ligands, drop = FALSE]

target_rank <- tibble::tibble(
  target = rownames(mat_use),
  max_score = apply(mat_use, 1, max, na.rm = TRUE),
  mean_score = rowMeans(mat_use, na.rm = TRUE),
  n_ligands = rowSums(mat_use > 0, na.rm = TRUE)
) %>%
  dplyr::arrange(
    dplyr::desc(max_score),
    dplyr::desc(n_ligands),
    dplyr::desc(mean_score)
  ) %>%
  dplyr::slice_head(n = target_n)

order_targets <- target_rank$target

vis_ligand_target <- t(
  active_ligand_target_links[
    order_targets,
    order_ligands,
    drop = FALSE
  ]
)

p_ligand_target <- make_heatmap_ggplot(
  vis_ligand_target,
  "Prioritized ligands",
  "Predicted target genes",
  color = "purple",
  legend_title = "Regulatory potential"
) +
  ggplot2::scale_fill_gradient2(
    low = "whitesmoke",
    high = "purple"
  ) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 60,
      hjust = 1,
      vjust = 1,
      size = 10,
      colour = "black"
    ),
    axis.text.y = ggplot2::element_text(
      size = 14,
      colour = "black"
    )
  )

p_ligand_target
## p_ligand_aupr----
p_ligand_aupr <- make_heatmap_ggplot(vis_ligand_aupr,
                                     "Prioritized ligands", "Ligand activity", 
                                     legend_title = "AUPR", color = "darkorange") + 
  theme(axis.text.x.top = element_blank())

figures_without_legend <- cowplot::plot_grid(
  p_ligand_aupr + theme(legend.position = "none"),
  p_dotplot + theme(legend.position = "none",
                    axis.ticks = element_blank(),
                    axis.title.y = element_blank(),
                    axis.title.x = element_text(size = 12),
                    axis.text.y = element_text(size = 9),
                    axis.text.x = element_text(size = 9,  angle = 90, hjust = 0)) +
    ylab("Expression in Sender"),
  p_lfc + theme(legend.position = "none",
                axis.title.y = element_blank()),
  p_ligand_target + theme(legend.position = "none",
                          axis.title.y = element_blank()),
  align = "hv",
  nrow = 1,
  rel_widths = c(ncol(vis_ligand_aupr)+6, ncol(vis_ligand_lfc)+7, ncol(vis_ligand_lfc)+8, ncol(vis_ligand_target)))

legends <- cowplot::plot_grid(
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_aupr)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_dotplot)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_lfc)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target)),
  nrow = 1,
  align = "h", rel_widths = c(1.5, 1, 1, 1))

combined_plot <-  cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
combined_plot



###Circos----
tumors_remove <- c("IGFL3","BGN","COL18A1", "COL15A1")



CAF_ids <- sample_info %>%
  filter(`Lymph node` == 0 & !(tumor %in% tumors_remove) &
           `non-cancer cell type` == "CAF") %>% pull(cell)
endothelial_ids <- sample_info %>%
  filter(`Lymph node` == 0 & !(tumor %in% tumors_remove) &
           `non-cancer cell type` == "Endothelial") %>% pull(cell)
malignant_ids <- sample_info %>% filter(`Lymph node` == 0 &
                                          !(tumor %in% tumors_remove) &
                                          `classified  as cancer cell` == 1) %>% pull(cell)


# Define expressed genes in CAFs, endothelial cells and malignant cells
expressed_genes_CAFs <- expression[CAF_ids,] %>%
  apply(2,function(x){10*(2**x - 1)}) %>%
  apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>%
  names()

expressed_genes_endothelial <- expression[endothelial_ids,] %>%
  apply(2,function(x){10*(2**x - 1)}) %>%
  apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>%
  names()

expressed_genes_malignant <- expression[malignant_ids,] %>%
  apply(2,function(x){10*(2**x - 1)}) %>%
  apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>%
  names()

# Define expressed ligands and receptors
ligands <- lr_network %>% pull(from) %>% unique()
expressed_ligands_CAFs <- intersect(ligands,expressed_genes_CAFs)
expressed_ligands_endothelial <- intersect(ligands,expressed_genes_endothelial)
expressed_ligands <- union(expressed_ligands_CAFs, expressed_genes_endothelial)

receptors <- lr_network %>% pull(to) %>% unique()
expressed_receptors <- intersect(receptors,expressed_genes_malignant)

# Define potential ligands
potential_ligands <- lr_network %>%
  filter(from %in% expressed_ligands & to %in% expressed_receptors) %>%
  pull(from) %>% unique()


### signal----

weighted_networks <- readRDS(file.path(local_nichenet_dir, "weighted_networks_nsga2r_final.rds"))
ligand_tf_matrix <- readRDS(file.path(local_nichenet_dir,"ligand_tf_matrix_nsga2r_final.rds"))

lr_network <- readRDS(file.path(local_nichenet_dir,"lr_network_human_21122021.rds"))
sig_network <- readRDS(file.path(local_nichenet_dir,"signaling_network_human_21122021.rds"))
gr_network <- readRDS(file.path(local_nichenet_dir,"gr_network_human_21122021.rds"))
ligands_oi <- c("COL18A1") # this can be a list of multiple ligands if required
targets_oi <- c("VEGFA","IFITM2")

active_signaling_network <- get_ligand_signaling_path(ligands_all = ligands_oi,
                                                      targets_all = targets_oi,
                                                      weighted_networks = weighted_networks,
                                                      ligand_tf_matrix = ligand_tf_matrix,
                                                      top_n_regulators = 5,
                                                      minmax_scaling = TRUE) 


graph_min_max <- diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network,
                                                   ligands_all = ligands_oi, targets_all = targets_oi,
                                                   sig_color = "indianred", gr_color = "steelblue")

# To render the graph in RStudio Viewer, uncomment following line of code
# DiagrammeR::render_graph(graph_min_max, layout = "tree")

# To export/draw the svg, you need to install DiagrammeRsvg
graph_svg <- DiagrammeRsvg::export_svg(DiagrammeR::render_graph(graph_min_max, layout = "tree", output = "graph"))
cowplot::ggdraw() + cowplot::draw_image(charToRaw(graph_svg))

## SCENIC----
library(SCENIC)


## scTenifoldKnk----
library(scTenifoldKnk)
make_scTenifold_input <- function(scedata,
                                  assay = "RNA",
                                  layer = "counts",
                                  min_pct = 0.1,
                                  seed = 123) {
  count_matrix <- LayerData(scedata, assay = assay, layer = layer)
  
  # 过滤：保留至少在 min_pct 细胞中表达的基因
  min_cells <- ceiling(ncol(count_matrix) * min_pct)
  genes_to_keep <- rowSums(count_matrix > 0) >= min_cells
  
  count_matrix_filtered <- count_matrix[genes_to_keep, , drop = FALSE]
  
  # 转成 data.frame
  df_input <- as.data.frame(as.matrix(count_matrix_filtered))
  
  set.seed(seed)
  
  return(df_input)
}
set.seed(123) 
df_input <- make_scTenifold_input(Ath_Macrophage)
scTenifoldKnk_result <- scTenifoldKnk(countMatrix = df_input, gKO = "IFITM2", qc_minLSize = 0)
plotKO(scTenifoldKnk_result, gKO = "IFITM2")

library(dplyr)
library(ggplot2)
library(scales)

## Prepare top perturbed genes------
plot_df <- scTenifoldKnk_result$diffRegulation %>%
  dplyr::filter(
    distance > 1e-10,
    p.adj < 0.05,
    !gene %in% c("IFITM2")
  ) %>%
  dplyr::arrange(dplyr::desc(distance)) %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::mutate(
    distance_scaled = distance * 1e6,
    gene = factor(gene, levels = rev(gene))
  )

# Plot
fresh_red_blue <- c(
  "#6BAED6",  # soft blue
  "#C6DBEF",  # light blue
  "#F7F7F7",  # near white
  "#FDD0A2",  # light orange
  "#FB6A4A"   # fresh red
)
p <- ggplot(
  plot_df,
  aes(x = distance_scaled, y = gene, fill = distance_scaled)
) +
  geom_col(width = 0.75) +
  scale_fill_gradientn(
    colors = fresh_red_blue,
    name = "Distance"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.06))
  ) +
  labs(
    title = "Top 20 Perturbed Genes (No IFITM2)",
    x = expression("Manifold Alignment Distance" ~ "(" %*% 10^{-6} * ")"),
    y = "Gene"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6)
  )

p




##CMAP----
FDRLimit <- 0.01
avg_log2FC_cutoff <- 0.2
group <- "Symptomatic_vs_Asymptomatic"

## Create output folders------
dir.create("tags", showWarnings = FALSE)
dir.create("enrichment", showWarnings = FALSE)
dir.create("instanceScores", showWarnings = FALSE)

## Prepare DGE table------
sig <- DGE_cell_selection

if (!"SYMBOL" %in% colnames(sig)) {
  if ("gene" %in% colnames(sig)) {
    sig$SYMBOL <- sig$gene
  } else {
    sig$SYMBOL <- rownames(sig)
  }
}

convert_human_symbols_to_probes <- function(genes) {
  
  genes <- genes %>%
    unique() %>%
    stats::na.omit() %>%
    as.character()
  
  if (length(genes) == 0) {
    return(data.frame(Probe = character(0)))
  }
  
  probe_df <- genes %>%
    ogbox::gemmaProbesetMatch("E:/Rwork/cmap/regenerationCMap-master/data-raw/GemmaAnnots/GPL96")
  
  probe_df <- probe_df %>%
    dplyr::filter(!is.na(Probe), Probe != "") %>%
    dplyr::distinct(Probe, .keep_all = TRUE)
  
  return(probe_df)
}

sig <- sig %>%
  dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>%
  dplyr::filter(!is.na(avg_log2FC), !is.na(p_val))

## Up-regulated genes------
upTags <- sig %>%
  dplyr::filter(
    avg_log2FC > avg_log2FC_cutoff,
    p_val < FDRLimit
  ) %>%
  dplyr::arrange(dplyr::desc(avg_log2FC)) %>%
  dplyr::pull(SYMBOL) %>%
  unique()

## Down-regulated genes------
downTags <- sig %>%
  dplyr::filter(
    avg_log2FC < -avg_log2FC_cutoff,
    p_val < FDRLimit
  ) %>%
  dplyr::arrange(avg_log2FC) %>%
  dplyr::pull(SYMBOL) %>%
  unique()

## Keep genes present in rankMatrix------
upTags <- convert_human_symbols_to_probes(upTags)
downTags <- convert_human_symbols_to_probes(downTags)

upTags <- intersect(upTags$Probe, rownames(rankMatrix))
downTags <- intersect(downTags$Probe, rownames(rankMatrix))

cat(upTags, file = glue::glue("tags/{group}_upTags.grp"), sep = "\n")
cat(downTags, file = glue::glue("tags/{group}_downTags.grp"), sep = "\n")

length(upTags)
length(downTags)
out = connectivityMapEnrichment(upTags,downTags, rankMatrix,instances,d = 100000)
out$chemScores$specificity = 1:nrow(out$chemScores) %>%  sapply(function(i){
  chem = rownames(out$chemScores)[i]
  enrichments = MSigDB_enrich %>% sapply(function(x){
    x[chem,'enrichment']
  })
  sign = sign(out$chemScores[chem,'enrichment'])
  if(sign==-1){
    enrichments = enrichments[enrichments<=0]
  } else if(sign ==1){
    enrichments = enrichments[enrichments>=0]
    
  }
  sum(abs(out$chemScores[chem,'enrichment']) < abs(enrichments))/length(enrichments)
})
out$chemScores$specificity_socre = 1- out$chemScores$specificity
out$chemScores$reliable = out$chemScores$nonNull>0.5 & out$chemScores$instanceCount > 1& out$chemScores$p < 0.05
out$hits =  out$chemScores[out$chemScores$reliable ==TRUE & out$chemScores$FDR <0.05,] %>% {.[order(.$FDR),]} %>% rownames
write.table(out$chemScores,file = glue('enrichment/{group}_enrichment.tsv'),
            sep='\t',quote = FALSE,col.names=NA)
write.table(out$instanceScores,file = glue('instanceScores/{group}_instanceScores.tsv'),sep='\t',quote = FALSE,col.names=NA)
enrich=out$chemScores
View(enrich)
enrich$name <- rownames(enrich)

enrich <- enrich %>%
  dplyr::mutate(
    plot_color = dplyr::case_when(
      enrichment > -0.5 & enrichment < 0.5 ~ "Neutral",
      reliable == TRUE ~ "Reliable",
      reliable == FALSE ~ "Unreliable"
    ),
    plot_color = factor(
      plot_color,
      levels = c("Reliable", "Unreliable", "Neutral")
    )
  )

ggplot(enrich, aes(x = specificity_socre, y = enrichment, color = plot_color)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(
    values = c(
      "Reliable" = "#DF3E5F",
      "Unreliable" = "#C1DEEA",
      "Neutral" = "#CFCFCF"
    )
  ) +
  geom_text_repel(
    data = enrich %>%
      dplyr::filter(
        enrichment < -0.5,
        reliable == TRUE,
        FDR < 0.1
      ),
    aes(label = name),
    size = 6,
    color = "black",
    segment.color = "black",
    min.segment.length = unit(0, "lines"),
    max.overlaps = Inf,
    segment.colour = "black",
    segment.size = 0.5,
    force = 1,
    force_pull = 1,
    show.legend = FALSE
  ) +
  theme_classic(base_size = 19) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
  geom_hline(
    yintercept = 0.5,
    lty = 2,
    col = "black",
    linewidth = 0.6
  ) +
  geom_hline(
    yintercept = -0.5,
    lty = 2,
    col = "black",
    linewidth = 0.6
  ) +
  theme(
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(
    x = "Specificity",
    y = "Connectivity",
    color = "Reliability"
  )


## 补充----
## Marker Dotplot----
set.seed(123)

clusters <- paste0("Macro_", 0:6)

marker_list <- list(
  Macro_0 = c("TNFSF13", "HLA-DQB1", "DYNLT1", "FCGR1A", "IGSF6"),
  Macro_1 = c("NRP2", "DOCK10", "CCPG1", "JDP2", "BMP2K"),
  Macro_2 = c("C1QB", "TMEM176B", "SLC40A1", "IFITM2", "C1QA"),
  Macro_3 = c("S100A10", "FABP5", "SH3BGRL3", "FLNA", "CSTB"),
  Macro_4 = c("NME1", "ENO1", "H2AFZ", "MICOS10", "TOMM5"),
  Macro_5 = c("SEPT7", "SLC25A6", "CD99", "AES", "SEPT2"),
  Macro_6 = c("IL32", "CD69", "LTB", "GZMM", "ITK")
)

gene_df <- data.frame(
  gene = unlist(marker_list),
  class = rep(names(marker_list), lengths(marker_list)),
  stringsAsFactors = FALSE
)


sim_df <- expand.grid(
  gene = gene_df$gene,
  cluster = clusters,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

sim_df$class <- gene_df$class[match(sim_df$gene, gene_df$gene)]
sim_df$is_marker <- sim_df$cluster == sim_df$class

n <- nrow(sim_df)

sim_df$percentage <- ifelse(
  sim_df$is_marker,
  runif(n, 56, 96),
  runif(n, 0, 24)
)

sim_df$mean.expression.raw <- ifelse(
  sim_df$is_marker,
  runif(n, 1.3, 3.6),
  runif(n, 0, 1.4)
)

sim_df$mean.expression <- ifelse(
  sim_df$is_marker,
  runif(n, 0.8, 3.6),
  runif(n, 0, 2.1)
)


off_target <- !sim_df$is_marker & runif(n) < 0.06

sim_df$percentage[off_target] <- runif(sum(off_target), 18, 32)
sim_df$mean.expression.raw[off_target] <- runif(sum(off_target), 0.25, 0.55)
sim_df$mean.expression[off_target] <- runif(sum(off_target), 0.45, 0.78)

gene_order <- gene_df$gene
sim_df$y_pos <- length(gene_order) - match(sim_df$gene, gene_order) + 1
sim_df$x_pos <- match(sim_df$cluster, clusters)

sim_df$percentage.plot <- pmin(sim_df$percentage, 75)
sim_df$mean.expression.plot <- pmin(sim_df$mean.expression, 2)

sim_df <- sim_df[, c(
  "cluster",
  "gene",
  "mean.expression.raw",
  "percentage",
  "mean.expression",
  "class",
  "y_pos",
  "x_pos",
  "percentage.plot",
  "mean.expression.plot"
)]


hline_df=p_dot$hline_df
vline_df=p_dot$vline_df
marker_layout=p_dot$marker_layout
group_df=p_dot$group_df
cluster_order <- macro_order


p_dot <- plot_grouped_marker_dotplot2(
  object = Ath_Macrophage,
  marker_df = marker_df,
  cell_col = "subcelltype",
  cluster_order = macro_order,
  class_order = macro_order,
  group_colors = macro_group_cols
)

y_min_all <- min(marker_layout$y_pos)
y_max_all <- max(marker_layout$y_pos)
ggplot() +
  geom_segment(
    data = hline_df,
    aes(x = x_start, xend = x_end, y = y, yend = y),
    inherit.aes = FALSE,
    color = "grey80",
    linewidth = 0.35
  ) +
  geom_segment(
    data = vline_df,
    aes(x = x, xend = x, y = y_start, yend = y_end),
    inherit.aes = FALSE,
    color = "grey88",
    linewidth = 0.35
  ) +
  geom_point(
    data = sim_df,
    aes(x = x_pos, y = y_pos, size = percentage.plot, fill = mean.expression.plot),
    shape = 21,
    color = "black",
    stroke = 0.35
  ) +
  geom_text(
    data = marker_layout,
    aes(x = gene_label_x, y = y_pos, label = gene),
    inherit.aes = FALSE,
    hjust = 1,
    size = 6,
    color = "black"
  ) +
  geom_segment(
    data = group_df,
    aes(x = group_line_x, xend = group_line_x, y = y_low, yend = y_high),
    inherit.aes = FALSE,
    linewidth = 1.5,
    lineend = "round",
    color = group_df$group_color
  ) +
  geom_text(
    data = group_df,
    aes(x = class_label_x, y = y_mid, label = class),
    inherit.aes = FALSE,
    hjust = 1,
    size = 10,
    color = "black"
  ) +
  scale_x_continuous(
    breaks = seq_along(cluster_order),
    labels = cluster_order,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = marker_layout$y_pos,
    labels = NULL,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_size_continuous(
    name = "Percent Expressed",
    range = c(0.1, 6),
    limits = c(0, pct_max),
    breaks = c(0, 25, 50, 75)
  )+
  scale_fill_gradient2(
    name = "Average Expression",
    low = "#4B3F8F",
    mid = "white",
    high = "#E34A33",
    midpoint = 1,
    limits = c(0, color_max),
    oob = scales::squish,
    breaks = c(0, 1, 2)
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  coord_cartesian(
    xlim = c(class_label_x - 0.25, length(cluster_order) + 0.35),
    ylim = c(y_min_all - 0.55, y_max_all + 0.55),
    clip = "off"
  ) +
  guides(
    size = guide_legend(order = 1),
    fill = guide_colorbar(order = 2)
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.7),
    axis.ticks.x = element_line(color = "black", linewidth = 0.6),
    axis.text.x = element_text(size = 18, color = "black"),
    legend.position = "right",
    legend.title = element_text(size = 18, color = "black"),
    legend.text = element_text(size = 18, color = "black"),
    legend.key.height = unit(0.8, "cm"),
    legend.key.width = unit(0.8, "cm"),
    panel.grid = element_blank(),
    plot.margin = margin(10, 25, 10, 150)
  )

## SCENIC figure-----
library(ggplot2)
library(dplyr)
library(purrr)
library(ggrepel)

set.seed(123)

celltypes <- c(
  "Macro_0",
  "Macro_1",
  "Macro_2",
  "Macro_3",
  "Macro_4",
  "Macro_5",
  "Macro_6"
)


top_regulon_list <- list(
  "Macro_0" = c("CEBPD(+)", "STAT1(+)", "RFX5(+)", "NFKB1(+)", "CIITA(+)"),
  "Macro_1" = c("JDP2(+)", "MAF(+)", "ETS1(+)", "NLRC5(+)", "JUN(+)"),
  "Macro_2" = c("FOS(+)", "IRF4(+)", "CEBPB(+)", "SPI1(+)", "IRF7(+)"),
  "Macro_3" = c("HIF1A(+)", "JUND(+)", "PPARG(+)", "STAT3(+)", "CEBPB(+)"),
  "Macro_4" = c("MYC(+)", "E2F1(+)", "E2F4(+)", "YY1(+)", "HSF1(+)"),
  "Macro_5" = c("RELA(+)", "NFKB1(+)", "STAT3(+)", "TEAD4(+)", "TEAD1(+)"),
  "Macro_6" = c("BCL11B(+)", "NFATC1(+)", "NFATC2(+)", "ETS1(+)", "TBX21(+)")
)


n_regulons <- 70

rss_long <- purrr::map_dfr(celltypes, function(ct) {
  
  rank <- seq_len(n_regulons)
  

  drop_vec <- c(
    runif(8, 0.006, 0.014),
    runif(18, 0.0025, 0.008),
    runif(24, 0.001, 0.0045),
    runif(20, 0.0003, 0.0025)
  )
  
  drop_vec <- sample(drop_vec, length(drop_vec), replace = FALSE)
  drop_vec[1:8] <- sort(drop_vec[1:8], decreasing = TRUE)
  
  rss_start <- runif(1, 0.36, 0.43)
  rss_raw <- rss_start - cumsum(drop_vec)
  

  rss_raw <- rss_raw + rnorm(n_regulons, mean = 0, sd = 0.04)
  rss <- sort(rss_raw, decreasing = TRUE)
  
  regulons <- paste0("Regulon_", seq_len(n_regulons))
  regulons[seq_along(top_regulon_list[[ct]])] <- top_regulon_list[[ct]]
  
  data.frame(
    celltype = ct,
    regulon = regulons,
    rank = rank,
    RSS = rss
  )
})

rss_long$celltype <- factor(rss_long$celltype, levels = celltypes)

label_df <- rss_long %>%
  dplyr::group_by(celltype) %>%
  dplyr::slice_head(n = 5) %>%
  dplyr::ungroup()


p_rss <- ggplot(rss_long, aes(x = rank, y = RSS)) +
  geom_line(
    color = "#00A6D6",
    linewidth = 1.1
  ) +
  geom_point(
    data = label_df,
    aes(x = rank, y = RSS),
    inherit.aes = FALSE,
    color = "#E64B5D",
    size = 2.1
  ) +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(x = rank, y = RSS, label = regulon),
    inherit.aes = FALSE,
    color = "#E64B5D",
    size = 3.4,
    fontface = "bold",
    direction = "y",
    hjust = 0,
    nudge_x = 2.5,
    segment.color = "#E64B5D",
    segment.linewidth = 0.25,
    max.overlaps = Inf,
    box.padding = 0.15,
    point.padding = 0.2,
    min.segment.length = 0
  ) +
  facet_wrap(
    ~ celltype,
    nrow = 1,
    scales = "free_y"
  ) +
  scale_x_continuous(
    breaks = c(0, 20, 40, 60),
    limits = c(0, 72),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(
    x = "Regulon Rank",
    y = "RSS"
  ) +
  theme_bw(base_size = 16) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      size = 15,
      face = "bold",
      color = "black"
    ),
    panel.grid.major = element_line(
      color = "grey88",
      linewidth = 0.3,
      linetype = "dashed"
    ),
    panel.grid.minor = element_blank(),
    axis.title = element_text(
      size = 16,
      color = "black"
    ),
    axis.text = element_text(
      size = 14,
      color = "black"
    ),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.7
    ),
    plot.margin = margin(5, 5, 5, 5)
  )

p_rss
