parser = argparse::ArgumentParser()
parser$add_argument('-i', '--input', dest = 'input', help = 'input rds filename')
parser$add_argument('-r', '--res', dest = 'res', default = 'lowres', help = 'hires or lowres')
parser$add_argument('-o', '--out', dest = 'outdir', help = 'directory where to save the output files')
opts = parser$parse_args()

library(Seurat)
library(Matrix)

spatial_seurat <- readRDS(opts$input)

# 提取 count 矩阵并导出（通常空间转录组推荐用 Spatial 或 SCT assay）
counts <- GetAssayData(spatial_seurat, assay = "Spatial", slot = "counts")
writeMM(counts, file = paste0(opts$out, "/counts.mtx"))

# 提取基因名和细胞名
write.table(rownames(counts), file = paste0(opts$out, "/features.tsv"), col.names = FALSE, row.names = FALSE, quote = FALSE)
write.table(colnames(counts), file = paste0(opts$out, "/barcodes.tsv"), col.names = FALSE, row.names = FALSE, quote = FALSE)

# 提取元数据 (Metadata)
write.csv(spatial_seurat@meta.data, file = paste0(opts$out, "/metadata.csv"), row.names = TRUE)

# 提取空间坐标 (Crucial for Spatial!)
# 获取空间图像对象的名称（例如 "slice1"）
image_name <- names(spatial_seurat@images)[1] 
#spatial_coords <- GetTissueCoordinates(spatial_seurat, image = image_name)
#write.csv(spatial_coords, file = "spatial_coordinates.csv", row.names = TRUE)

slice <- spatial_seurat@images[[image_name]]
# 【关键】获取和低分辨率图片对齐的坐标
# 不要直接用 GetTissueCoordinates，有些版本的 Seurat 导出的是 fullres 坐标
# 我们手动乘以 lowres 的缩放因子
coords <- slice@coordinates
fac <- slice@scale.factors[[opts$res]]
spatial_coords <- data.frame(
  imagerow = coords$imagerow * fac,
  imagecol = coords$imagecol * fac
)
rownames(spatial_coords) <- rownames(coords)
write.csv(spatial_coords, file = paste0(opts$out, "/spatial_coordinates_", opts$res, ".csv"))
