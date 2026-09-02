library(Seurat)
library(png)
library(ggplot2)
library(jsonlite)

parser = argparse::ArgumentParser()
parser$add_argument('-r', '--rds', dest = 'rds', help = 'input rds filename')
parser$add_argument('--im_dir', dest = 'image', help = 'image dir')
parser$add_argument('-o', '--out', dest = 'outdir', help = 'directory where to save the output files')
opts = parser$parse_args()

obj <- readRDS(opts$rds)
spatial_dir <- opts$image

### read hires image
img <- png::readPNG(file.path(spatial_dir, "tissue_hires_image.png"))

### read coords
pos <- read.csv(
    file.path(spatial_dir, "tissue_positions_list.csv"),
    header = FALSE,
    colClasses = "character"
)
colnames(pos) <- c(
    "barcode",
    "tissue",
    "y_fullres",
    "x_fullres",
    "y2_fullres",
    "x2_fullres"
)

sf <- fromJSON(file.path(spatial_dir, "scalefactors_json.json"))
hires_scale <- sf$tissue_hires_scalef
print(hires_scale)

pos$x <- as.numeric(pos$x_fullres) * hires_scale
pos$y <- as.numeric(pos$y_fullres) * hires_scale

pos <- pos[match(colnames(obj), pos$barcode),]

stopifnot(all(pos$barcode == colnames(obj)))

hires_img <- Read10X_Image(
    image.dir = spatial_dir,
    image.name = "tissue_hires_image.png",
    assay = "Spatial",
    slice = "slice1",
    filter.matrix = FALSE
)

hires_img@coordinates <- data.frame(
    tissue = as.numeric(pos$tissue),
    imagerow = as.numeric(pos$y_fullres),
    imagecol = as.numeric(pos$x_fullres),
    row.names = pos$barcode
)

hires_img@scale.factors$lowres <- hires_img@scale.factors$hires
DefaultAssay(hires_img) <- "Spatial"
obj[["slice1"]] <- hires_img

SpatialFeaturePlot(obj, features='nCount_Spatial', stroke=NA, pt.size.factor=1) + theme(aspect.ratio = NULL)
ggsave(paste0(opts$out, '/test_plot.png'))

saveRDS(obj, paste0(opts$out, '/output_hires.rds'))
