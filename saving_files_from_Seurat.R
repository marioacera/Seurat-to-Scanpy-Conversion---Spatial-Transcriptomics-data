library('zellkonverter')
library('Seurat')
library(patchwork)
library(SingleCellExperiment)
library(png)
library('animation')
library(EBImage)
library(threadr)

data <- readRDS('<set_path>/seurat_list.rds')

data1 <- data$CTRL_1
SpatialFeaturePlot(data1, features = "nCount_Spatial" )
write.csv(data1@images$CTRL_1@coordinates, '<set_path>/<sample_name>/spatial/tissue_position_list.csv')

toJSON <- as.list(data1@images$CTRL_1@scale.factors)
class(toJSON)<-"list"
write_json(toJSON, file = '<set_path>/<sample_name>/spatial/scalefactors_json.json')

## Matrices corresponding to red, green and blue color channels
r <- data1@images$CTRL_1@image[,,1]
g <- data1@images$CTRL_1@image[,,2]
b <- data1@images$CTRL_1@image[,,3]

## Construct an color Image object . If not R will save 3 images with one channel each
img <- rgbImage(r, g, b)
display(a)
writeImage(img,'<set_path>/<sample_name>/spatial/lowres_tissue_img_basic.png')

sce<- as.SingleCellExperiment(data1)

writeH5AD(sce, file = '<set_path>/<sample_name>/filtered_featured_bc_matrix.h5ad')


library(Seurat)
library(ggplot2)
library(iotools)
library(SeuratDisk)


Convert(
  source="<set_path>/<sample_name>/sample_results.h5ad",
  dest = "h5seurat",
  assay = "RNA",
  overwrite = FALSE,
  verbose = TRUE
)

data <- LoadH5Seurat("sample_results.h5ad.h5seurat")

