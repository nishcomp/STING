# "Segmenting Transcriptomics in Spatial Neighborhoods"

## Objectives
1. To find the spatially variable genes on the tissue based on spatial autocorrelation statistics (Moran's I)
2. To take the spatially variable genes and study their hotspots and immediate neighborhood

## Package Dependencies
- Seurat
- spdep
- sf
- doParallel
- KNN
- dplyr
- ggplot2

## Basic Workflow

```r
library(ggraph)
library(igraph)
library(tidyverse)
library(RColorBrewer)
library(dplyr)
library(STING)
```

```r
seurat_obj <- Load10X_Spatial(data.dir = sample_path)
seurat_obj <- NormalizeData(seurat_obj)
```

## caclulates the spatial autocorrelation of a gene based on its neighborhood (k, number of neighbors)

**The sample studied here is P1CRC, downloaded from [10XVisium](https://www.10xgenomics.com/products/visium-hd-spatial-gene-expression/dataset-human-crc).**

```r
svgs_moran_test <- calculate_moranSVGs(sample_name = "colon10x", seurat_obj=seurat_obj, layer = "data", k = 20,
                   output_path = output_dir, min_count = 100)
                   
# outpus a csv file with Moran's I statistics, p-value, and fdr-corrected p-value for each gene present in > min_count #cells.
head(svgs_moran_test)
```
![SVG_calculated](images/calculate_SVG_output.png)

## group the genes in different categories 
*for the classification, i used protein atlas, cytokines and chemokines, immunoglobulins etc list obtained from the internet.*
Or simply label "others" or "uncategorized" to make it compatible with the plotting function

![Categories Added](images/add_gene_category.png)


## plot the obtained genes and their statistics

1. Circular Plot

```r
plot_svgs(
    df = csvgs_moran_test], 
    title = "Your Title", 
    type = "circular", 
    genes = "gene", 
    range = c(1, 10), 
    output_file = paste0("./outputdir/", "sample", "circular_SVGS.png")
  )
```

![Cicular Plot](images/p1crccircular_SVGS.png)

2. Linear Plot

```r
plot_svgs(
    df = csvgs_moran_test], 
    title = "Your Title", 
    type = "linear", 
    genes = "gene", 
    range = c(1, 10), 
    output_file = paste0("./outputdir/", "sample", "linear_SVGS.png")
  )

```

![Cicular Plot](images/p1crclinear_SVGS.png)

