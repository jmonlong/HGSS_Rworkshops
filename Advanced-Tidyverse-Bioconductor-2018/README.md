# Advanced R: Tidyverse and Bioconductor

*March 26 2018*

## Install packages for the workshop

Run the following commands:

```r
## Install CRAN packages
install.packages(c('data.table','dplyr','ggplot2', 'magrittr', 'tidyr'))

## Install Bioconductor packages
source('http://bioconductor.org/biocLite.R')
biocLite(c('GenomicRanges', 'VariantAnnotation', 'clusterProfiler', 'EnrichedHeatmap', 'AnnotationHub', 'Gviz', 'org.Hs.eg.db'))
```

## Data

We'll use a simplified version of the Gencode annotation. 
Download it [here](https://github.com/jmonlong/HGSS_Rworkshops/raw/master/Advanced-Tidyverse-Bioconductor-2018/gencodeForWorkshop.tsv.gz) or by running the following command in R:

```r
download.file('https://github.com/jmonlong/HGSS_Rworkshops/raw/master/Advanced-Tidyverse-Bioconductor-2018/gencodeForWorkshop.tsv.gz','gencodeForWorkshop.tsv.gz')
```

## Slides

- `HGSS-Rworkshop2018-advanced.pdf` contains the slides.
- `HGSS-Rworkshop2018-advanced.tex` contains the LaTeX source.
