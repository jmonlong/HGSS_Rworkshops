# Advanced R: Tidyverse and Bioconductor

*March 26 2018*

## Install packages for the workshop

Run the following commands:

```r
## Install CRAN packages
install.packages(c('data.table','dplyr','ggplot2', 'parallel', 'magrittr', 'tidyr'))

## Install Bioconductor packages
source('http://bioconductor.org/biocLite.R')
biocLite(c('GenomicRanges', 'VariantAnnotation', 'clusterProfiler', 'EnrichedHeatmap', 'AnnotationHub', 'Gviz'))
```
