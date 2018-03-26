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

## Live script

The commands I typed during the workshop:

- `HGSS-Rworkshop2018-advanced-liveScript.Rmd`: the R Markdown source file.
- `HGSS-Rworkshop2018-advanced-liveScript.html`: the HTML output.
- [`HGSS-Rworkshop2018-advanced-liveScript.md`](HGSS-Rworkshop2018-advanced-liveScript.html): the Markdown output for GitHub.

## To go further

### Bioconductor resources

- [Workflows](http://bioconductor.org/help/workflows/) gives detailed hands-on tutorials for different bioinformatics analysis, such as gene expression, epigenetics, proteomics analysis and others.
- [Courses](http://bioconductor.org/help/course-materials/): slides from courses.
- [Videos](https://www.youtube.com/user/bioconductor/videos).
