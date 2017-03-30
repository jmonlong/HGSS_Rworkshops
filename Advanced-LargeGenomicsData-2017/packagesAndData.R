install.packages(c("data.table","dplyr","ggplot2","parallel", "BatchJobs", "parallel", "magrittr"))
source("http://bioconductor.org/biocLite.R")
biocLite(c("GenomicRanges","Rsamtools", "VariantAnnotation", "rtracklayer", "Gviz"))

download.file("ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz", "gencode.gtf.gz")

download.file("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz", "tgp22.vcf.gz")
download.file("ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi", "tgp22.vcf.gz.tbi")

