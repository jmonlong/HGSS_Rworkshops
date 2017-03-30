## Live script of the workshop

## fread
genc.df = read.table("gencode.gtf.gz", sep="\t", as.is=TRUE)

library(data.table)
genc.dt = fread("gunzip -c gencode.gtf.gz")
genc.df = as.data.frame(genc.dt)
unique(genc.df[,1])
sum(genc.df[,3] == "gene")

con = file("gencode.gtf.gz", "r")
read.table(con, nrows = 1)
read.table(con, nrows = 1)
close(con)

con = file("gencode.gtf.gz", "r")
while(length((chunk.df = read.table(con,nrows=1000)))>0){
  .. INSTRUCTIONS
}
close(con)

## Trick: assign new object AND use a function
length((t = 1:10))

## Use Bioconductor packages
library(rtracklayer)
genc = import("gencode.gtf.gz")
head(genc)
unique(genc$gene_type)
length(unique(genc$gene_type))
length(unique(genc$gene_name))

for(gt in unique(genc$gene_type)){
  genc.gt = subset(genc, gene_type==gt)
  message(gt, " ", length(unique(genc.gt$gene_name)))
}

countUnique = function(x) length(unique(x))
genc.df = as.data.frame(genc)
aggregate(gene_name ~ gene_type, data=genc.df, countUnique)

tapply(genc$gene_name, genc$gene_type, countUnique)

## Indexed files
library(GenomicRanges)
reg30 = GRanges(22, IRanges(30e6, 31e6))

library(VariantAnnotation)
tgp = readVcf("tgp22.vcf.gz", "hg19", reg30)
## or
tgp = readVcf(TabixFile("tgp22.vcf.gz"), "hg19", reg30)
head(tgp)
length(tgp)
rowRanges(tgp)

library(data.table)
genc.dt = fread("gunzip -c gencode.gtf.gz")
setkey(genc.dt, V1, V4)
write.table(genc.dt, file="gencode-ordered.gtf.gz", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

library(Rsamtools)
bgzip("gencode-ordered.gtf.gz")
indexTabix("gencode-ordered.gtf.bgz", format="gff")

## Avoid loops
## for loop BAD
nb.exon = rep(NA, 1000)
for(ii in 1:1000){
  type10 = sample(genc$type, 10)
  nb.exon.temp = sum(type10 == "exon")
  nb.exon[ii] = nb.exon.temp
}

## Sampling 10 elements
sample(genc$type, 10)

sample10count <- function(ii){
  type10 = sample(genc$type, 10)
  nb.exon = sum(type10 == "exon")
  data.frame(perm=ii, nb.exon=nb.exon)
}
perm.l = lapply(1:1000, sample10count)

perm.df = do.call(rbind , perm.l)
## or
perm.df = as.data.frame(rbindlist(perm.l))

head(perm.df)

countGenesPerGT <- function(gt){
  genc.gt = subset(genc, gene_type==gt)
  data.frame(gene_type=gt, genes=length(unique(genc.gt$gene_name)))
}
lapply(unique(genc$gene_type), countGenesPerGT)

## Parallel
library(parallel)
perm.l = mclapply(1:500, sample10count, mc.cores=3)


library(BatchJobs)

sample10countBatchJobs <- function(ii, genc){
  type10 = sample(genc$type, 10)
  nb.exon = sum(type10 == "exon")
  data.frame(perm=ii, nb.exon=nb.exon)
}

reg = makeRegistry("perm")

batchMap(reg, sample10countBatchJobs, 1:10, more.args=list(genc=genc))

submitJobs(reg, 1:4)
showStatus(reg)
perm.l = reduceResultsList(reg)
do.call(rbind, perm.l)


## dplyr
library(dplyr)
genc.df = as.data.frame(genc)

genc.nbexon = genc.df %>% group_by(gene_name) %>%
  summarize(nb.exon=sum(type=="exon"))
head(genc.nbexon)

genc.df %>% filter(type=="exon") %>% group_by(gene_name) %>% summarize(nb.exon=n())

genc.df %>% group_by(gene_name) %>% summarize(nb.exon=sum(type=="exon")) %>% arrange(desc(nb.exon))

genc.df %>% group_by(gene_name) %>% summarize(nb.exon=sum(type=="exon")) %>% arrange(nb.exon)

genc.df %>% group_by(gene_name, gene_type) %>% summarize(nb.exon=sum(type=="exon")) %>% group_by(gene_type) %>% summarize(nb.exon.ave=mean(nb.exon)) %>% arrange(desc(nb.exon.ave)) %>% head

genc.df %>% filter(type=="gene") %>% group_by(gene_type) %>% summarize(size.mean=mean(width), size.median=median(width)) %>% arrange(desc(size.mean))

genc.df %>% group_by(seqnames, gene_type) %>% summarize(nb.gene=sum(type=="gene"))

genc.df %>% filter(gene_type=='protein_coding', type=="exon") %>% mutate(exon_number=as.numeric(exon_number)) %>% group_by(exon_number) %>% summarize(mean=mean(width), med=median(width))

testSize <- function(df){
  lmo = lm(width~as.numeric(exon_number), data=df)
  lmos = summary(lmo)
  data.frame(pv=lmos$coefficients[2,4], coef=lmos$coefficients[2,1])
}
genc.df %>% filter(type=="exon") %>% group_by(gene_name) %>% filter(length(unique(exon_number))>1) %>% do(testSize(.))


tgp ## Coordinates + info + genotype
rowRanges(tgp) ## the coordinates of the variants
genc

## Change chr names
seqlevels(genc) = gsub("chr","",seqlevels(genc))


## Overlap
genc.tgp = subsetByOverlaps(genc, rowRanges(tgp))
genc.tgp %>% as.data.frame %>% group_by(gene_type) %>% summarize(gene=length(unique(gene_name)))

res.l = lapply(unique(genc$gene_type), function(gt){
  genc.gt = subset(genc, gene_type==gt)
  genc.tgp = subsetByOverlaps(genc.gt, rowRanges(tgp))
  data.frame(gene_type=gt, gene=length(unique(genc.tgp$gene_name)))
})
do.call(rbind, res.l)

ol = findOverlaps(genc, rowRanges(tgp))
ol %>% as.data.frame %>% mutate(gene_name=genc$gene_name[queryHits], gene_type=genc$gene_type[queryHits]) %>% group_by(gene_type) %>% summarize(gene=length(unique(gene_name)))

ol %>% as.data.frame %>% mutate(gene_name=genc$gene_name[queryHits], gene_type=genc$gene_type[queryHits]) %>% group_by(gene_name) %>% summarize(variant=length(unique(subjectHits)))

ol %>% as.data.frame %>% mutate(gene_name=genc$gene_name[queryHits], gene_type=genc$gene_type[queryHits]) %>% group_by(gene_name, gene_type) %>% summarize(variant=length(unique(subjectHits))) %>% group_by(gene_type) %>% summarize(var.mean=mean(variant))

info(tgp)
tgp.gr = rowRanges(tgp)
tgp.gr$af = sapply(info(tgp)$AF, "[", 1)


findOverlaps(genc, tgp.gr) %>% as.data.frame %>%
  mutate(gene_type=genc$gene_type[queryHits],
         af=tgp.gr$af[subjectHits]) %>% group_by(gene_type) %>% summarize(af=mean(af))

## ggplot2
library(ggplot2)
genc.df %>% filter(type=="gene", gene_type %in% c("protein_coding", "lincRNA")) %>% ggplot(aes(x=width, fill=gene_type)) + geom_histogram() + scale_x_log10() + theme_bw()


genc.df %>% filter(type=="gene", !grepl("IG", gene_type), !grepl("TR", gene_type), !grepl("RNA", gene_type)) %>% ggplot(aes(x=width, fill=gene_type)) + geom_histogram() + scale_x_log10() + theme_bw() + scale_fill_brewer(palette="Set1") + theme(legend.position="bottom")

genc.df %>% filter(gene_type=='protein_coding', type=="exon") %>% mutate(exon_number=as.numeric(exon_number)) %>% group_by(exon_number) %>% summarize(mean=mean(width), med=median(width)) %>% filter(exon_number<=30) %>% ggplot(aes(x=exon_number, y=mean)) + geom_point() + geom_line() + theme_bw() + xlab("exon number") + ylab("average size (bp)")


findOverlaps(genc, tgp.gr) %>% as.data.frame %>% mutate(gene_type=genc$gene_type[queryHits], af=tgp.gr$af[subjectHits]) %>% ggplot(aes(x=reorder(gene_type, af, median), y=af)) + geom_boxplot() + scale_y_log10() + coord_flip() + theme_bw() + xlab("gene type") + ylab("allele frequency")






