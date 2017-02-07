## HGSS R Workshop 2

## data.frame basics
myDF = data.frame(gene=c("A","B","C"),id=1:3)
myDF
myDF$gene
str(myDF)
myDF2 = data.frame(gene=c("A","B","C"),id=1:3, stringsAsFactors = FALSE)
myDF2
myDF2$gene
str(myDF2)
colnames(myDF2)[2] = "ID"
myDF2

## Melting matrices
load("dataWS2.RData")
str(mat.ge)
mat.ge[1:10,1:10]
mat.ge.10 = mat.ge[1:10,]
dim(mat.ge.10)
library(reshape)
mat.ge.10.df = melt(mat.ge.10)
head(mat.ge.10.df)
colnames(mat.ge.10.df) = c("gene","sample","expression")
head(mat.ge.10.df)
str(mat.ge.10.df)
mat.ge.10.df$gene = as.character(mat.ge.10.df$gene)
mat.ge.10.df$sample = as.character(mat.ge.10.df$sample)
str(mat.ge.10.df)

## Merging df
df1 = data.frame(colA=c(1,1,2), colB=c("a","b","c"))
df2 = data.frame(colA=c(2,1,3), colC=c(TRUE,FALSE,FALSE))
df1
df2
merge(df1,df2)
merge(df1,df2, all=TRUE)

head(metadata.df)
mat.ge.10.df = merge(mat.ge.10.df,metadata.df)
head(mat.ge.10.df)

## Subset
head(subset(mat.ge.10.df, gender=="male"))
subset(mat.ge.10.df, gene=="gene1")

## Histogram
library(ggplot2)

ggplot(mat.ge.10.df, aes(x=expression)) + geom_histogram()
ggplot(mat.ge.10.df, aes(x=expression)) + geom_histogram(fill="blue")
ggplot(mat.ge.10.df, aes(x=expression, fill=gene)) + geom_histogram()
ggplot(mat.ge.10.df, aes(x=expression, fill=gene)) + geom_histogram(position="dodge")
ggplot(mat.ge.10.df, aes(x=expression, fill=gender)) + geom_histogram(position="dodge")
ggplot(mat.ge.10.df, aes(x=expression, fill=gender)) + geom_histogram()
ggplot(mat.ge.10.df, aes(x=expression, fill=gender)) + geom_histogram(position="fill")

ggplot(mat.ge.10.df, aes(x=expression, fill=gene)) + geom_histogram() +
  xlab("gene expression") + ggtitle("Expression of 10 genes") + xlim(2,6) + 
  theme_bw() + theme(legend.position="top")

ggplot(mat.ge.10.df, aes(x=expression, fill=gene)) + geom_histogram() +
  xlab("gene expression") + ggtitle("Expression of 10 genes") + xlim(2,6) + 
  theme_bw() + theme(legend.position=c(.8,.8))

## Facets
## Wrap NO DOT
ggplot(mat.ge.10.df, aes(x=expression)) + geom_histogram(binwidth=.2) +
  theme_bw() + facet_wrap(~gene, ncol=3, scales="free")

ggplot(mat.ge.10.df, aes(x=expression)) + geom_histogram(binwidth=.2) +
  theme_bw() + facet_grid(gender~gene, scales="free")

ggplot(mat.ge.10.df, aes(x=expression)) + geom_histogram(binwidth=.2) +
  theme_bw() + facet_grid(gene~.)

ggplot(mat.ge.10.df, aes(x=expression)) + geom_histogram(binwidth=.2) +
  theme_bw() + facet_grid(gene~., scales="free")

## Reading big file
gc = read.table("gencode.v22.annotation.gtf",sep="\t",as.is=TRUE, nrows=10)
head(gc)
system.time({gc = read.table("gencode.v22.annotation.gtf",sep="\t",as.is=TRUE)}) ## Took me 167 s

colC = c("character", "NULL","character", "numeric","numeric","NULL","NULL","NULL","character")
system.time({gc = read.table("gencode.v22.annotation.gtf",sep="\t",colClasses=colC,nrows=3e6)})  ## Took me 142 s
head(gc)

library(data.table)
system.time({gc.DT = fread("gencode.v22.annotation.gtf", skip=5, select=c(1,3,4,5,9))}) ## Took me 12 s
gc = as.data.frame(gc.DT)
colnames(gc) = c("seqname","feature","start","end","attribute")
gc[1:5,1:4]

head(gc$feature)
unique(gc$feature)
table(gc$feature)

ggplot(gc, aes(x=feature)) + geom_bar()
ggplot(gc, aes(x=seqname)) + geom_bar()
ggplot(gc, aes(x=seqname, fill=feature)) + geom_bar()
ggplot(gc, aes(x=seqname, fill=feature)) + geom_bar(position="dodge") + theme_bw()

## Order X axis
## Use factors -> order specified in the data
chrs = unique(gc$seqname)
gc$chr.factor = factor(gc$seqname, levels=chrs)
ggplot(gc, aes(x=chr.factor, fill=feature)) + geom_bar()
## OR change it for the graph only
ggplot(gc, aes(x=seqname, fill=feature)) + geom_bar() + 
  scale_x_discrete(limits=chrs)

## Get attributes
getAtt <- function(attributes, att.name="gene_id"){
  sub(paste0(".*",att.name," \"([^\"]+)\";.*"), "\\1", attributes)
}
## For read.table users
getAtt <- function(attributes, att.name="gene_id"){
  sub(paste0(".*",att.name,"  ([^;]+);.*"), "\\1", attributes)
}
##
getAtt(head(gc$attribute))
getAtt(head(gc$attribute), att.name = "gene_type")

library(dplyr)
gc$attribute %>% head %>% getAtt ## Pipes for fun, not important

gc = gc %>% mutate(geneId=getAtt(attribute),
                   geneType=getAtt(attribute, "gene_type"))
str(gc)
save(gc, file="backupWS2.RData")
## load("backupWS2.RData")

ggplot(subset(gc, feature=="gene"), aes(x=geneType)) + geom_bar() + 
  coord_flip()
## Many different gene types with almost no genes

## Top 10 longest pseudogenes
gc %>% filter(feature=="gene" & geneType=="pseudogene") %>% 
  mutate(length=end-start+1) %>% arrange(desc(length)) %>%
    select(geneType,length, geneId) %>% head(10)

## Top 10 shortest protein_coding exons
gc %>% filter(feature=="exon" & geneType=="protein_coding") %>% 
  mutate(length=end-start+1) %>% arrange(length) %>%
    select(geneType,length, geneId) %>% head(10)

## Gene size distribution
ggplot(subset(gc, feature=="gene"), aes(x=end-start+1, fill=geneType)) + 
  geom_histogram() + scale_x_log10() + xlab("gene size")

ggplot(subset(gc, feature=="gene"), aes(x=end-start+1, fill=geneType)) + 
  geom_density() + scale_x_log10()  + xlab("gene size") ## Beautiful

ggplot(subset(gc, feature=="gene"), aes(x=end-start+1, colour=geneType)) + 
  geom_density() + scale_x_log10()  + xlab("gene size")



##
## Grouping rows
##

## Number of genes per gene type
nbg.df = gc %>% filter(feature=="gene") %>% group_by(geneType) %>%
  summarize(nbGenes=n()) %>% arrange(desc(nbGenes))

nbg.df
head(nbg.df$geneType)
gc$geneType[which(!(gc$geneType %in% head(nbg.df$geneType)))] = "others"

ggplot(subset(gc, feature=="gene"), aes(x=end-start, fill=geneType)) + 
  geom_histogram() + scale_x_log10() + scale_fill_brewer(palette="Set2")
  
## Number of exons per gene
nbex.df = gc %>% filter(feature=="exon") %>% group_by(geneId, geneType) %>%
  summarize(n=n())

ggplot(nbex.df, aes(x=n)) + geom_histogram() + theme_bw() + xlab("number of exons")

ggplot(nbex.df, aes(x=n, fill=geneType)) + geom_histogram() + theme_bw() +
  xlab("number of exons") + xlim(0,100)

nbex.df %>% arrange(desc(n)) %>% head  ## Doesn't work because still grouped
nbex.df %>% ungroup %>% arrange(desc(n)) %>% head  ## works

## Bonus : Average exon size per gene per gene type
exon.size = gc %>% filter(feature=="exon") %>% group_by(geneId,geneType) %>%
  summarize(exon.size=mean(end-start))

ggplot(exon.size, aes(x=exon.size, fill=geneType)) + geom_histogram() + theme_bw() +
  scale_x_log10() + theme(legend.position=c(1,1),legend.justification=c(1,1))

exon.size %>% group_by(geneType) %>% summarize(exon.size=mean(exon.size)) %>% 
  arrange(desc(exon.size)) %>% head ## Average per gene type

##
## Scatterplots
##

## Number of transcript per gene
nbtr.df = gc %>% filter(feature=="transcript") %>% group_by(geneId, geneType) %>%
  summarize(nbtr=n())

ggplot(nbtr.df, aes(x=nbtr, fill=geneType)) + geom_histogram()

## Size of genes
gsize = gc %>% filter(feature=="gene") %>% mutate(gsize=end-start+1)
head(gsize)

## Merge and scatterplot
nbtr.df = merge(nbtr.df, gsize)

ggplot(nbtr.df, aes(x=nbtr, y=gsize)) + geom_point()

ggplot(nbtr.df, aes(x=nbtr, y=gsize, colour=geneType, shape=geneType)) +
  geom_point(alpha=.7) + scale_y_log10() + theme_bw() +
    theme(legend.position="bottom") 
## Actually too many shapes asked for...

ggplot(nbtr.df, aes(x=nbtr, y=gsize, colour=geneType)) +
  geom_point(alpha=.7) + scale_y_log10() + theme_bw() +
    theme(legend.position="bottom") 

ggplot(nbtr.df, aes(x=nbtr, y=gsize)) + facet_wrap(~geneType, scales="free") + 
  geom_point(alpha=.5) + scale_y_log10() + theme_bw() + xlab("number of transcripts") +
    ylab("gene size (bp)") + theme(legend.position="bottom") 


##
## Genomic Ranges
##
library(GenomicRanges)

## Exon density : practicing on one gene and trying GR functions
gene1 = subset(gc, geneId=="ENSG00000223972.5" & feature=="exon")
dim(gene1)
gr1 = with(gene1, GRanges(seqname, IRanges(start, end)))
gr1
start(gr1)
end(gr1)
width(gr1)
range(gr1)
gr1.red = reduce(gr1)
gr1.red

ex.dens1 = sum(width(gr1)) / width(range(gr1))
ex.dens1 ## 90% of the gene is covered by exons

## Exon density : function and apply to blocks
exonDens <- function(inDF){
  gr = with(inDF, GRanges(seqname, IRanges(start, end)))
  gr.red = reduce(gr)
  ex.dens = sum(width(gr.red)) / width(range(gr))
  return(data.frame(exon.dens=ex.dens)) ## RETURN A DATA.FRAME
}

gc %>% filter(feature=="exon" & seqname=="chr13") %>% group_by(geneId) %>%
  do(exonDens(.))

#### This would run the previous column sequentially on four chromosomes
## ex.dens.list = lapply(c("chr13","chr14","chr15","chr16"), function(seqn){
##   res.df = gc %>% filter(feature=="exon" & seqname==seqn) %>% group_by(geneId) %>% do(exonDens(.))
##   return(res.df)
## })
####

## Changing it to run on 2 processors in parallel would be :
library(parallel)
ex.dens.list = mclapply(c("chr13","chr14","chr15","chr16"), function(seqn){
  res.df = gc %>% filter(feature=="exon" & seqname==seqn) %>% group_by(geneId) %>% do(exonDens(.))
  return(res.df)
}, mc.cores=2)
## The result is a list of data.frames, let's bind everything together
ex.dens.df = do.call(rbind, ex.dens.list)
head(ex.dens.df)

ggplot(ex.dens.df, aes(x=exon.dens)) + geom_histogram() + theme_bw()


##
## AnnotationHub
##
library(AnnotationHub)
ah = AnnotationHub()
hist.prom = ah$goldenpath.hg19.encodeDCC.wgEncodeBroadHistone.wgEncodeBroadHistoneGm12878H3k4me3StdPk.broadPeak_0.0.1.RData
head(hist.prom) ## GRanges with histone marks
mean(width(hist.prom))

## Quick ggplot when you don't have a data.frame. Adding a dotted line at 10kbp.
qplot(x=width(hist.prom)) + geom_histogram() + theme_bw() + scale_x_log10() +
  geom_vline(xintercept=1e4, linetype=2)

## Remove peaks larger than 10Kbp
hist.prom = hist.prom[which(width(hist.prom)<1e4)]

##
## Test overlaps
##
gc.gr = with(gc, GRanges(seqname, IRanges(start, end)))
gc$overlap.prom = overlapsAny(gc.gr, hist.prom)  ## Works because gc and gc.gr are in the same order; gc.gr is just the GRanges version

## Number of genes in each gene types overlapping or not histone peaks
ggplot(subset(gc, feature=="gene"), aes(x=geneType,fill=overlap.prom)) +
  geom_bar(position="dodge") + theme_bw() + coord_flip()
## -> Most genes don't overlap peaks.
## -> Some non-protein coding genes do overlap histone marks.

## Proportion of each feature overlapping histone peaks
ggplot(gc, aes(x=feature,fill=overlap.prom)) + geom_bar(position="fill") +
  theme_bw() 

## Same information but precomputing the information with dplyr (faster)
gc %>% group_by(feature) %>% summarize(ol.prop=mean(overlap.prom)) %>%
  ggplot(aes(x=feature, y=ol.prop)) + geom_bar(stat="identity") +
    theme_bw() + ylim(0,1) + ylab("proportion overlapping histone marks")

##
## Find overlaps
##
ol = findOverlaps(gc.gr, hist.prom)
## Overlapping 2.5 millions regions with 50000 regions took a few seconds !
ol
ol.1kb = findOverlaps(gc.gr, hist.prom, maxgap=1000) ## +/- 1 kbp
ol.1kb

## Create a new data.frame to merge histone peak score
gc.hscore = gc[queryHits(ol.1kb),]
gc.hscore$hscore = hist.prom$score[subjectHits(ol.1kb)]
str(gc.hscore) ## Checking

ggplot(gc.hscore, aes(x=hscore, fill=geneType)) + geom_histogram() + theme_bw()

ggplot(gc.hscore, aes(x=hscore)) + geom_histogram() + theme_bw() +
  facet_wrap(~geneType, scales="free")

##
## Distance to nearest
##
gc.gene = subset(gc, feature=="gene")
gc.gene.gc =  with(gc.gene, GRanges(seqname, IRanges(start, end)))
dhprom = distanceToNearest(gc.gene.gc, hist.prom)
gc.gene$dist.hprom = mcols(dhprom)$distance
str(gc.gene) ## Checking

ggplot(gc.gene, aes(x=dist.hprom, colour=geneType)) + geom_density() +
  theme_bw() + scale_x_log10()
## Caution log scale remove the 0s because log(0) is not defined

ggplot(gc.gene, aes(x=dist.hprom+1, colour=geneType)) + geom_density() +
  theme_bw() + scale_x_log10()
## Indeed here are all the genes that overlaped histhon marks (dist.hprom=0)



##
## Boxplots
##

## Remember the gene expression data.frame ?
head(mat.ge.10.df)

ggplot(mat.ge.10.df, aes(x=gender, y=expression, fill=group)) + geom_boxplot() +
  theme_bw()

ggplot(mat.ge.10.df, aes(x=gene, y=expression, fill=group)) + geom_boxplot() +
  theme_bw() + facet_grid(gender~.)
