## Vector construction
luckyNumbers = c(4,8,15,16,23,42)
luckyNumbers
oneToTen = 1:10
tenOnes = rep(1,10)
samples = c("sampA","sampB","sampC")
samples

## Vector manipulation
## Using indexes
luckyNumbers[3]
luckyNumbers[3] = 79
luckyNumbers[2:4]
luckyNumbers[2:4] = c(14,3,9)
luckyNumbers[c(1,6)]
luckyNumbers[-2]
## Using names
names(luckyNumbers)
names(luckyNumbers) = c("frank","henry","philip","steve","tom","francis")
names(luckyNumbers) = NULL
names(luckyNumbers) = c("frank","henry","philip","steve","tom","francis")
luckyNumbers["philip"]
luckyNumbers["philip"] = 79
luckyNumbers[c("henry","philip","steve")]
luckyNumbers[c("henry","philip","steve")] = c(14,3,9)
luckyNumbers[c("frank","francis")]
names(luckyNumbers)[5]
names(luckyNumbers)[5] = "tomtom"
names(luckyNumbers)[1:3]
names(luckyNumbers)[1:3] = c("newName1","newName2","newName3")
luckyNumbers = c(luckyNumbers, 65)
names(luckyNumbers)[7] = "simon"
## Size of the vector
length(luckyNumbers)
luckyNumbers[length(luckyNumbers)]

## More manipulation
sort(luckyNumbers)
luckyNumbers[order(luckyNumbers)]
sort(luckyNumbers,decreasing=TRUE)
c(luckyNumbers,1:10,tenOnes)
sort(c(luckyNumbers,1:10,tenOnes))
rev(1:10)
sample(1:10)
sample(1:10,2)

## Exploration
head(luckyNumbers)
tail(luckyNumbers)
summary(luckyNumbers)
min(luckyNumbers)
max(luckyNumbers)
mean(luckyNumbers)
mean(luckyNumbers[1:3])
mean(luckyNumbers[c(1,6)])

## Operations
luckyNumbers * 4
luckyNumbers - luckyNumbers
luckyNumbers * 1:length(luckyNumbers)

## Get my favorite number
n = c(3,12,458,-449826.67)
((n*6+21)/3-1)/2 -n

## Matrix
neo = matrix(1:12,3,4)
neo
neo[2,4]
neo[2,4] = 0
neo[,1]
neo[,1] = rep(21,3)
neo[3,]
sample(neo[3,])
neo[3,] = sample(neo[3,])
neo[1:2,1:3]
neo[1:2,1:3] = matrix(rev(1:6),2,3)
neo[c(1,3),c(2,4)]

## Manipulating more
onesMat = matrix(rep(1,6),2,3)
dim(onesMat)
onesMat
2 * onesMat
rbind(onesMat, 2*onesMat)
cbind(onesMat, 2*onesMat)
cbind(neo, rep(100,3))
neo = rbind(neo, rep(100,4))
cbind(neo[,1:2],rep(1,4),neo[,3:4])
neo = cbind(neo[,1:2],rep(1,4),neo[,3:4])
neo = neo[1:3,c(1,2,4,5)]
## Names
colnames(neo) = c("gene1","gene2","gene3","gene4")
rownames(neo) = c("sample1","sample2","sample3")
neo
neo["sample2","gene3"]
neo["sample2","gene3"] = 33
neo[c("sample1","sample3"),c("gene2","gene4")]
neo[c("sample1","sample3"),c("gene2","gene4")] = matrix(82:85,2,2)

## Some functions
head(neo)
mean(neo)
sum(neo) / length(neo)
neo
neo - 1
neo - neo
neo + neo

## Data structure exercice.
mat = matrix(runif(100*4),100,4)
colnames(mat) = c("sampleA","sampleB","sampleC","sampleD")
mat[,1] = mat[,1] + 2
mat[,2] = mat[,2] * 4
head(mat)
mean(mat[,1])
mean(mat[,2])
mean(mat[,3])
mean(mat[,4])
colnames(mat)[3] ## If the highest mean was found for mat[,3]
max(mat[,1])
max(mat[,2])
max(mat[,3])
max(mat[,4])
colnames(mat)[1] ## If the highest max was found for mat[,1]

## Apply
matExample = matrix(runif(100*100),100,100)
medCols = apply(matExample,2,median)
min(medCols)
max(medCols)

## Import/Export data
## In text format
write.table(matrixToWrite,file="textFile.txt", col.names=TRUE,row.names=FALSE, quote=FALSE, sep="\t")
input.data = read.table(file="textFile.txt", as.is=TRUE, header=TRUE, sep="\t")
## Exercice
mat.ge = read.table("dataForBasicPlots.tsv", as.is=TRUE, sep="\t", header=TRUE)
nrow(mat.ge) ## number of genes
ncol(mat.ge) ## number of samples
mat.ge[1:5,1:5]
## In R format
save(neo, matExample, file="someStuff.RData")
save.image(file="allTheStuff.RData")
load("someStuff.RData")
print(load("someStuff.RData"))  ## Load and print what was loaded !
load("allTheStuff.RData")
ls()

## Plots
hist(mat.ge[,1])
plot(mat.ge[,1],mat.ge[,2])

par(mfrow=c(1,2)) ## To print different plots in the same windows
hist(mat.ge[,1])
plot(mat.ge[,1],mat.ge[,2])
dev.off() ## close the window

hist(mat.ge[,1])
x11() ## To open a new window for a second plot
plot(mat.ge[,1],mat.ge[,2])

## Plot exercise
geneMed = apply(mat.ge,1,median)
hist(geneMed,xlab="gene expression",ylab="number of genes",main="Average gene expression across the samples")
abline(v=mean(geneMed),lty=2)

sampMed = apply(mat.ge,2,median)
hist(sampMed,xlab="gene expression",ylab="number of samples",main="Average gene expression across the genes")
mat.ge.noOutlier = mat.ge[,-which.min(sampMed)]
sampMed.noOutlier = apply(mat.ge.noOutlier,2,median)
hist(sampMed.noOutlier,xlab="gene expression",ylab="number of samples",main="Average gene expression across the genes")

mat.ge = as.matrix(mat.ge)
plot(mat.ge["gene333",],mat.ge["gene666",],xlab="gene333",ylab="gene666",main="Correlated genes")
lines(mat.ge["gene333",],mat.ge["gene667",],type="p",col=2,pch=2)

boxplot(mat.ge[,1:10],xlab="sample",ylab="expression",main="Box plot example")


## Conditions examples
test = 2 + 2 == 4
test
!test
2>3 & 5<10
2>3 | 5<10
5:10 > 6
which(5:10 > 6)
5:10 >= 6 & 5:10<=7
which(5:10 >= 6 & 5:10<=7)

## Exercise
randVec = sample(0:10, 30, TRUE)
randVec = round(runif(30, -0.5, 10.5))
randVec = randVec[which(randVec>=3)]
randVec[which(randVec>8)] = 8

## Remove genes with median lower than 1
mat.ge = mat.ge[which(geneMed>=1),]

## If block
classGene1 = "unknown"
if(median(mat.ge[1,]) > 4){
  classGene1 = "high"
} else {
  classGene1 = "low"
}

classGene1 = ifelse(median(mat.ge[1,]) > 4, "high", "low") ## Other way

## Functions
## Example
clean.vec.fun <- function(x){
  x = x[which(x>=3)]
  x[which(x>8)] = 8
  return(x)
}
vec110 = 1:10
vec110.cleaned = clean.vec.fun(vec110)

## Exercise
classify <- function(v,low.th=3,high.th=7){
    v.med = median(v)
    if(v.med < low.th){
        return("small")
    } else if(v.med > high.th){
        return("high")
    } else {
        return("medium")
    }
}
classify(1:3)
classify(1:10)
classify(7:10)

aveMinMax <- function(input){
  return(mean(c(min(input), max(input))))
}
higher3 <- function(input){
  return(length(which(input>3)))
}
## Other way
higher3 <- function(input){
  return(sum(input>3))
}
aveMinMax(1:3)
higher3(1:10)
higher3(7:10)

geneAveMinMax = apply(mat.ge, 1, aveMinMax)

## Final exercise
load("metadata.RData")
head(groups)
geDiff <- function(geneExp){
  geCases = geneExp[which(groups=="case")]
  geControls = geneExp[which(groups=="control")]
  return(geCases - geControls)
}
mat.ge.diff = apply(mat.ge, 1, geDiff)
hist(mat.ge.diff)

## One-liner quizz
mean(apply(mat,2,mean)>0)
colnames(mat)[which.max(apply(mat,2,max))]
mat[which(apply(mat,1,function(r)all(r>0))),]
