### to import quant.sf file from Salmon analysis. 

library(tximport)
library(readr)

#create table with transcript number to IDs
genecode <- "C:/Users/Pierre/Desktop/IWP2 Paper/Third submission/RNAseq analysis/iTRanscriptome Full Dataset/gencode.vM17.metadata.MGI"
tx2gene <- read.csv(genecode, h=F, sep="\t") #on the server
colnames(tx2gene) <- c("TXNAME","GENEID")

#select all files

dir <- "C:/Users/Pierre/Desktop/IWP2 Paper/Third submission/RNAseq analysis/iTRanscriptome Full Dataset/Salmon/quants/"
setwd(dir)
files <- list.files(dir)
all(file.exists(files))

#import data
txi <- tximport(files = files, type="salmon", tx2gene = tx2gene, dropInfReps=T)
colnames(txi$counts) <- files
