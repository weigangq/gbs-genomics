# adegenet tutorials
# https://github.com/thibautjombart/adegenet/wiki/Tutorials
# DAPC:
# https://raw.githubusercontent.com/thibautjombart/adegenet/master/tutorials/tutorial-dapc.pdf

library(tidyverse)
library(readxl)
#setwd("/Users/weiga/Dropbox/gbs-genomics/")
setwd("Dropbox/gbs-genomics/")
library(adegenet)
library(vcfR)

# read vcf and convert to genlight object
vcf <- read.vcfR("data/snps.vcf")
x.gl <- vcfR2genlight(vcf)

# add SNP positions
pos <- scan("data/snp_pos.txt")
x.gl@position <- as.integer(pos)

# add cut-tree 0.1 as ecological populations
grps <- read_tsv("data/cut-groups.tsv", col_names = F)
pop <- data.frame(clade = grps$X2, row.names = grps$X1)
x.gl@pop <- as.factor(pop[x.gl$ind.names,1])

# attach phenotype
groups <- read_xlsx("data/gbs-pheno.xlsx")
pheno <- groups %>% 
  filter(node != 'B111') %>%
  select(1,2)
x.gl@other <- list(data.frame(virulence = pheno$virulence, row.names = pheno$node))
g.size = 2150631

# plot SNP density
temp <- density(position(x.gl), bw=10)
plot(temp, type="n", xlab="Position in the alignment",
     main="Location of the SNPs")
polygon(c(temp$x,rev(temp$x)), c(temp$y, rep(0,length(temp$x))),
        col=transp("blue",.3))
points(position(x.gl), rep(0, nLoc(x.gl)), pch="|", col="blue")

# plot allele frequencies
myFreq <- glMean(x.gl)
hist(myFreq, proba=TRUE, col="gold", xlab="Allele frequencies",
     main="Distribution of (second) allele frequencies" , ylim = c(0,5))
temp <- density(myFreq)
lines(temp$x, temp$y*1.8,lwd=3)

# remove SNPs with missing data
snp.na <- apply(as.matrix(x.gl), 2, function(x) sum(is.na(x)))
snp.na2 <- apply(as.matrix(x.gl[-121,]), 2, function(x) sum(is.na(x))) # doesn't help
iso.na <- apply(as.matrix(x.gl), 1, function(x) sum(is.na(x)))

x.na <- glNA(x.gl) # count the number of NA for each SNPs
with.na <- x.na[which(x.na>0)] # with NA
no.na <- x.na[which(x.na == 0)] # no NA
plot(density(no.na))
plot(density(with.na))
# use only NA-free loci
y.gl <- x.gl[, which(x.na==0)]

# pca analysis
pca1 <- glPca(y.gl)
scatter(pca1, posi = "bottomright")


