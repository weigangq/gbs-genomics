# gbs snp analysis
library(tidyverse)
library(readxl)
#setwd("~/Dropbox/gbs-genomics/")
setwd("/Users/weiga/Dropbox/gbs-genomics/")

library(adegenet)
# ref/tutorial: https://adegenet.r-forge.r-project.org/files/tutorial-genomics.pdf
# Step 1.  Build a genlight object
x <- read_tsv("data/snp_mat.tsv", col_names = F)
ids <- scan("data/snp_sample_ids.txt", what = "c")
pos <- scan("data/snp_pos.txt")
x.gen <- new("genlight", t(x))
x.gen@ind.names <- ids
x.gen@position <- as.integer(pos)
g.size = 2150631

# step 2. plot SNP density & allele freq distribution
x.den <- density(position(x.gen), bw=10)
plot(x.den, type = "n", xlab = "Position", main = "SNP density")
polygon(c(x.den$x,rev(x.den$x)), c(x.den$y, rep(0,length(x.den$x))),
        col=transp("blue",.3))
points(position(x.gen), rep(0, nLoc(x.gen)), pch="|", col="blue")

# the following is nicer but takes too long
#snpposi.plot(position(x.gen), genome.size = g.size,  codon=FALSE)

# plot allele frequency distribution
x.freq <- glMean(x.gen)
hist(x.freq, proba=TRUE, col="gold", xlab="Allele frequencies",
     main="Distribution of alternative allele frequencies")
temp <- density(x.freq)
lines(temp$x, temp$y*1.8,lwd=3)

# PCA analysis; doesn't work. stuck for too long
pca1 <- glPca(x.gen, useC = T, nf = 30)

# discrimintant principal component analysis (DAPC)
# ref: https://github.com/thibautjombart/adegenet/wiki/Tutorials
grp <- find.clusters(x.gen, max.n.clust=15)

########################################
#smartsnp
# Reference: https://christianhuber.github.io/smartsnp/articles/mallard_smartpca_analysis.html
library(smartsnp)
groups <- read_xlsx("data/gbs-pheno.xlsx")
samples <- groups %>% 
  filter(node != 'B111') %>%
  pull(2)

pcaR <- smart_pca(snp_data = "data/snp_mat.tsv", 
                    sample_group = samples,
                    missing_value = NA)

# extract eigenvalues (PCA1 and PC2 axes)
pcaR_eigen <- pcaR$pca.eigenvalues 

# extract principal coefficients (high SNP loadings indicate loci with stronger variation across individuals)
pcaR_load <- pcaR$pca.snp_loadings

pcaR_load %>% ggplot(aes(PC1, PC2, label = rownames(pcaR_load))) +
  geom_point() +
  theme_bw()


# extract principal components (position of individuals in PCA space used to generate the ordination)
pcaR_coord <- pcaR$pca.sample_coordinates 

# plot pca
pcaR_coord %>% ggplot(aes(PC1, PC2, color = Group)) +
  geom_point(alpha = 0.5, size = 5) +
  theme_bw()

# make NJ tree
library(ape)
#x.mat <- read_tsv("data")
tre <- nj(dist(as.matrix(x.gen)))
plot(tre, type = "fan", cex = 0.5)

myCol <- colorplot(pcaR_coord[,3:4], pcaR_coord[,3:4], transp=T)
plot(tre, type = "fan", show.tip = F)
tiplabels(pch=20, col=myCol, cex=2)
