setwd("C:/Users/weiga/Dropbox/QiuDi/GBS-Wu/")
#library(ape)
#library(genoPlotR)
library(tidyverse)

#####################
# 5/12/2023, n=256 isolates
#####################

library(ggtree)
library(readxl)
#tree <- treeio::read.tree("contree-257.dnd4")
tree <- treeio::read.tree("contree-v2.dnd")
ggtree(tree, layout = "fan")  +  geom_treescale() 

dd <- read_xlsx("gbs-pheno.xlsx", col_types = "text")
tips <- tibble(id=tree$tip.label)
tips <- tips %>% left_join(dd, c("id" = "node")) # direct reading didn't work
col <- c("green", "red", "blue")

p1 <- ggtree(tree, layout = "fan")  %<+% tips +
  geom_treescale() + 
  geom_tippoint(aes(color = virulence)) +
#  geom_tiplab(align = TRUE, size = 2) +
  geom_tiplab(align = TRUE, size = 2.7, aes(color = factor(wgs))) +
  scale_color_manual(values = col)  +
  theme(legend.position = c(0.6,0.4), 
        legend.title = element_blank(), # no title
        legend.key = element_blank()) # no keys

p1



####################
snp <- read_tsv("snp-locus.tsv")
anno <- annotation(x1 = snp$position, text = str_replace(snp$locus, "CDH81_", ""), rot=rep(45, 43))

#anno.others <- annotation(x1=snp$position, text=NA)

annos <- list(anno, anno, anno, anno, anno, anno, anno)

# change original file (with awk) to fit with tree order, to ensure comparison properly aligned:
bbone <- read_mauve_backbone("g7-mauve.tsv3")
names <- c("B111", "B507", "B509", "A909", "B508", "B105", "COH1")
names(bbone$dna_segs) <- names
tree <- newick2phylog("((B111:0.34744,(B507:0.03096,(B509:0.01022,A909:0.1094):0.03032):0.08275):0.05592,(B508:0.18249,(B105:0.01654, COH1:0.00228):0.31197)1.000:0.01669);")

#tree <- newick2phylog("((B508:0.18249,(B105:0.01654, COH1:0.00228):0.31197):0.01669, (B111:0.34744,(B507:0.03096,(B509:0.01022,A909:0.1094):0.03032):0.08275):0.05592);")

plot(tree)


bbone$dna_segs <- bbone$dna_segs[names(tree$leaves)]

 bbone$comparisons[[1]]$col <- "gray"
 bbone$comparisons[[2]]$col <- "gray"
 bbone$comparisons[[3]]$col <- "gray"
 bbone$comparisons[[4]]$col <- "gray"
 bbone$comparisons[[5]]$col <- "gray"
 bbone$comparisons[[6]]$col <- "gray"

#color by length of comparison
for (i in 1:length(bbone$comparisons)){
cmp <- bbone$comparisons[[i]]
bbone$comparisons[[i]]$length <- abs(cmp$end1 - cmp$start1) + abs(cmp$end2 - cmp$start2)
}

# color LCB (local colinear block): frame black
for (i in 1:length(bbone$dna_segs)) {
  bbone$dna_segs[[i]]$fill <- bbone$dna_segs[[i]]$col
  bbone$dna_segs[[i]]$col <- 'black'
}

# change color for an LCB, e.g., first one to orange 
# for (i in 1:length(bbone$dna_segs)) { # for each genome
#  bbone$dna_segs[[i]][1,'fill'] <- "orange"  # change color
#  bbone$dna_segs[[i]][3,'col'] <- "black" # change frame color
#}
 
plot_gene_map(dna_segs=bbone$dna_segs, comparisons=bbone$comparisons, tree=tree, global_color_scheme=c("length", "increasing", "gray", 0.7),  override_color_schemes=TRUE, annotations = anno, annotation_height=3, dna_seg_scale=TRUE, tree_scale = T)

#plot_gene_map(dna_segs=bbone$dna_segs, comparisons=bbone$comparisons,tree=tree)


#################################
t <- read.tree(text = "((((B112:0.0,B111:0.0):0.22579,B503:0.11224)1.000:0.12165,((B507:0.00355,B506:0.06971)1.000:0.02741,((B509:0.01022,(B126:0.04505,A909:0.05069)1.000:0.05871)1.000:0.00486,(B217:0.00205,(B504:0.00399,B505:0.00478)1.000:0.00322)1.000:0.01804)1.000:0.02546)1.000:0.08275)1.000:0.05592,((B136:0.02937,(B514:0.01688,(B508:0.04824,(B124:0.00418,((B211:0.00461,B501:0.00262)1.000:0.00979,((B517:0.00280,(B522:0.00426,(B510:0.00508,B518:0.00115)1.000:0.00105)1.000:0.00248)0.999:0.00067,(B502:0.00315,B516:0.00418)0.829:0.00030)1.000:0.00083)0.996:0.00094)1.000:0.03390)1.000:0.04866)1.000:0.03401)1.000:0.05158,((B521:0.10950,(B511:0.01204,B512:0.00304)1.000:0.04145)1.000:0.19599,(COH1:0.00228,(B523:0.01270,(((B120:0.00327,(B128:0.00268,((B105:0.00190,B106:0.00086)1.000:0.00277,B138:0.00222)0.742:0.00043)0.994:0.00071)1.000:0.00209,(B520:0.00399,(B205:0.00239,B515:0.00356)0.998:0.00056)0.794:0.00018)1.000:0.00288,(B519:0.00211,((B116:0.00177,((B109:0.00208,B110:0.00257)1.000:0.00194,(B103:0.00293,B104:0.00287)1.000:0.00080)1.000:0.00117)1.000:0.00086,(B114:0.00291,((B118:0.00162,(B209:0.00205,(B122:0.00138,B134:0.00240)0.688:0.00037)1.000:0.00177)0.991:0.00067,(B130:0.00177,(B140:0.00161,B219:0.00142)0.791:0.00050)0.457:0.00015)0.250:0.00032):0.00018)0.889:0.00045)0.993:0.00221)1.000:0.00328)1.000:0.00248)1.000:0.23079)1.000:0.08118)1.000:0.01669);")
t <- read.tree("gbs50.dnd2")
tip.color = c(2,2,1,1,1,1,2,3,2,1,
              1,2,1,1,2,2,1,1,1,1,
              1,1,1,1,1,1,4,1,2,2,
              2,2,2,1,2,1,1,2,2,2,
              2,2,2,2,2,2,2,2,2,2)
plot(t, type = "p", font =2, tip.color = tip.color, no.margin = T, edge.width = 2)
add.scale.bar(lwd=2)

gbs.table <- read.csv(header = T, file = "gbs48.csv", row.names = 1)
name.tree <- scan("id-list.txt", what = "c")
gbs.table <- gbs.table[name.tree,]

t2 <- read.tree("batch-cc17/cc17.dnd")
t2 <- root(t2, "COH1")
tip.color <- c(2,2,2,2,2,2,2,
               2,2,1,2,1,1,2,
               1,2,2,2,2,2,2,
               2,2,2
               )
plot(t2, type = "p", font =2, no.margin = T, tip.color = tip.color, edge.width = 2)

