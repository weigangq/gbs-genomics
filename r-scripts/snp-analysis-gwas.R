# gbs GWAS analysis with SNPs
library(tidyverse)
library(readxl)
setwd("/Users/weiga/Dropbox/gbs-genomics/")

# read pheno
groups <- read_xlsx("data/gbs-pheno.xlsx")
pheno <- groups %>% 
  filter(node != 'B111') %>%
  select(1,2) %>% rename(iso = node)

# read geno
traw <- read_tsv("data/region_genotype.traw")
traw.long <- traw %>% 
  select(4,7:261) %>% 
  pivot_longer(2:256, names_to = "iso", values_to = "geno", values_drop_na = T)
traw.long <- traw.long %>% 
  mutate(iso = str_remove(iso, "^.+_")) %>% 
  mutate(geno = if_else(geno == 2, 1, 0))

allele <- traw %>% select(4:6) %>% 
  pivot_longer(2:3, names_to = "type", values_to = "nt") %>% 
  mutate(geno = if_else(type == 'COUNTED', 0, 1))

traw.long <- traw.long %>% left_join(allele, c("POS", "geno")) 
traw.long <- traw.long %>% select(-4)
x <- traw.long %>% left_join(pheno, "iso")

# add csq, containing redundant positions!!
csq <- read_tsv("data/snps-anno.tsv", col_names = F)
csq <- csq %>% mutate(X2 = str_remove(X2, '\\*'))
csq <- csq %>% select(-4)
colnames(csq) <- c("POS", "mut.type", "locus", "strand", "aa", "nt")
# make unique positions
ids <- csq %>% group_by(POS) %>% count()
duplicated_ids <- ids %>% filter(n > 1)
csq.distinct <- distinct(csq, POS, .keep_all = T)

# consistency index & homoplasy
library(ape)
library(phangorn)
library(TreeSearch)

m <- traw.long %>% select(1,2,3) %>% mutate(geno = as.character(geno)) %>% 
  pivot_wider(names_from = POS, values_from = geno, values_fill = "?")
mat <- as.matrix(m[,2:44619])
rownames(mat) <- m$iso
tree <- read.tree("data/snp-tree.dnd2")
mat2 <- mat[tree$tip.label,]
my_data <- phyDat(mat2, type = "USER", levels = c('0','1'))
ci <- CI(tree, my_data, sitewise = T)
df.ci <- tibble(POS = as.integer(colnames(mat)), ci = ci)
df.ci %>% ggplot(aes(ci)) +
  geom_histogram(bins = 100) + 
  theme_bw()

df.ci %>% ggplot(aes(POS, ci)) +
  geom_line(alpha = 0.1) +
#  geom_point(shape = 1, alpha = 0.1) + 
  theme_bw()

# contigency test (to get homoplay scores)
snp_ct <- x %>%   group_by(POS, nt, virulence) %>%   count()

# get pos with < 4 entries:
bad_pos <- snp_ct %>% group_by(POS) %>% count() %>% filter(n < 4)

library(broom)

# batch fisher's exact test
snp_fisher <- x %>% 
  filter(!POS %in% bad_pos$POS) %>%   
  group_by(POS) %>%   
  do(tidy(xtabs(~ nt + virulence, data = .) %>% 
            fisher.test(simulate.p.value = T)))

# takes ~10 min
y <- snp_fisher %>% left_join(csq.distinct, "POS")

y <- y %>% mutate(log.p = -log10(p.value))
y <- y %>% left_join(df.ci, "POS")
y <- y %>% mutate(homoplasy = if_else(ci <= 0.25, "0-Hi-hmp",
                                      if_else(ci < 1, "1-Med-hmp", "2-Consistent"))
                  )
#save(y, file = "RData-saved/snps-y.RData")
# load("data/snps.RData")
# save(x, file = "RData-saved/snps-long.RData")

# plot vocano
y %>% filter(mut.type %in% c("missense", "synonymous")) %>% 
  ggplot(aes(x = estimate, y = log.p, color = homoplasy)) +   
  geom_point(shape = 1, alpha = 0.1) +   
  scale_x_log10() + 
  geom_vline(xintercept = 1, linetype = 2) +  
  theme_bw() +  
  xlab("odds ratio (log10)") +  
  ylab("signficance (-log10[p])") + 
  facet_grid(rows =  vars(mut.type), cols =  vars(homoplasy)) + 
  theme(legend.position = "none")

y %>% filter(mut.type %in% c("missense", "synonymous")) %>% 
  group_by(mut.type, homoplasy) %>% count()

# Manhattan plots (significance)
y %>% filter(mut.type %in% c("missense", "synonymous")) %>% 
  ggplot(aes(x = POS, y = log.p, color = homoplasy)) +   
  geom_point(shape = 1, alpha = 0.1) +
  geom_hline(yintercept = -log10(1.25e-7), linetype = 2) +
  xlab("SNP position") +  
  ylab("signficance (-log10[p])") + 
  theme_bw() +
  facet_grid(rows =  vars(mut.type), cols =  vars(homoplasy)) + 
  theme(legend.position = "none") 

# Manhattan plots (homoplasy)
y %>% filter(mut.type %in% c("missense", "synonymous")) %>% 
  ggplot(aes(x = POS, y = ci, color = homoplasy)) +   
  geom_line(alpha = 0.5) +   
  xlab("SNP position") +  
  ylab("Consistency") + 
  theme_bw() +
  facet_grid(rows =  vars(mut.type), cols =  vars(homoplasy)) + 
  theme(legend.position = "none") 

# Gene effects
y %>% filter(!is.na(locus)) %>% 
  filter(mut.type == 'missense') %>% 
  group_by(locus) %>% count() %>% arrange(desc(n)) %>% pull(n) %>% stem()

##########################################
# binary PGLMM
# Ref: https://rdrr.io/cran/ape/man/binaryPGLMM.html
y %>% filter(log.p > 15, homoplasy == '2-Consistent')

pos = 293668  # positive control
pos = 163757 # negative control (hi homoplasy)
pos = 135 # negative control (consistent)

dat <- x %>% filter(POS == pos) %>% 
  select(2,3,5) %>% mutate(Y = if_else(virulence == 'colonizing', 0, 1)) %>% 
  as.data.frame()

rownames(dat) <- dat$iso 
dat <- dat[tree$tip.label,]

# Fit model, p = 0.1128 (ns), works well
out <- binaryPGLMM(Y ~ geno, phy = tree, data = dat)

# glm only, p < 1e-16, highly biased
summary(glm(Y ~ geno, data = dat))

# compare with phylolm, negative control: p = 1.4e-5, less biased
library(phylolm)
summary(phyloglm(Y ~ geno, phy = tree, data = dat))

#########################
# Trait plots, with ggtree (not working; better with heatmap & customized row cluster with phylogenetic tree (see above)
library(ggtree)

pos <- c(135, 163757, 293668)

snp_dat <- x %>% filter(POS %in% pos) %>% 
  select(2,3,5) %>% mutate(Y = if_else(virulence == 'colonizing', 0, 1))

# Create ggtree object
p <- ggtree(tree) 
tips <- tibble(id = tree$tip.label)
tips <- tips %>% left_join(pheno, c("id" = "iso"))
col <- c("green", "red")
p <- p %<+% tips +
  geom_tippoint(aes(color = virulence)) +
  scale_color_manual(values = col)
p
###########################
# PGLMM for all significant SNPs
# Trait plot with pheatmap
y.can1 <- y %>% filter(log.p > -log10(1.25e-7) & homoplasy == '0-Hi-hmp') # cutoff: BH correction: 0.05/39858
y.can2 <- y %>% filter(log.p > -log10(1.25e-7) & homoplasy == '1-Med-hmp') 
y.can3 <- y %>% filter(log.p > -log10(1.25e-7) & homoplasy == '2-Consistent') 

geno.long1 <- x %>% filter(POS %in% y.can1$POS)
l.can1 <- split(geno.long1, geno.long1$POS)

geno.long2 <- x %>% filter(POS %in% y.can2$POS)
l.can2 <- split(geno.long2, geno.long2$POS)

geno.long3 <- x %>% filter(POS %in% y.can3$POS)
l.can3 <- split(geno.long3, geno.long3$POS)

test.pglmm <- function(x) {
  cat("processing ...", x$POS[1], "\t")
  #head(x)
  dat <- x %>% 
    select(2,3,5) %>% 
    mutate(Y = if_else(virulence == 'colonizing', 0, 1)) %>% 
    as.data.frame()
  
  taxa.full <- tree$tip.label
  taxa.present <- dat$iso
  
  if(length(taxa.present) < length(taxa.full)) {
    taxa.absence <- tree$tip.label[!(tree$tip.label %in% x$iso)]
    tr <- drop.tip(tree, taxa.absence)
    cat("missing isolates", length(taxa.absence))
  } else {
    tr <- tree
  }
  cat("\n")  
  rownames(dat) <- dat$iso 
  dat <- dat[tr$tip.label,]

  # Fit model, p = 0.1128 (ns), works well
  binaryPGLMM(Y ~ geno, phy = tr, data = dat)
}  

pglmm.out1 <- lapply(l.can1, function(x) test.pglmm(x))
#save(pglmm.out1, file = "pglmm-output-hi-homoplasy.RData")
#pglmm.out2 <- lapply(l.can2, function(x) test.pglmm(x))
#save(pglmm.out2, file = "pglmm-output-med-homoplasy.RData")
load("RData-saved/snps.RData")
load("RData-saved/snps-long.RData")

# all consistent SNPs, not worth doing
# pglmm.out3 <- lapply(l.can3, function(x) test.pglmm(x))

pglmm2df <- function(x) {
  df <- data.frame(B = x$B['geno',],
                   Bse = x$B.se['geno',],
                   Bzscore = x$B.zscore['geno',],
                   Bpval = x$B.pvalue['geno',]
                    )
  df
}

# pglmm output: take long: ~400 positions per hour
#pg1.out <- lapply(pglmm.out1, function(x) pglmm2df(x))
#pg.out <- lapply(pglmm.out1, function(x) pglmm2df(x))

############
# combine PGLMM data & save object
#load("RData-saved/pglmm-output-hi-homoplasy.RData")
#pg1.out <- lapply(pglmm.out1, function(x) pglmm2df(x))
#pg1.df <- bind_rows(pg1.out)
#pg1.df <- pg1.df %>% mutate(pos = as.integer(names(pglmm.out1)))
#rm(pglmm.out1) # off-load, too large
#save(pg1.df, file = "RData-saved/pglmm-sig-hi-homoplasy-df.RData")

#load("RData-saved/pglmm-output-med-homoplasy.RData")
#pg2.out <- lapply(pglmm.out2, function(x) pglmm2df(x))
#pg2.df <- bind_rows(pg2.out)
#pg2.df <- pg2.df %>% mutate(pos = as.integer(names(pglmm.out2)))
#rm(pglmm.out2) # off-load, too large
#save(pg2.df, file = "RData-saved/pglmm-sig-med-homoplasy-df.RData")

#load("RData-saved/pglmm-output-consistent.RData")
#pg3.out <- lapply(pglmm.out3, function(x) pglmm2df(x))
#pg3.df <- bind_rows(pg3.out)
#pg3.df <- pg3.df %>% mutate(pos = as.integer(names(pglmm.out3)))
#rm(pglmm.out3) # off-load, too large
#save(pg3.df, file = "RData-saved/pglmm-sig-consistent-df.RData")

################################
# from wallace "gbs-pglmm" folder (calculated from r_env)
load("RData-saved/pglmm-sig-hi-homoplasy-df.RData") # containing "pg1.df"
load("RData-saved/pglmm-sig-med-homoplasy-df.RData") # containing "pg2.df"
load("RData-saved/pglmm-sig-consistent-df.RData") # containing "pg3.df"
load("RData-saved/pglmm.ns.hi-hmp-df.RData") # containing "pg1.ns.df"
load("RData-saved/pglmm.ns.med-hmp-df.RData") # containing "pg2.ns.df"
load("RData-saved/pglmm.ns.consistent-df.RData") # containing "pg3.ns.df"

pg.df <- bind_rows(pg1.df, pg2.df, pg3.df, pg1.ns.df, pg2.ns.df, pg3.ns.df)

load("RData-saved/snps-y.RData") # get ""y"
y <- y %>% left_join(pg.df, c("POS" = "pos"))
y <- y %>% mutate(log.p.lmm = -log10(Bpval))
y <- y %>% mutate(mut.type = if_else(is.na(mut.type), "igs", mut.type))
#save(y, file = "RData-saved/snps-pglmm.RData")

# Plot 1. p val comparisons
library(ggrepel)
y.sig <- y %>% filter(mut.type %in% c('missense', 'synonymous', 'igs')) %>%
  filter(log.p.lmm > 3)

y %>% filter(mut.type %in% c('missense', 'synonymous', 'igs')) %>% 
  ggplot(aes(log.p, log.p.lmm)) +
  geom_point(shape =1, alpha = 0.5 ) +
#  facet_wrap(~homoplasy) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_hline(yintercept = 3, linetype = 2, color = "red") +
  geom_text_repel(data = y.sig, aes(label = POS), max.overlaps = 100) +
  facet_grid(rows =  vars(mut.type), cols =  vars(homoplasy)) + 
  theme_bw()

# Plot 2. Manhattan plots (homoplasy)
y %>% filter(mut.type %in% c("missense", "synonymous", "igs")) %>% 
  ggplot(aes(x = POS, y = log.p.lmm, color = mut.type)) +   
  geom_point(alpha = 0.5, color = "gray") +   
  xlab("SNP position") +  
#  ylab("") + 
  theme_bw() +
  geom_hline(yintercept = 3, linetype = 2, color = "red") +
  geom_point(data = y.sig) +
  geom_text_repel(data = y.sig, aes(label = POS), max.overlaps = 100) +
  scale_color_manual(values = c("black", "red", "blue")) +
#  facet_grid(rows =  vars(mut.type), cols =  vars(homoplasy)) + 
  theme(legend.position = "bottom") 

# Plot 3. vocano
y %>% filter(mut.type %in% c("missense", "synonymous", "igs")) %>% 
  ggplot(aes(x = estimate, y = log.p)) +   
  geom_point(shape = 1, alpha = 0.1) +   
  scale_x_log10() + 
  geom_vline(xintercept = 1, linetype = 2) +  
  theme_bw() +  
  xlab("odds ratio (log10)") +  
  ylab("signficance (-log10[p])") + 
  facet_grid(rows =  vars(mut.type), cols =  vars(homoplasy)) + 
  theme(legend.position = "none")

y %>% filter(mut.type %in% c("missense", "synonymous", "igs")) %>% 
  ggplot(aes(x = estimate, y = log.p.lmm)) +   
  geom_point(shape = 1, alpha = 0.1) +   
  scale_x_log10() + 
  geom_vline(xintercept = 1, linetype = 2) +  
  theme_bw() +  
  xlab("odds ratio (log10)") +  
  ylab("signficance (-log10[p])") + 
  facet_grid(rows =  vars(mut.type), cols =  vars(homoplasy)) + 
  theme(legend.position = "none")

###################

pg.df <- pg.df %>% arrange(Bpval)
pg.df <- pg.df %>% left_join(y, c("pos" = "POS"))
pg.df <- pg.df %>% mutate(log.p.lmm = -log10(Bpval))

# p.vals (with and w/o correction) are not very correlated
pg.df %>% ggplot(aes(log.p, log.p.lmm)) +
  geom_point(shape =1 ) +
  theme_bw()

# Manhattan plot
library(ggrepel)
sig <- pg.df %>% filter(log.p.lmm > -log10(0.05))
pg.df %>% ggplot(aes(pos, log.p.lmm, color = mut.type)) +
  geom_jitter(alpha = 0.5) +
  geom_text_repel(data = sig, aes(pos, log.p.lmm, label = pos, color = mut.type)) +
  theme_bw() +
  geom_hline(yintercept = -log10(0.05), linetype = 2) + 
  theme(legend.position = "bottom")

# plot effect sizes, waterfall plot
# 78 positions are duplicated (complex SNPs?). Make a unique factor:
#pg.df <- pg.df %>% mutate(pos.uniq = as.character(pos))
#pg.df$pos.uniq <- make.unique(pg.df$pos.uniq)

#t <- pg.df %>% arrange(B) %>% pull(pos.uniq)
#pg.df$pos.uniq <- factor(pg.df$pos.uniq, levels = t)

pg.df %>% 
  ggplot(aes(y = B, x = pos, color = mut.type)) + 
  geom_errorbar(aes(ymin = B - Bse, ymax = B + Bse), alpha = 0.2) + 
  geom_point(shape = 1) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw() + 
  coord_flip()

# heatmap of top SNPs
library(pheatmap)
top <- pg.df %>% filter(log.p.lmm > 3)
#top.df <- as.data.frame(mat2[,top$pos.uniq])
top.df.long <- x %>% filter(POS %in% top$pos) %>% select(POS, iso, geno)

#top.df.long <- top.df %>% 
#  mutate(iso = rownames(top.df)) %>% 
#  pivot_longer(1:61, names_to = "pos", values_to = "state")
#top.df.long  <- top.df.long %>% mutate(pos = as.integer(pos), state = as.integer(state))

top.wide <- top.df.long %>% 
  pivot_wider(names_from = "POS", values_from = "geno")

# pad 0's for easy sorting
colnames(top.wide)[2:61] <- str_pad(colnames(top.wide)[2:61], width = 7, side = "left", pad = "0")
top.df <- data.frame(top.wide[,2:61], row.names = top.wide$iso)
#col_name <- sort(colnames(top.df))
#top.df <- top.df[,col_name]

# use custom row tree
# ref: https://stackoverflow.com/questions/62915975/using-heatmaply-dendrogram-sorting-in-pheatmap

library(phytools)
library(seriation)
library(dendextend)

row_dist <- as.dist(cophenetic.phylo(tree))
row_hc <- as.hclust(multi2di(force.ultrametric(tree)))
row_dend <- as.dendrogram(row_hc)
row_dend <- seriate_dendrogram(row_dend, row_dist, method="OLO")

anno_row <- data.frame(iso_type = factor(pheno$virulence), row.names = pheno$iso)

pos <- as.integer(str_remove(colnames(top.df), "^X0*" ))
anno_col1 <- data.frame(sites = pos, names = colnames(top.df)) %>% left_join(csq.distinct, c("sites" = "POS"))

anno_col1 <- anno_col1 %>% mutate(mut.type = if_else(str_detect(mut.type, 'NA'), "igs", mut.type))

anno_col <- data.frame(mut_type = anno_col1$mut.type, row.names = anno_col1$names)

anno_color = list(iso_type = c(colonizing = "cyan", invasive = "magenta"),
                  mut_type = c(synonymous = "green", missense = "red", igs = "gray" )
                  )


pheatmap(top.df, 
         cluster_rows = as.hclust(row_dend), 
         cluster_cols = F, show_rownames = F,
         annotation_row = anno_row,
         annotation_col = anno_col,
         annotation_colors = anno_color,
         cutree_rows = 7,
         color = c("lightpink", "lightblue"), breaks = c(-1, 0.5, 1.5),
         cellwidth = 12, cellheight = 2.5, 
         legend = F
         )

# replot vocano and manhattan (to compare GWAS)
# plot vocano
pg.df %>% filter(mut.type %in% c("missense", "synonymous")) %>% 
  ggplot(aes(x = estimate, y = log.p)) +   
  geom_point(shape = 1, alpha = 0.2) +   
  scale_x_log10() + 
  geom_vline(xintercept = 1, linetype = 2) +
  geom_hline(yintercept = c(3,7), linetype = 2) +
  geom_point(aes(x = estimate, y = log.p.lmm), color = "red", alpha = 0.2) +
  geom_text_repel(data = pg.df %>% filter(log.p.lmm > 3), aes(x = estimate, y = log.p.lmm, label = pos), color = "red", alpha = 0.5) +
  theme_bw() +  
  xlab("odds ratio (log10)") +  
  ylab("signficance (-log10[p])") + 
  theme(legend.position = "none")

# output locus info
top %>% arrange(pos) %>% select(pos, locus, mut.type, aa, nt)
write_tsv(top, file = "top-snps.tsv")

##################################
# consistent & significant SNPs
#################################

top.con <- y.can3 %>% filter(!is.na(locus)) %>% 
  filter(mut.type %in% c("missense", "synonymous")) %>% 
  group_by(locus, mut.type) %>% count()

top.con.wide <- top.con %>% 
  pivot_wider(names_from = mut.type, values_from = n, values_fill = 0)

top.con.loci <- top.con %>% filter(mut.type == 'missense' & n >= 5 ) %>% pull(locus)
top.con.pos <- y.can3 %>% filter(locus %in% top.con.loci) %>% pull(POS)
  
top.con %>% ggplot(aes(locus, n)) +
  geom_col() +
  facet_wrap(~mut.type, nrow = 1) + 
  theme_bw() +
  coord_flip()

top.con %>% ggplot(aes(locus, n)) +
  geom_col() +
  facet_wrap(~mut.type, nrow = 1) + 
  theme_bw() +
  coord_flip()

top.con.wide %>% ggplot(aes(synonymous, missense, label = locus)) +
  geom_jitter(shape = 1) +
  geom_text_repel() +
  geom_hline(yintercept = 3, linetype = 2) +
  theme_bw()

top.con.long <- x %>% filter(POS %in% top.con.pos) %>% select(POS, iso, geno)

# heatmap to verify
top.con.wide <- top.con.long %>% 
  pivot_wider(names_from = "POS", values_from = "geno")

# pad 0's for easy sorting
colnames(top.con.wide)[2:(ncol(top.con.wide)-1)] <- str_pad(colnames(top.con.wide)[2:(ncol(top.con.wide)-1)], width = 7, side = "left", pad = "0")
top.con.df <- data.frame(top.con.wide[,2:(ncol(top.con.wide)-1)], row.names = top.con.wide$iso)

pos <- as.integer(str_remove(colnames(top.con.df), "^X0*" ))

anno_col1 <- data.frame(sites = pos, names = colnames(top.con.df)) %>% left_join(csq.distinct, c("sites" = "POS"))

anno_col1 <- anno_col1 %>% mutate(mut.type = if_else(str_detect(mut.type, 'NA'), "igs", mut.type))

anno_col <- data.frame(mut_type = anno_col1$mut.type, row.names = anno_col1$names)

anno_color = list(iso_type = c(colonizing = "cyan", invasive = "magenta"),
                  mut_type = c(synonymous = "green", missense = "red", igs = "gray", stop_gained = "gray" )
)


pheatmap(top.con.df, 
         cluster_rows = as.hclust(row_dend), 
         cluster_cols = F, show_rownames = F,
         annotation_row = anno_row,
         annotation_col = anno_col,
         annotation_colors = anno_color,
         cutree_rows = 7,
         color = c("lightpink", "lightblue"), breaks = c(-1, 0.5, 1.5),
         #cellwidth = 5, cellheight = 2.5, 
         legend = F
)

######
load("RData-saved/pglmm-output-consistent.RData")
pg.df <- bind_rows(pglmm.out3)
pg.df <- pg.df %>% mutate(pos = as.integer(names(pglmm.out3)))
pg.df <- pg.df %>% arrange(Bpval)
pg.df <- pg.df %>% left_join(y, c("pos" = "POS"))
pg.df <- pg.df %>% mutate(log.p.lmm = -log10(Bpval))

# Manhattan plot
library(ggrepel)
sig <- pg.df %>% filter(log.p.lmm > 3)
pg.df %>% ggplot(aes(pos, log.p.lmm, color = mut.type)) +
  geom_point(alpha = 0.5) +
  geom_text_repel(data = sig, aes(pos, log.p.lmm, label = pos, color = mut.type)) +
  theme_bw() +
  geom_hline(yintercept = 3, linetype = 2) + 
  theme(legend.position = "bottom")

# plot effect sizes, waterfall plot
# 78 positions are duplicated (complex SNPs?). Make a unique factor:
pg.df <- pg.df %>% mutate(pos.uniq = as.character(pos))
pg.df$pos.uniq <- make.unique(pg.df$pos.uniq)

t <- pg.df %>% arrange(B) %>% pull(pos.uniq)
pg.df$pos.uniq <- factor(pg.df$pos.uniq, levels = t)

pg.df %>% 
  ggplot(aes(y = B, x = pos.uniq, color = mut.type)) + 
  geom_errorbar(aes(ymin = B - Bse, ymax = B + Bse), alpha = 0.2) + 
  geom_point(shape = 1) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme_bw() + 
  coord_flip()

# heatmap of top SNPs
library(pheatmap)
top <- pg.df %>% filter(log.p.lmm > -log10(0.05))
top.df.long <- x %>% filter(POS %in% top$pos) %>% select(POS, iso, geno)

top.con.wide <- top.df.long %>% 
  pivot_wider(names_from = "POS", values_from = "geno")

# pad 0's for easy sorting
colnames(top.con.wide)[2:(ncol(top.con.wide)-1)] <- str_pad(colnames(top.con.wide)[2:(ncol(top.con.wide)-1)], width = 7, side = "left", pad = "0")
top.con.df <- data.frame(top.con.wide[,2:(ncol(top.con.wide)-1)], row.names = top.con.wide$iso)


# use custom row tree
# ref: https://stackoverflow.com/questions/62915975/using-heatmaply-dendrogram-sorting-in-pheatmap

library(phytools)
library(seriation)
library(dendextend)

row_dist <- as.dist(cophenetic.phylo(tree))
row_hc <- as.hclust(multi2di(force.ultrametric(tree)))
row_dend <- as.dendrogram(row_hc)
row_dend <- seriate_dendrogram(row_dend, row_dist, method="OLO")

anno_row <- data.frame(iso_type = factor(pheno$virulence), row.names = pheno$iso)

pos <- as.integer(str_remove(colnames(top.con.df), "^X0*" ))
anno_col1 <- data.frame(sites = pos, names = colnames(top.con.df)) %>% left_join(csq.distinct, c("sites" = "POS"))

anno_col1 <- anno_col1 %>% mutate(mut.type = if_else(str_detect(mut.type, 'NA'), "igs", mut.type))

anno_col <- data.frame(mut_type = anno_col1$mut.type, row.names = anno_col1$names)

anno_color = list(iso_type = c(colonizing = "cyan", invasive = "magenta"),
                  mut_type = c(synonymous = "green", missense = "red", igs = "gray" )
)

library(pheatmap)
pheatmap(top.con.df, 
         cluster_rows = as.hclust(row_dend), 
         cluster_cols = F, show_rownames = F,
         annotation_row = anno_row,
         annotation_col = anno_col,
         annotation_colors = anno_color,
         cutree_rows = 7,
         color = c("lightpink", "lightblue"), breaks = c(-1, 0.5, 1.5),
         cellwidth = 12, cellheight = 2.5, 
         legend = F
)

# replot vocano and manhattan (to compare GWAS)
# plot vocano
pg.df %>% filter(mut.type %in% c("missense", "synonymous")) %>% 
  ggplot(aes(x = estimate, y = log.p)) +   
  geom_point(shape = 1, alpha = 0.2) +   
  scale_x_log10() + 
  geom_vline(xintercept = 1, linetype = 2) +
  geom_hline(yintercept = c(3,7), linetype = 2) +
  geom_point(aes(x = estimate, y = log.p.lmm), color = "red", alpha = 0.2) +
  geom_text_repel(data = pg.df %>% filter(log.p.lmm > 3), aes(x = estimate, y = log.p.lmm, label = pos), color = "red", alpha = 0.5) +
  theme_bw() +  
  xlab("odds ratio (log10)") +  
  ylab("signficance (-log10[p])") + 
  theme(legend.position = "none")

# output locus info
top %>% arrange(pos) %>% select(pos, locus, mut.type, aa, nt)
write_tsv(top, file = "top-snps.tsv")

########################
# combined plots (all 3 homoplasy levels)
# plot vocano
y %>% filter(mut.type %in% c("missense", "synonymous")) %>% 
  ggplot(aes(x = estimate, y = log.p)) +   
  geom_point(shape = 1, alpha = 0.2, color = "black") +   
  scale_x_log10() + 
  geom_vline(xintercept = 1, linetype = 2) +  
  theme_bw() +  
  geom_point(aes(x = estimate, y = -log10(Bpval)), color = "red", alpha = 0.5) +
  geom_hline(yintercept = 2, linetype = 2, color = "red") + 
  xlab("odds ratio (log10)") +  
  ylab("signficance (-log10[p])") + 
  facet_grid(rows =  vars(mut.type), cols =  vars(homoplasy)) +
  theme(legend.position = "none")


# Manhattan plots (significance)
y.tmp <- y %>% filter(mut.type %in% c("missense", "synonymous")) 
y.tmp <- y.tmp %>% mutate(part = if_else(POS < 7.5e5, "front-end", if_else(POS < 1.5e6, 'middle', 'back-end')))

y.tmp %>% 
  ggplot(aes(x = POS, y = -log10(Bpval), color = mut.type)) +   
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = -log10(1e-2), linetype = 2) +
  xlab("SNP position") +  
  ylab("signficance (-log10[p])") + 
  theme_bw() +
  geom_text_repel(data = y.tmp %>% filter(-log10(Bpval) > 2), aes(x = POS, y = -log10(Bpval), label = POS), alpha = 0.5) +
#  facet_grid(rows =  vars(mut.type), cols =  vars(homoplasy)) + 
  theme(legend.position = "none") +
  facet_wrap(~part, nrow = 3, scales = "free_x")

y.tmp2 <- y.tmp %>% filter(-log10(Bpval) > 2)

y <- y %>% mutate(p.val.lmm = -log10(Bpval))
y %>% filter(!is.na(Bpval)) %>% 
  ggplot(aes(x = log.p, y = p.val.lmm, color = mut.type)) +
  geom_point(shape = 1) +
  geom_hline(yintercept = 3, linetype = 2) +
  theme_bw() +
  theme(legend.position = "none")
  

##################################
# Compare with MCMCglm (Bayesian)
# doesn't work
library(MCMCglmm)

n <- 100
V <- vcv(tree)
V <- V/max(V)
detV <- exp(determinant(V)$modulus[1])
V <- V/detV^(1/n)

invV <- Matrix(solve(V),sparse=T)
rownames(invV) <- dat$iso

nitt <- 43000
thin <- 10
burnin <- 3000

prior <- list(R=list(V=1, fix=1), G=list(G1=list(V=1, nu=1000, alpha.mu=0, alpha.V=1)))

MCMCglmm(Y ~ geno, random = ~iso, ginvers = list(iso = invV),
                 data = dat, slice=F, nitt=nitt, thin=thin, burnin=burnin,
                 family="categorical", prior=prior, verbose = T)
########################################################
