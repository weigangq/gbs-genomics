# simulation study of bacterial GWAS
# 1. validation of binaryPGLMM
# 2. power analysis. with respect to: 
#    (1) sample size: fix at GBS size of n = 255
#    (2) tree shape [neutral vs gbs tree]
#    (3) effect size: p.case.minor & p.case.common
#    (4) homoplasy (SNP rates): high rates dilute phylo effect; only effective to remove high-consistent SNPs

library(tidyverse)
library(ape)
library(phangorn)
library(broom)

setwd("/Users/weiga/Dropbox/gbs-genomics/")

##########
# Model 0. neutral tree
# simulate neutral tree, 100 OTUs
#ntips <- 200
nsites <- 100
ntests <- 100
#tr <- rcoal(n <- ntips, br = runif)

# GBS tree
tree <- read.tree("data/snp-tree.dnd2")

# simulate lineage-associated SNPs, equal rates
# rownames(X) are the tip labels
# rate = 1e-3 for high consitency, 1e-2 for median, 1e-1 for hi homoplasy
rate_transition <- 1
list.X <- vector("list", length = ntests)
for(i in 1:100) {
  X <- replicate(nsites, rTraitDisc(tree, states = c(0,1), rate = rate_transition , model = "ER"))
# add colnames (SNP ids)
  colnames(X) <- paste0("rs", str_pad(1:nsites, side = "left", pad = "0", width = 3), sep = "")
  list.X[[i]] <- X
}

list.X <- setNames(list.X, paste("test", 1:100, sep = "_"))

# calculate CI for verification
cal.ci <- function(x, tree) {
  #x: tree-simulated SNPs
  my_data <- phyDat(x, type = "USER", levels = c('0','1'))
  ci <- CI(tree, my_data, sitewise = T)
  df.ci <- tibble(POS = colnames(x), ci = ci)
  return(df.ci)
}

df.ci %>% filter(!is.na(ci)) %>% 
  ggplot(aes(x = POS, y = ci)) +
  geom_col() + 
  theme_bw() + 
  coord_flip()

df.ci %>% filter(!is.na(ci)) %>% 
  ggplot(aes(ci)) +
  geom_histogram() + 
  theme_bw() 
# simulate causal SNP & binary phenotype
# Ref. Earle et al (2016) "Testing power by simulating pheonotypes". "To assess the performance of the method for controlling population structure, we performed 100 simulations per species. 
# In each simulation, a biallelic SNP was chosen randomly (from those SNPs with minor allele frequency above 20%) to be the causal SNP. 
# Binary phenotypes (case or control) were then simulated for each genome with case probabilities of 0.25 and 0.5, respectively, in individuals with the common and rare allele at the causal SNP (an odds ratio of 3). 
# For each simulated data set, we tested for locus effects at every biallelic SNP, and for lineage effects at every principal component, as described above. 
# The power to detect locus effects was defined as the proportion of simulations in which the causal SNP was found to have a significant locus effect. This was compared to a theoretically optimum power computed as the proportion of simulations in which the causal SNP was found to have a significant locus effect when population structure and multiple testing were not controlled for."

assign_causal_snp <- function(trait_mat) {
  # trait_mat are outputs from rTraitDisc
  # n: number of causal SNPs to simulate
  # with "0" as ancestral allele and "1" as derived/minor allele
  tips <- rownames(trait_mat)
  snps <- colnames(trait_mat)
  minor.ct <- apply(trait_mat, 2, function(x) table(x)['1'])
  minor.pool <- which(minor.ct >= 0.2 * nrow(trait_mat))
  try(if(length(minor.pool) < 10) stop("<10 SNPs to select; rerun rTraitDisc with higher rates"))
  pick.ind <- sample(minor.pool, 1)
  pick.hap <- trait_mat[,pick.ind]
  
  pheno <- character(length = length(tips))
  p.case.common <- 0.25
  p.case.minor <- 0.75
  for(i in 1:length(tips)) {
    if(pick.hap[i] == '1') { # minor allele
      pheno[i] <- rbinom(1, 1, p.case.minor)
    } else { # common allele
      pheno[i] <- rbinom(1, 1, p.case.common)
    }
  }
  colnames(trait_mat)[pick.ind] <- "pick"
  names(pheno) <- tips
  bind_cols(iso = rownames(trait_mat), phe = pheno, trait_mat)
  #return(list(mat = trait_mat, phe = pheno))
}

list.sims <- lapply(names(list.X)[1:2], function(name) {
  x.mat <- list.X[[name]]
  df.ci <- cal.ci(x.mat, tree)
  sim.out <- assign_causal_snp(x.mat)
  sim.long <- sim.out %>% 
    pivot_longer(3:(nsites+2), names_to = "pos", values_to = "geno")

  # contigency test without phylogenetic correction
  snp_ct <- sim.long %>%  group_by(pos, geno, phe) %>%   count()

  # get pos with < 4 entries:
  bad_pos <- snp_ct %>% group_by(pos) %>% count() %>% filter(n < 4)

# batch fisher's exact test
  snp_fisher <- sim.long %>% 
    filter(!pos %in% bad_pos$pos) %>%   
    group_by(pos) %>%   
    do(tidy(xtabs(~ geno + phe, data = .) %>% 
            fisher.test(simulate.p.value = T)))

  y <- snp_fisher %>% mutate(log.p = -log10(p.value))
  y <- y %>% mutate(causal = if_else(pos == 'pick', 1, 0))
  pick.pos <- colnames(x.mat)[!colnames(x.mat) %in% unique(snp_ct$pos)]
  y <- y %>% mutate(pos = if_else(pos == 'pick', pick.pos, pos))
  y <- y %>% left_join(df.ci, c("pos" = "POS"))
  y <- y %>% mutate(homoplasy = if_else(ci <= 0.25, "0-Hi-hmp",
                                        if_else(ci < 1, "1-Med-hmp", "2-Consistent")))
  
  sim.long <- sim.long %>% filter(!pos %in% bad_pos$pos)
  l.snps <- split(sim.long, sim.long$pos)
  pglmm.out <- lapply(l.snps, function(x) test.pglmm(x, tree))
  pg.out <- lapply(pglmm.out, function(x) pglmm2df(x))
  pg.df <- bind_rows(pg.out)
  pg.df <- pg.df %>% mutate(pos = names(pglmm.out))
  pg.df <- pg.df %>% mutate(pos = if_else(pos == 'pick', pick.pos, pos))
  y <- y %>% left_join(pg.df, "pos")
  y <- y %>% mutate(log.p.lmm = -log10(Bpval), test = name)
  return(y)
  } 
)

# plot vocano
y %>% 
  ggplot(aes(x = estimate, y = log.p, color = as.factor(causal))) +   
  geom_point(size = 3, alpha = 0.5) +   
  scale_x_log10() + 
  geom_vline(xintercept = 1, linetype = 2) +  
  theme_bw() +  
  xlab("odds ratio (log10)") +  
  ylab("signficance (-log10[p])") +
#  geom_point(data = y %>% filter(homoplasy == 'causal'), aes(estimate, log.p), color = "black", size = 5) +
  facet_wrap(~homoplasy)
#  facet_grid(rows =  vars(mut.type), cols =  vars(homoplasy)) + 
#  theme(legend.position = "none")

# with PGLMM
sim.long <- sim.long %>% filter(!pos %in% bad_pos$pos)
l.snps <- split(sim.long, sim.long$pos)

# input: long table "list of long tables" & tree
test.pglmm <- function(x, tree) {
  cat("processing ...", x$pos[1], "\t")
  #head(x)
  cat("\n")  
  dat <- data.frame(iso = x$iso, phe = as.integer(x$phe), pos = x$pos, geno = as.integer(x$geno))
  rownames(dat) <- x$iso
  binaryPGLMM(phe ~ geno, phy = tree, data = dat)
}  

# screen out incomplete contigency tables:
#l.snps2 <- l.snps[!names(l.snps) %in% bad_pos$pos]
pglmm.out <- lapply(l.snps, function(x) test.pglmm(x, tree))

pglmm2df <- function(x) {
  df <- data.frame(B = x$B['geno',],
                   Bse = x$B.se['geno',],
                   Bzscore = x$B.zscore['geno',],
                   Bpval = x$B.pvalue['geno',]
  )
  df
}

pg.out <- lapply(pglmm.out, function(x) pglmm2df(x))
pg.df <- bind_rows(pg.out)
pg.df <- pg.df %>% mutate(pos = names(pglmm.out))
pg.df <- pg.df %>% mutate(pos = if_else(pos == 'pick', pick.pos, pos))
y <- y %>% left_join(pg.df, "pos")
y <- y %>% mutate(log.p.lmm = -log10(Bpval))

# not very effective
y %>% filter(!is.na(log.p.lmm)) %>% 
  ggplot(aes(log.p, log.p.lmm, color = as.factor(causal))) +
  geom_point(size = 3, alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  theme_bw() + 
  facet_wrap(~homoplasy)

load("RData-saved/gwas-sim.RData")

list.sims2 <- lapply(list.sims, function(x) {
  ci.causal <- x %>% filter(causal == '1') %>% pull(ci) %>% nth(1);
  x <- x %>% mutate(ci.causal = ci.causal);
  x <- x %>% mutate(hmp.causal = if_else(ci.causal <= 0.25, "0-Hi-hmp",
                                         if_else(ci.causal < 1, "1-Med-hmp", "2-Consistent")));
  return(x)
  })

y.sim <- bind_rows(list.sims2)
# re-plot vocano
#library(ggrepel)

# by homoplasy
y.sim %>%
  ggplot(aes(log.p, log.p.lmm, color = as.factor(causal))) +
  geom_point(size = 3, alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  theme_bw() + 
  scale_color_manual(values = c("gray", "brown")) +
  geom_hline(yintercept = 4, linetype = 2, color = "red") +
  facet_wrap(~homoplasy) +
  theme(legend.position = "bottom") 



# FDR
false.fisher <- y.sim %>% 
  filter(log.p >= 4) %>% 
  xtabs(~ causal + homoplasy, data = .) %>% 
  as.data.frame() %>% 
  pivot_wider(names_from = 'causal', values_from = "Freq") 

false.pglmm <- y.sim %>% 
  filter(log.p.lmm >= 4) %>% 
  xtabs(~ causal + homoplasy, data = .) %>% 
  as.data.frame() %>% 
  pivot_wider(names_from = 'causal', values_from = "Freq") 


colnames(false.fisher) <- c("homoplasy", "lineage.snps", "causal.snps")
colnames(false.pglmm) <- c("homoplasy", "lineage.snps", "causal.snps")

fp.rate <- function(x) {
  fp <- x[1]/(x[1] + x[2])
  t <- binom.test(x[1], x[1] + x[2])
  lo <- t$conf.int[1]
  hi <- t$conf.int[2]
  return(c(fp, lo, hi))
}

false.fisher <- apply(false.fisher[,2:3], 1, function(x) fp.rate(x))
colnames(false.fisher) <- c('0-Hi', '1-Med', '2-Consistent')
rownames(false.fisher) <- c('fp', 'lo', 'hi')
df1 <- t(false.fisher) %>% as.data.frame() %>% mutate(test = "fisher.test", homoplasy = colnames(false.fisher))

false.pglmm <- apply(false.pglmm[,2:3], 1, function(x) fp.rate(x))
colnames(false.pglmm) <- c('0-Hi', '1-Med', '2-Consistent')
rownames(false.pglmm) <- c('fp', 'lo', 'hi')
df2 <- t(false.pglmm) %>% as.data.frame() %>% mutate(test = "binaryPGLMM", homoplasy = colnames(false.pglmm))

fp.df <- bind_rows(df1, df2)

fp.df %>% ggplot(aes(x = homoplasy, y = fp, color = test, group = test)) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2) +
  geom_line() +
  geom_point(size = 3) +
  ylab("false positive rate") +
  scale_color_manual(values = c("red", "black")) +
  theme_bw() +
  theme(legend.position = "bottom")

# why some are corrected by pglmm and other are not?
y.sim %>% filter(log.p > 4) %>% 
  ggplot(aes(y = log.p.lmm, x = ci, color = as.factor(causal))) + 
  geom_jitter(shape = 1) +
  geom_hline(yintercept = 4, linetype = 4) +
  facet_wrap(~test) +
  theme_bw()

y.sim %>% 
  ggplot(aes(x = estimate, y = log.p, color = as.factor(causal))) +   
  geom_point(size = 2, alpha = 0.5, shape = 1) +   
  scale_x_log10() + 
  geom_vline(xintercept = 1, linetype = 2) +
#  geom_hline(yintercept = c(3,7), linetype = 2) +
  geom_point(aes(x = estimate, y = log.p.lmm, color = as.factor(causal)), alpha = 0.5, size = 2) +
#  geom_text_repel(data = y, aes(x = estimate, y = log.p.lmm, label = pos), color = "red", alpha = 0.5) +
  theme_bw() +  
  xlab("odds ratio (log10)") +  
  ylab("signficance (-log10[p])") + 
  geom_hline(yintercept = 3, linetype = 2) +
  theme(legend.position = "none")

# by causal homoplasy
y.sim %>%
  ggplot(aes(log.p, log.p.lmm, color = as.factor(causal))) +
  geom_point(size = 3, alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  theme_bw() + 
  scale_color_manual(values = c("gray", "brown")) +
  geom_hline(yintercept = 4, linetype = 2, color = "red") +
  facet_wrap(~hmp.causal) +
  theme(legend.position = "bottom") 

y.sim %>% group_by(hmp.causal, causal) %>% filter(log.p >= 4 & log.p.lmm >= 4) %>% count()
y.sim %>% group_by(hmp.causal, causal) %>% filter(log.p >= 4 & log.p.lmm < 4) %>% count()

library(pROC)
y.sim <- y.sim %>% mutate(success.correction = if_else(log.p <= 4 | causal == '1', NA, if_else(log.p.lmm <= 4, '1', '0')))

y.sim %>% roc(causal, log.p, plot = T, smooth = T, percent = T, partial.auc = c(100,90))
y.sim %>% roc(causal, log.p.lmm, plot = T, smooth = T, percent = T, add = T, partial.auc = c(100,90))

# bootstrap test
roc.test(response = y.sim$causal, predictor1 = y.sim$log.p, predictor2 = y.sim$log.p.lmm, partial.auc = c(100,90), percent = T)

# group plots
roc.fisher1 <- y.sim %>% filter(hmp.causal == '0-Hi-hmp') %>% 
  roc(causal, log.p) %>% 
  coords(transpose = F) %>% 
  filter(sensitivity > 0.9, 
         specificity > 0.9) %>% 
  mutate(test = 'fisher', homoplasy = '0-Hi')

roc.lmm1 <- y.sim %>% filter(hmp.causal == '0-Hi-hmp') %>% 
  roc(causal, log.p.lmm) %>% 
  coords(transpose = F) %>% 
  filter(sensitivity > 0.9, 
         specificity > 0.9) %>% 
  mutate(test = 'pglmm', homoplasy = '0-Hi')

roc.fisher2 <- y.sim %>% filter(hmp.causal == '1-Med-hmp') %>% 
  roc(causal, log.p) %>% 
  coords(transpose = F) %>% 
  filter(sensitivity > 0.9, 
         specificity > 0.9) %>% 
  mutate(test = 'fisher', homoplasy = '1-Med')

roc.lmm2 <- y.sim %>% filter(hmp.causal == '1-Med-hmp') %>% 
  roc(causal, log.p.lmm) %>% 
  coords(transpose = F) %>% 
  filter(sensitivity > 0.9, 
         specificity > 0.9) %>% 
  mutate(test = 'pglmm', homoplasy = '1-Med')

roc.fisher3 <- y.sim %>% filter(hmp.causal == '2-Consistent') %>% 
  roc(causal, log.p) %>% 
  coords(transpose = F) %>% 
  filter(sensitivity > 0.9, 
         specificity > 0.9) %>% 
  mutate(test = 'fisher', homoplasy = '2-Consistent')

roc.lmm3 <- y.sim %>% filter(hmp.causal == '2-Consistent') %>% 
  roc(causal, log.p.lmm) %>% 
  coords(transpose = F) %>% 
  filter(sensitivity > 0.9, 
         specificity > 0.9) %>% 
  mutate(test = 'pglmm', homoplasy = '2-Consistent')


x.roc <- bind_rows(roc.fisher1, roc.lmm1, roc.fisher2, roc.lmm2,roc.fisher3, roc.lmm3)

x.roc %>% 
  ggplot(aes(threshold, specificity, color = homoplasy, group = homoplasy)) +
#  geom_point(shape =1) +
  geom_line(linewidth = 1.5) +
  facet_wrap(~test) +
  theme_bw() +
  theme(legend.position = "bottom") 

####################################

y.sim %>% filter(hmp.causal == '0-Hi-hmp') %>% 
  roc(causal, log.p, plot = T, percent = T, smooth = T) 

y.sim %>% filter(hmp.causal == '1-Med-hmp') %>% 
  roc(causal, log.p, plot = T, percent = T, smooth = T, add = T, color = 'red') 

y.sim %>% filter(hmp.causal == '2-Consistent') %>% 
  roc(causal, log.p, plot = T, percent = T, smooth = T, add = T) 

y.sim %>% 
  roc(success.correction, ci.causal, plot = T, percent = T, smooth = T) 

y.sim %>% 
  roc(success.correction, ci, plot = T, percent = T, smooth = T, add = T) 

#roc(causal ~ log.p, dat = y.sim %>% filter(hmp.causal == '2-Consistent'), plot = T, percent = T)

###############################
# Discard below
# simulate a causal SNPs with r=0.2, r=0.5, and r=0.9
# ref: https://stats.stackexchange.com/questions/12857/generate-random-correlated-data-between-a-binary-and-a-continuous-variable
# output: 
# method 1. matrix multiplication/Copulas, only 1 SNP?
sim.causal.snps.r <- function(r, size = 100, num = 10) {
  # r: desired correlation coefficient
  # size: number of OTUs; default 100
  # num: number of causal SNPs; 
  sigma <- matrix(c(1,r,r,1), ncol=2) # var-covariance matrix
  s <- chol(sigma) # choleski decomposition
  #n <- size # number of random deviates (data points)
  z <- s %*% matrix(rnorm(size * 2), nrow = 2) # 100 correlated normally distributed deviates with cor(x,y)=r
  u <- pnorm(z) # get probabilities for each deviates, 2 rows of correlated rnorms
  snp.states <- qbinom(u[1,], 1, 0.5) # discretize the 1st vector of probabilities into 0/1 with Bernoulli trial, as SNPs
  return(snp.states)
}

# method 2. logistic regression
# π(x)=exp(β0+β1x)/{1+exp(β0+β1x)}
sim.causal.snps <- function(beta, size = 100, num = 10) {
  # beta: effect size of logistic regression
  # size: OTU size
  # num: num of SNPs
  #n <- 10
  #beta0 <- -1.6
  #beta <- 0.03
  x <- runif(n = size) # a dependent variable
  pi_x <- exp(beta * x) / (1 + exp(beta * x)) # logistic prob
  y <- rbinom(n = length(x), size = 1, prob = pi_x) # binary outcomes
  data <- data.frame(x, pi_x, y)
  names(data) <- c("age", "pi", "y")
}

# use "bindata"
# https://stackoverflow.com/questions/16089178/how-to-simulate-correlated-binary-data-with-r

# use SimCorNutRes
# https://cran.r-project.org/web/packages/SimCorMultRes/vignettes/SimCorMultRes.html
# Example 3.6 (Simulation of clustered binary responses under a conditional marginal logit model without utilizing the NORTA method)
# Pr(Yit=1|xit)=F(0.2xi) where F is the cdf of the standard logistic distribution (mean =0 and variance =π2/3),
library(SimCorMultRes)
set.seed(123)
# sample size
sample_size <- 100
# cluster size
cluster_size <- 10
# intercept
beta_intercepts <- 0
# pseudo-covariate
# x <- rep(0, each = cluster_size * sample_size)
# regression parameter associated with the covariate
beta_coefficients <- 0.2
# time-stationary covariate
x <- rep(rnorm(sample_size), each = cluster_size)# simulation of clustered binary responses
# correlation matrix for the NORTA method
latent_correlation_matrix <- toeplitz(c(1, rep(0.9, cluster_size - 1)))
simulated_binary_dataset <- rbin(clsize = cluster_size, 
                                 intercepts = beta_intercepts,
                                 betas = beta_coefficients, 
                                 xformula = ~x, 
                                 cor.matrix = latent_correlation_matrix,
                                 link = "probit")
# simulated marginal probabilities
colMeans(simulated_binary_dataset$Ysim)
library("gee")
# fitting a GEE model
binary_gee_model <- gee(y ~ x, family = binomial("probit"), id = id, data = simulated_binary_dataset$simdata %>% filter(time == 3))

y.list <- split(simulated_binary_dataset$simdata, simulated_binary_dataset$simdata$time)

