library(readr)
library(qqman)
library(tidyverse)

# read in gemma linear mixed model results
results <- read_table("../gbs_lmm_clean.txt")

# compute bonferroni threshold, can be added to the manhattan plot
bonferroni <- log10(0.05/nrow(results))

# manhattan plot
manhattan(results, chr='chr', bp = 'ps', p="af", snp= "rs", annotateTop = TRUE)

