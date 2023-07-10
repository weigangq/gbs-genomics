# gbs snp analysis 

library(ggplot2)
library(tidyverse)
library(dplyr)

setwd('/Users/brandonely/Desktop/gbs-genomics/data')

#read vcftools outputs
pi.site <- read.table('snps2.sites.pi', header = T)
pi.windows <- read.table('snps2_10k.windowed.pi', header = T)
taj <- read.table('snps2_10k.Tajima.D', header = T)
tstv <- read.table('snps2.tstv.summary', header = T)


### snp frequency across 10K windows of genome ###
pi.windows %>%
  ggplot(aes(x = BIN_START, y = N_VARIANTS)) +
  geom_bar(stat = 'identity',color = 'red3') +
  #formating
  theme_classic() +
  theme(title = element_text(size = 12), axis.text=element_text(size=12), axis.title=element_text(size=12), panel.border = element_rect(color = "black",fill = NA, linewidth = 1)) +
  labs(title = 'SNP Frequency', subtitle = '10K sliding windows',x = 'Genome Position', y = '# of SNPs')  


### ts/tv counts ###
tstv_2 <- tstv[1:6, ]
tstv_2 <- tstv_2 %>% mutate(type = if_else(MODEL == 'AG' | MODEL == 'CT', 'Transition', 'Transversion'))
tstv_2$MODEL <- factor(tstv_2$MODEL, levels = c("AG", "CT", "AT", "AC", "GT", "CG"))

tstv_2 %>%
  ggplot(aes(x = MODEL, y = COUNT, fill = type)) +
  geom_bar(stat = 'identity') +
  #formating
  scale_y_continuous(limits = c(0, 18000), expand = c(0, 0)) +
  theme_classic() +
  theme(title = element_text(size = 12), axis.text=element_text(size=12), axis.title=element_text(size=12), panel.border = element_rect(color = "black",fill = NA, linewidth = 1)) +
  labs(title = 'SNP Transition/Transversion',x = '', y = 'Count', fill = '')  


### plots for nucleotide diversity ###
hist(pi.windows$PI, br = 20)
boxplot(pi.windows$PI,ylab="diversity")

pi.windows %>%
  ggplot(aes(x = BIN_START, y = PI)) +
  geom_line(color = 'steelblue') +
  #formating
  theme_classic() +
  theme(title = element_text(size = 12), axis.text=element_text(size=12), axis.title=element_text(size=12), panel.border = element_rect(color = "black",fill = NA, linewidth = 1)) +
  labs(title = 'SNP Nucleotide Diversity', subtitle = '10K sliding windows',x = 'Genome Position', y = 'Diveristy (Pi)')


### plots for tajima's D ###
hist(taj$TajimaD,br=20)

taj %>%
  ggplot(aes(x = BIN_START, y = TajimaD)) +
  geom_line(color = 'purple4') +
  #formating
  theme_classic() +
  theme(title = element_text(size = 12), axis.text=element_text(size=12), axis.title=element_text(size=12), panel.border = element_rect(color = "black",fill = NA, linewidth = 1)) +
  labs(title = 'Tajima D', subtitle = '10K sliding windows',x = 'Genome Position', y = 'TajimaD')



