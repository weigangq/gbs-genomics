# gbs snp analysis 

library(ggplot2)
library(tidyverse)
library(dplyr)
library(reshape2)

setwd('/Users/brandonely/Desktop/gbs-genomics/data')

##### snp frequency across 10K windows of genome #####

#read data
pi.windows <- read.table('snps2_10k.windowed.pi', header = T)

#plot
pi.windows %>%
  ggplot(aes(x = BIN_START, y = N_VARIANTS)) +
  geom_line(color = 'red3') +
  #formating
  theme_classic() +
  theme(title = element_text(size = 12), axis.text=element_text(size=12), axis.title=element_text(size=12), panel.border = element_rect(color = "black",fill = NA, linewidth = 1)) +
  labs(title = 'SNP Frequency', subtitle = '10K sliding windows',x = 'Genome Position', y = '# of SNPs')  


##### ts/tv counts #####

#read data and prepare dataframe
tstv <- read.table('snps2.tstv.summary', header = T)
tstv <- tstv[1:6, ]
tstv <- tstv %>% mutate(type = if_else(MODEL == 'AG' | MODEL == 'CT', 'Transition', 'Transversion'))

#order the bars for barplot
tstv$MODEL <- factor(tstv$MODEL, levels = c("AG", "CT", "AT", "AC", "GT", "CG"))

#plot
tstv %>%
  ggplot(aes(x = MODEL, y = COUNT, fill = type)) +
  geom_bar(stat = 'identity') +
  #formating
  scale_y_continuous(limits = c(0, 18000), expand = c(0, 0)) +
  theme_classic() +
  theme(title = element_text(size = 12), axis.text=element_text(size=12), axis.title=element_text(size=12), panel.border = element_rect(color = "black",fill = NA, linewidth = 1)) +
  labs(title = 'SNP Transition/Transversion',x = '', y = 'Count', fill = '')  


##### nucleotide diversity #####

#plot histogram
hist(pi.windows$PI, br = 20)

#plot line graph
pi.windows %>%
  ggplot(aes(x = BIN_START, y = PI)) +
  geom_line(color = 'steelblue') +
  #formating
  theme_classic() +
  theme(title = element_text(size = 12), axis.text=element_text(size=12), axis.title=element_text(size=12), panel.border = element_rect(color = "black",fill = NA, linewidth = 1)) +
  labs(title = 'SNP Nucleotide Diversity', subtitle = '10K sliding windows',x = 'Genome Position', y = 'Diveristy (Pi)')


##### plots for tajima's D #####

#read data
taj <- read.table('snps2_10k.Tajima.D', header = T)

#plot histogram
hist(taj$TajimaD,br=20)

#plot line graph
taj %>%
  ggplot(aes(x = BIN_START, y = TajimaD)) +
  geom_line(color = 'green4') +
  #formating
  theme_classic() +
  theme(title = element_text(size = 12), axis.text=element_text(size=12), axis.title=element_text(size=12), panel.border = element_rect(color = "black",fill = NA, linewidth = 1)) +
  labs(title = 'Tajima D', subtitle = '10K sliding windows',x = 'Genome Position', y = 'TajimaD')


##### fst analysis #####

#read data
cut0.win <- read.table('cut_0.windowed.weir.fst', header = T) %>%
  mutate('population' = 'cut_0')
cut1.win <- read.table('cut_1.windowed.weir.fst', header = T) %>%
  mutate('population' = 'cut_1')
cut2.win <- read.table('cut_2.windowed.weir.fst', header = T) %>%
  mutate('population' = 'cut_2')
cut3.win <- read.table('cut_3.windowed.weir.fst', header = T) %>%
  mutate('population' = 'cut_3')

#combine dataframes 
fst <- do.call("rbind", list(cut0.win, cut1.win, cut2.win, cut3.win))
fst <- fst %>% select(2,5,6,7)

#plot
fst %>%
  ggplot(aes(x = BIN_START, y = WEIGHTED_FST, color = population)) +
  geom_hline(yintercept = 0) +
  geom_line() + 
  facet_grid(rows = vars(population))


##### pi, tajima D, fst on same plot #####

#merge data
pi2 <- pi.windows %>% select(2,5)
taj2 <- taj %>% select(2,4)
taj2$BIN_START <- taj2$BIN_START + 1
cut1.fst2 <- cut1.win %>% select(2,5)

df_list <- list(pi2, taj2, cut1.fst2)
df <- df_list %>% reduce(full_join, by='BIN_START')
df2 <- melt(df, id.var="BIN_START")

#plot
ggplot(df2, aes(x = BIN_START, y = value)) + geom_line(aes(color = variable)) +
  labs(title = 'GBS snp analysis', subtitle = '10k sliding windows' ,x = 'genome position') +
  facet_grid(variable ~ ., scales = "free_y") + theme(legend.position = "none")

##### within group pi analysis #####

#read data
cut0_pi.win <- read.table('cut_0.windowed.pi', header = T) %>%
  mutate('population' = 'cut_0')
cut1_pi.win <- read.table('cut_1.windowed.pi', header = T) %>%
  mutate('population' = 'cut_1')
cut2_pi.win <- read.table('cut_2.windowed.pi', header = T) %>%
  mutate('population' = 'cut_2')
cut3_pi.win <- read.table('cut_3.windowed.pi', header = T) %>%
  mutate('population' = 'cut_3')

#combine data
pop_pi.win <- do.call("rbind", list(cut0_pi.win, cut1_pi.win, cut2_pi.win, cut3_pi.win))
pop_pi.win <- pop_pi.win %>% select(2,5,6)

#plot
pop_pi.win %>%
  ggplot(aes(x = BIN_START, y = PI, color = population)) +
  geom_line() + 
  facet_grid(rows = vars(population))

