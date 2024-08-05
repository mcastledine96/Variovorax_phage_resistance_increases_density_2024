#load packages
library(ggplot2)
library(dplyr)
library(tidyverse)
library(lme4)
library(patchwork)
library(emmeans)

#data
b <- read.csv('data/vario_r_comm_010724.csv',header=T)

b <- mutate(b,
            prop = count / total)

m1 <- glmer(cbind(count, total-count) ~ treat * bact + (1|com), data = b, family = binomial)
m2 <- glmer(cbind(count, total-count) ~ treat + bact + (1|com), data = b, family = binomial)
anova(m1,m2)
emmeans::emmeans(m1, pairwise ~ bact|treat)
emmeans::emmeans(m1, pairwise ~ treat|bact, type = 'response') #p does worst with phage resistant V
emmeans::emmeans(m1, pairwise ~ treat|bact)

rmen <- emmeans::emmeans(m1, pairwise ~ treat|bact)$contrasts
rmen <- data.frame(rmen)

rmen <- rmen[,c(2,1,3,6,7)]
rmen$bact <- gsub('o', 'O', rmen$bact)
rmen$bact <- gsub('p', 'P', rmen$bact)
rmen$bact <- gsub('v', 'V', rmen$bact)
rmen$contrast <- gsub('anc', 'Ancestral', rmen$contrast)
rmen$contrast <- gsub(' s', ' Susceptible', rmen$contrast)
rmen$contrast <- gsub('r ', 'Resistant ', rmen$contrast)
rmen$contrast <- gsub(' r', ' Resistant ', rmen$contrast)

rmen <- mutate(rmen,
              estimate = round(estimate, digits = 3),
              z.ratio = round(z.ratio, digits = 3),
              p.value = round(p.value, digits = 3))

rmen$p.value[rmen$p.value == 0] <- '<0.001'

rmen_tab <- flextable(rmen) %>%
  set_header_labels(bact = "Species", contrast = "Contrast", time = 'Week', estimate = 'Estimate', z.ratio = "z-ratio", p.value = "p-value") %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 12, part = 'all') %>%
  autofit() %>%
  align(j = c(1:5), align = 'center', part = 'all') %>%
  align(j = 1, align = 'center', part = 'header') %>%
  hline(i = c(3,6,9), border = officer::fp_border(color="black")) 

cols <- c("white","lightgrey", "black", "white","lightgrey", "black")

ggplot(b, aes(x = bact, y = prop, group = interaction(bact, treat))) +
  geom_boxplot(aes(fill = treat),col = 'darkgrey') +
  geom_point(position = position_dodge(0.75), shape = 21, fill = 'white', size =2, col = "black") +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black'), axis.title = element_text(size = 12, colour = 'black'), legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 12), strip.background = element_blank(), strip.text = element_text(hjust = 0, size = 12)) +
  scale_fill_manual(values = cols, labels = c('Ancestral', 'Phage\nresistant', 'Phage\nsusceptible')) +
  labs(fill = expression(italic("Variovorax")~`phage resistance`)) +
  ylab("Relative proportion") +
  xlab('Species') +
  scale_x_discrete(labels = c(expression(italic("Ochrobactrum")), expression(italic("Pseudomonas")), expression(italic("Variovorax"))))
