#Tidy R script containing all PVO analyses relevant to Vario script

library(ggplot2)
library(dplyr)
library(tidyverse)
library(lme4)
library(patchwork)
library(emmeans)

den_dat <- read.csv("data/pvodensity_time_121120.csv", header = T)

#calculate population densities

den_dat <- mutate(den_dat,
                  cfus = count * dil * dil_fact * froz,
                  logcfus = log10(cfus))

#take averages

den_dat2 <- group_by(den_dat, phage, sp, time, type) %>%
  summarise(., mean_logden = mean(logcfus),
            sd = sd(logcfus),
            se = sd/sqrt(length(logcfus)))

#tidy script
den_dat <- den_dat[, c(1:5, 10:12)]

#look at all species for the supplement

com_ph <- filter(den_dat, type == "comm")
com_ph2 <- filter(den_dat2, type == "comm")
mono_ph <- filter(den_dat, type == "mono")
mono_ph2 <- filter(den_dat2, type == "mono")

labels = c("2", "4", "6", "8")

bact_labs <- c("(i) O", "(ii) P", "(iii) V")
names(bact_labs) <- c("O", "P", "V")

bact_den1 <- ggplot() +
  geom_point(data = mono_ph, aes(x = time, y = logcfus, group = interaction(evo_line, rep), col = phage), alpha = 0.3) +
  facet_wrap(~sp, labeller = labeller(sp = bact_labs)) +
  geom_line(data = mono_ph, aes(x = time, y = logcfus, group = interaction(evo_line, rep), col = phage), alpha = 0.3) +
  geom_point(data = mono_ph2, aes(x = time, y = mean_logden, group = phage, col = phage), size = 3) +
  geom_line(data = mono_ph2, aes(x = time, y = mean_logden, group = phage, col = phage), size = 1) +
  geom_errorbar(data = mono_ph2, aes(x = time, ymin = mean_logden - se, ymax = mean_logden+se, col = phage), width = 0.35, size = 0.8) +
  theme_bw() +
  xlab("Week") +
  ylab("Cell density log10(CFU/mL)") +
  theme(axis.text = element_text(size = 11, colour = 'black'), axis.title = element_text(size = 12, colour = 'black'), strip.background = element_blank(), strip.text = element_text(size = 11, hjust = 0), legend.text = element_text(size = 12), legend.position = "none", legend.title = element_text(size = 12)) +
  labs(title = "(a) Monoculture", col = "Phage")+
  palettetown::scale_color_poke(pokemon = 'blastoise', spread = 4, labels = c("Absent", "Present")) +
  scale_x_discrete(labels = labels)

bact_den2 <- ggplot() +
  geom_point(data = com_ph, aes(x = time, y = logcfus, group = interaction(evo_line, rep), col = phage), alpha = 0.3) +
  facet_wrap(~sp, labeller = labeller(sp = bact_labs)) +
  geom_line(data = com_ph, aes(x = time, y = logcfus, group = interaction(evo_line, rep), col = phage), alpha = 0.3) +
  geom_point(data = com_ph2, aes(x = time, y = mean_logden, group = phage, col = phage), size = 3) +
  geom_line(data = com_ph2, aes(x = time, y = mean_logden, group = phage, col = phage), size = 1) +
  geom_errorbar(data = com_ph2, aes(x = time, ymin = mean_logden - se, ymax = mean_logden+se, col = phage), width = 0.35, size = 0.8) +
  theme_bw() +
  xlab("Week") +
  ylab("Cell density log10(CFU/mL)") +
  theme(axis.text = element_text(size = 11, colour = 'black'), axis.title = element_text(size = 12, colour = 'black'), strip.background = element_blank(), strip.text = element_text(size = 11, hjust = 0), legend.position = "bottom") +
  labs(title = "(b) Polyculture", col = "Phage")+
  palettetown::scale_color_poke(pokemon = 'blastoise', spread = 4) +
  scale_x_discrete(labels = labels)

bact_den1 + bact_den2 + plot_layout(nrow = 2)

#grab out vario data

v_den <- filter(den_dat, sp == "V")
v_den2 <- filter(den_dat2, sp == "V")

v_den$type[v_den$type == "comm"] <- "(b) Polyculture"
v_den$type[v_den$type == "mono"] <- "(a) Monoculture"
v_den2$type[v_den2$type == "comm"] <- "(b) Polyculture"
v_den2$type[v_den2$type == "mono"] <- "(a) Monoculture"

ggplot() +
  geom_point(data = v_den, aes(x = time, y = logcfus, group = interaction(evo_line, rep), col = phage), alpha = 0.3) +
  facet_wrap(~type) +
  geom_line(data = v_den, aes(x = time, y = logcfus, group = interaction(evo_line, rep), col = phage), alpha = 0.3) +
  geom_point(data = v_den2, aes(x = time, y = mean_logden, group = phage, col = phage), size = 3) +
  geom_line(data = v_den2, aes(x = time, y = mean_logden, group = phage, col = phage), size = 1) +
  geom_errorbar(data = v_den2, aes(x = time, ymin = mean_logden - se, ymax = mean_logden+se, col = phage), width = 0.2, size = 0.8) +
  theme_bw() +
  xlab("Week") +
  ylab("Density log10(CFU/mL)") +
  theme(axis.text = element_text(size = 11, colour = 'black'), axis.title = element_text(size = 12, colour = 'black'), strip.background = element_blank(), strip.text = element_text(size = 12, hjust = 0), legend.text = element_text(size = 12), legend.position = "bottom", legend.title = element_text(size = 12)) +
  palettetown::scale_color_poke(pokemon = 'blastoise', spread = 4, labels = c("Absent", "Present")) +
  scale_x_discrete(labels = labels) + 
  labs(col = "Phage")

#phage densities
phage_den <- read.csv("data/phage_densities_101220.csv", header = T)

#calculate densities
phage_den <- mutate(phage_den,
                    pfus = count * dilution * dil_fact,
                    logpfus = log10(pfus))

#sort out infinite values in logged data
phage_den$logpfus[phage_den$logpfus == -Inf] <- 0

#mean of phage densities

phage_den2 <- na.omit(phage_den)
phage_means <- group_by(phage_den2, type, species, time) %>%
  summarise(mean_den = mean(logpfus),
            se = sd(logpfus)/sqrt(length(logpfus)))

#plot for individual evolution lines
phage_means$time <- substr(phage_means$time, 2, 2)
phage_den$time <- substr(phage_den$time, 2, 2)

phage_lab <- c("(a) ORM_20", "(b) CHF7MC", "(c) VAC_51")
names(phage_lab) <- c("O", "P", "V")

phageplot <- ggplot() +
  geom_point(data = phage_den, aes(x = time, y = logpfus, group = line, shape = type, group = type), alpha = 0.3) +
  facet_wrap(~species, labeller = labeller(species = phage_lab)) +
  geom_line(data = phage_den, aes(x = time, y = logpfus, group = line, shape = type, linetype = type), alpha = 0.3) +
  palettetown::scale_color_poke(pokemon = 'bulbasaur', spread = 4, labels = c("Polyculture", "Monoculture")) +
  theme_bw() +
  ylab("Phage density log10(pfu/mL)") +
  xlab("Week") +
  theme(axis.text = element_text(size = 11, colour = 'black'), axis.title = element_text(size = 12, colour = 'black'), strip.background = element_blank(), strip.text = element_text(size = 11, hjust = 0), legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
  geom_point(data = phage_means, aes(x = time, y = mean_den, group = type, shape = type, linetype = type), size = 2.5, position = position_dodge(0.5)) +
  geom_errorbar(data = phage_means, aes(x = time, ymin = mean_den - se, ymax = mean_den + se, group = type), position = position_dodge(0.5), width = 0.5, size = 0.8) +
  geom_line(data = phage_means, aes(x = time, y = mean_den, group = type, shape = type, linetype = type), position = position_dodge(0.5), size = 0.8) +
  labs(shape = "Treatment", linetype = "Treatment") +
  scale_shape_discrete(labels = c("Polyculture", "Monoculture")) +
  scale_linetype_discrete(labels = c("Polyculture", "Monoculture"))

v_phage <- filter(phage_den, species == "V")
v_phage <- filter(v_phage, time == "2")

#phage density analyses
p1 <- lm(logpfus ~ type, data = v_phage)
p2 <- lm(logpfus ~ 1, data = v_phage)
anova(p1,p2,test="F")
v_phage2 <- filter(phage_means, species == "V")

#bacterial density analyses
v_den$rando <- interaction(v_den$evo_line, v_den$rep)

v_mod1 <- lmer(logcfus ~ time * type * phage + (1|rando), data = v_den)
v_mod1.1 <- lmer(logcfus ~ time + type + phage + time:type + time:phage + phage:type + (1|rando), data = v_den)
anova(v_mod1, v_mod1.1) #non-sig 3 way interaction
v_mod2 <- lmer(logcfus ~ time + type + phage + time:phage + type:phage + (1|rando), data = v_den) #removing time:type interaction 
anova(v_mod1, v_mod2) #sig time:type interaction
v_mod3 <- lmer(logcfus ~ time + type + phage + time:type + type:phage + (1|rando), data = v_den) #removing time:phage interaction
anova(v_mod1, v_mod3) #non-sig time:phage interaction
v_mod4 <- lmer(logcfus ~ time + type + phage + time:phage + time:type + (1|rando), data = v_den) #remove type:phage interaction
anova(v_mod1, v_mod4) #non-sig type:phage interaction. This one is less sig so drop first
v_mod5 <- lmer(logcfus ~ time + type + phage + time:type + (1|rando), data = v_den) #drop time:phage
anova(v_mod4, v_mod5) #still time:phage non-sig
v_mod6 <- lmer(logcfus ~ time + type + phage + time:phage + (1|rando), data = v_den) #drop time:type
anova(v_mod4, v_mod6) #time:type significant
v_mod6.1 <- lmer(logcfus ~ time + type + phage + (1|rando), data = v_den) #drop time:type
anova(v_mod5, v_mod6.1)
v_mod7 <- lmer(logcfus ~ time + type + time:type + (1|rando), data = v_den)
anova(v_mod5, v_mod7) 

emmeans::emmeans(v_mod5, pairwise ~ type|time) #difference in poly to mono decreases over time
emmeans::emmeans(v_mod5, pairwise ~ phage) #presence of phage overall increases V density

library(flextable)
mens <- data.frame(emmeans::emmeans(v_mod5, pairwise ~ type|time)$contrasts)
mens <- mens[,c(1:3,6,7)]

mens$contrast <- "Poly. - Mono."

mens <- mutate(mens,
               estimate = round(estimate, digits = 3),
               t.ratio = round(t.ratio, digits = 3),
               p.value = round(p.value, digits = 3))

mens$p.value[mens$p.value == 0] <- "<0.001"
mens$time <- substring(mens$time, 2)

men_flex <- flextable(mens) %>%
  set_header_labels(time = "Week", contrast = "Contrast", estimate = "Estimate", t.ratio = "t-ratio", p.value = "p-value") %>%
  font(fontname = 'Times', part = 'all') %>%
  fontsize(size = 16, part = 'all') %>%
  autofit() 

#clonal diversity

vclones <- read.csv("data/v_clonal_diversity_101121.csv", header = T)

vclones <- mutate(vclones,
                  den = count * dil * dil_fact * 2,
                  logden = log10(den))

ggplot(vclones, aes(x = time, y = logden, group = interaction(time, evo_line), col = evo_line, fill = evo_line)) +
  MicrobioUoE::geom_pretty_boxplot() +
  geom_point(shape = 21, fill = "white", col = "black", position = position_dodge(0.75), size = 2) +
  theme_bw() +
  palettetown::scale_color_poke(pokemon = "blastoise", spread = 4, labels = c("Phage Absent", "Phage Present")) +
  palettetown::scale_fill_poke(pokemon = "blastoise", spread = 4, labels = c("Phage Absent", "Phage Present")) +
  theme_bw() +
  theme(axis.text = element_text(size = 11, colour = 'black'), axis.title =  element_text(size = 12, colour = 'black'), legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
  labs(fill = "Evolutionary history", col = "Evolutionary history") +
  ylab("Density log10(CFU/mL)") +
  xlab("Week") +
  scale_x_discrete(labels = c("2", "4"))

#analyses
vclones$rep2 <- interaction(vclones$evo_line, vclones$rep)

mod1 <- lmer(logden ~ evo_line * time + (1|rep2) + (1|block), data = vclones)
plot(mod1)
mod2 <- lmer(logden ~ evo_line + time + (1|rep2) + (1|block), data = vclones)
anova(mod1, mod2)
mod3 <- lmer(logden ~ time + (1|rep2) + (1|block), data = vclones) #sig effect of evo_line
anova(mod2, mod3) 
mod4 <- lmer(logden ~ evo_line + (1|rep2) + (1|block), data = vclones)
anova(mod2, mod4) #sig effect of time

emmeans::emmeans(mod2, pairwise ~ evo_line)
emmeans::emmeans(mod2, pairwise ~ time)

#phage resistance

res <- read.csv("data/vario_rcomp_260323.csv",header=T)

res <- mutate(res,
              count = count_t - count_c,
              den = count * dil * plated * froz,
              den_c = count_c * dil * plated * froz,
              t_0 = ((t0 * 10 * 40 * 2)*20)/6,
              t0_c = ((167 * 10 * 40 * 2)*20)/6)

res <- mutate(res,
              g_rate = log(den / t_0),
              g_ratec = log(den_c / t0_c),
              m = g_rate / g_ratec)

m1 <- lm(m ~ treat, data = res)
m2 <- lm(m ~ 1, data = res)
anova(m1,m2,test="F")
plot(m1)

res$treat[res$treat == "r"] <- "VP"
res$treat[res$treat == "s"] <- "V"

ggplot(res, aes(x = treat, y = m)) +
  MicrobioUoE::geom_pretty_boxplot(aes(fill = treat, col = treat)) +
  geom_point(position = position_jitter(0.05), shape = 21, fill = "white", size = 3)+
  xlab("Clone") +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black'), axis.title = element_text(size = 14, colour = 'black'), strip.background = element_blank(), legend.position = "none") +
  scale_color_manual(values = c("black", "#4080c0", "#684739")) +
  scale_fill_manual(values = c("black", "#4080c0", "#684739")) +
  ylab("Relative fitness") +
  scale_x_discrete(labels = c("Ancestor", "Phage susceptible", "Phage resistant")) +
  geom_hline(yintercept = 1, linetype = "dashed")

#growth curve data

dat <- read.csv("data/vario_daily_growth.csv",header=T)

dat$rep <- interaction(dat$Treat, dat$rep)

mod <- lmer(log(OD) ~ day + I(day^2) * Treat + (1|rep), data = dat)
qqnorm(resid(mod))
plot(mod)
mod2 <- lmer(log(OD) ~ day + I(day^2) + Treat + (1|rep), data = dat)
anova(mod,mod2)

library(ciTools)

newdata <- with(dat, expand.grid(Treat=unique(Treat), rep=unique(rep), day = seq(min(dat$day),max(dat$day),length.out=100)))

new <- add_ci(newdata, mod, alpha = 0.05, type = "boot", includeRanef = FALSE, nSims = 2000)

new_m <- group_by(new, Treat, day) %>%
  summarise(men = mean(pred),
            men_l = mean(LCB0.025),
            men_u = mean(UCB0.975))
us <- filter(new_m, day == 7)
us2 <- filter(new_m, day == 4.5)

new_m <- mutate(new_m, 
                m_2 = exp(men),
                m_l2 = exp(men_l),
                m_u2 = exp(men_u))

#back transform estimates
new2 <- mutate(new, 
                pred2 = exp(pred),
                LCB0.0252 = exp(LCB0.025),
                UCB0.9752 = exp(UCB0.975))


phage_lab <- c("(a) Phage susceptible", "(b) Phage resistant")
names(phage_lab) <- c("V", "VP")

ggplot() +
  geom_point(data = dat, aes(x = day, y = log(OD)), position = position_jitter(0.1), size = 2) + 
  geom_line(data = new, aes(x = day, y = pred)) +
  geom_ribbon(data = new, aes(x = day, ymin = LCB0.025, ymax = UCB0.975, fill = Treat), alpha = 0.3) +
  facet_wrap(~Treat, labeller = labeller(Treat = phage_lab)) +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour = 'black'), axis.title = element_text(size = 12, colour = 'black'), legend.position = "none", strip.background = element_blank(), legend.text = element_text(size = 12), legend.title = element_text(size = 12),  strip.text = element_text(hjust = 0, size = 12)) +
  xlab("Day") +
  scale_x_continuous(n.breaks = 7) +
  labs(y = expression(Optical~Density~log(OD["600"]))) +
  scale_fill_manual(values = c("#4080c0", "#684739"))

#vario growth over 24hrs

dat <- read.csv("data/growth_part1_vario_081222.csv",header=T)

dat <- mutate(dat,
              time_hr = time / 3600)

#1 - 6 = V (minus rep 3). 7-12 = VP
#C = rep 1, D = rep 2, E = rep 3

dat <- dat %>%
  separate(well, c("row", "col"))

dat2 <- filter(dat, row == "B" | row == "C" | row == "D" | row == "E")

dat2$treat <- as.numeric(dat2$col)

dat2$treat[dat2$row == "C" & dat2$treat == 1 | dat2$row == "C" & dat2$treat == 2 | dat2$row == "C" & dat2$treat == 3 | dat2$row == "C" & dat2$treat == 4 | dat2$row == "C" & dat2$treat == 5 | dat2$row == "C" & dat2$treat == 6] <- "V"

dat2$treat[dat2$row == "D" & dat2$treat == 1 | dat2$row == "D" & dat2$treat == 2 | dat2$row == "D" & dat2$treat == 3 | dat2$row == "D" & dat2$treat == 4 | dat2$row == "D" & dat2$treat == 5 | dat2$row == "D" & dat2$treat == 6] <- "V"

dat2$treat[dat2$row == "E" & dat2$treat == 1 | dat2$row == "E" & dat2$treat == 2 | dat2$row == "E" & dat2$treat == 3 | dat2$row == "E" & dat2$treat == 4 | dat2$row == "E" & dat2$treat == 5 | dat2$row == "E" & dat2$treat == 6] <- "V"

dat2$treat[dat2$row == "C" & dat2$treat == 7 | dat2$row == "C" & dat2$treat == 8 | dat2$row == "C" & dat2$treat == 9 | dat2$row == "C" & dat2$treat == 10 | dat2$row == "C" & dat2$treat == 11 | dat2$row == "C" & dat2$treat == 12] <- "VP"

dat2$treat[dat2$row == "D" & dat2$treat == 7 | dat2$row == "D" & dat2$treat == 8 | dat2$row == "D" & dat2$treat == 9 | dat2$row == "D" & dat2$treat == 10 | dat2$row == "D" & dat2$treat == 11 | dat2$row == "D" & dat2$treat == 12] <- "VP"

dat2$treat[dat2$row == "E" & dat2$treat == 7 | dat2$row == "E" & dat2$treat == 8 | dat2$row == "E" & dat2$treat == 9 | dat2$row == "E" & dat2$treat == 10 | dat2$row == "E" & dat2$treat == 11 | dat2$row == "E" & dat2$treat == 12] <- "VP"

dat2$treat[dat2$row == "B"] <- "control"

dat2$rep1 <- interaction(dat2$treat, dat2$col) #biological rep
dat2$rep2 <- interaction(dat2$treat, dat2$row) #technical rep

ggplot(dat2, aes(x = time_hr, y = od, col = rep1)) +
  geom_point() +
  facet_wrap(~treat)

dat2 <- filter(dat2, ! treat == "control")

ggplot(dat2, aes(x = time_hr, y = od, col = rep1)) +
  geom_point() +
  facet_wrap(~treat)

dat2 <- filter(dat2, ! col == 4) #failed to grow

dat2$od2 <- log(dat2$od)

#rolling reg

library(tidyverse) #install.packages(tidyverse)
library(zoo) #install.packages(zoo)
library(broom) #install.packages(broom)
library(growthcurver) # install.packages(growthcurver)
library(nls.multstart)

roll_regress <- function(x){
  temp <- data.frame(x)
  mod <- lm(temp)
  temp <- data.frame(slope = coef(mod)[[2]],
                     slope_lwr = confint(mod)[2, ][[1]],
                     slope_upr = confint(mod)[2, ][[2]],
                     intercept = coef(mod)[[1]],
                     rsq = summary(mod)$r.squared, stringsAsFactors = FALSE)
  return(temp)
}

num_points = ceiling(1.5*60/(60*0.167)) 

models <- dat2 %>%
  group_by(rep1, rep2, treat) %>%
  do(cbind(model = select(., od2, time_hr) %>% 
             zoo::rollapplyr(width = num_points, roll_regress, by.column = FALSE, fill = NA, align = 'center'),
           time = select(., time_hr),
           ln_od = select(., od2))) %>%
  rename_all(., gsub, pattern = 'model.', replacement = '')

preds <- models %>%
  filter(., !is.na(slope)) %>%
  group_by(time_hr) %>%
  do(data.frame(time2 = c(.$time_hr - 2, .$time_hr + 2))) %>%
  left_join(., models) %>%
  mutate(pred = (slope*time2) + intercept)

preds <- models %>%
  filter(., !is.na(slope)) %>%
  left_join(., models) 

growth_rate <- filter(models, slope == max(slope, na.rm = TRUE))

grows <- group_by(growth_rate, rep1, treat) %>%
  summarise(m.grow = mean(slope),
            slope_lwr = mean(slope_lwr),
            slope_upr = mean(slope_upr))

mod1 <- lm(m.grow ~ treat, data = grows)
mod2 <- lm(m.grow ~ 1, data = grows)
anova(mod1,mod2,test="F")
plot(mod1)

### Supernatent analysis

sup <- read.csv("data/sup_assay_1w.csv",header=T)

sup <- mutate(sup,
              cfu = count * dil * plated,
              cful = log10(cfu))

m1 <- lmer(cful ~ gen * sup + (1|ID), data = sup)
m2 <- lmer(cful ~ gen + sup + (1|ID), data = sup)
anova(m1,m2)
m3 <- lmer(cful ~ gen + (1|ID), data = sup)
anova(m2,m3)
m4 <- lmer(cful ~ sup + (1|ID), data = sup)
anova(m2,m4)
m5 <- lmer(cful ~ 1 + (1|ID), data = sup)
anova(m3,m5)

ggplot() +
  geom_line(data = sup, aes(x = sup, y = cful, group = ID)) +
  geom_point(data = sup, aes(x = sup, y = cful, group = ID), shape = 21, fill = "white", size = 2) +
  facet_wrap(~gen, labeller = labeller(gen = phage_lab)) +
  theme_bw() +
  theme(axis.text = element_text(size = 16, colour = 'black'), axis.title = element_text(size = 16, colour = 'black'), legend.position = "none", strip.background = element_blank(), legend.text = element_text(size = 16), legend.title = element_text(size = 16),  strip.text = element_text(hjust = 0, size = 16)) +
  xlab("Supernatant source population") +
  ylab("Density log10(CFU/mL)") +
  scale_x_discrete(labels = c("Susceptible", "Resistant", "Susceptible", "Resistant"))
