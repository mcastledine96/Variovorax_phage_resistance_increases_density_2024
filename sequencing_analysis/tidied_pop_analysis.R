library(vcfR)
library(ggplot2)
library(tidyverse)
library(vegan)

mut <- read.csv('wrangled_mutations.csv',header=T)

mut$Sample[270] <- 'V3'

mut$freq<- gsub('[%]','',mut$freq)
mut$freq <- as.numeric(mut$freq)

nrow(mut[mut$freq < 50,]) #1710

nophage <- filter(mut, phage == 'no')
phage <- filter(mut, phage == "yes")

summary(as.factor(mut$Sample)) 
(128 + 142 + 141 + 171 + 155 + 147 + 143 + 171 + 139 + 139 + 131 + 160)/12
(128 + 142 + 141 + 171 + 155 + 147)/6
(143 + 171 + 139 + 139 + 131 + 160)/6
#position - number of variants

length(unique(mut$position)) #unique variants

phage_pos <- phage$position
nophage_pos <- nophage$position

nophage_unique <- filter(nophage, ! position %in% phage_pos)
phage_unique <- filter(phage, ! position %in% nophage_pos)
length(unique(nophage_unique$position))
length(unique(phage_unique$position))

470-114-117 #239 shared variants

pos_phage_summary <- group_by(phage, position) %>%
  summarise(prevalence = n(),
            freq = mean(freq),
            .groups= 'drop') %>%
  mutate(unique_position = ifelse(position %in% phage_unique$position, 'phage_only', 'both'))

library(viridis)  

phage_unique$gene <- gsub('[ABPPCLCD_]','',phage_unique$gene)
phage_unique$gene <- gsub('[←]','',phage_unique$gene)
phage_unique$gene <- gsub('[→]','',phage_unique$gene)
phage_unique$gene <- gsub('[ ]','',phage_unique$gene)

phage_unique$gene <- gsub('[ABPPCLCD_]','',phage_unique$gene)
phage_unique$gene <- gsub('[←]','',phage_unique$gene)
phage_unique$gene <- gsub('[→]','',phage_unique$gene)
phage_unique$gene <- gsub('[ ]','',phage_unique$gene)

ggplot(phage_unique, aes(x = Sample, y = gene, fill = freq)) +
  geom_tile(width = 0.75) +
  scale_fill_viridis(option = "D") +
  theme_bw() + 
  ylab("Gene") +
  scale_x_discrete(labels = c("1", "2", "3", "4", "5", "6")) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 12, hjust = 0), axis.text = element_text(size = 6, colour = 'black'), axis.title = element_text(size = 12, colour = 'black'), 
        panel.grid.minor = element_blank()) +
  labs(fill = "Frequency") +
  geom_vline(aes(xintercept = 1.5), size = 0.5) +
  geom_vline(aes(xintercept = 2.5), size = 0.5) +
  geom_vline(aes(xintercept = 3.5), size = 0.5) +
  geom_vline(aes(xintercept = 4.5), size = 0.5) +
  geom_vline(aes(xintercept = 5.5), size = 0.5)

### This looks at different genes - reread in data and then go from here. Gathers some descriptive data on genes not found to be mutated in each population (ignores variant level differences which are examined above)

mut <- read.csv('wrangled_mutations.csv',header=T)

mut$Sample[270] <- 'V3'

mut$freq<- gsub('[%]','',mut$freq)
mut$freq <- as.numeric(mut$freq)

nophage <- filter(mut, phage == 'no')
phage <- filter(mut, phage == "yes")

nophage_2 <- group_by(nophage, gene, description, phage) %>%
  summarise(m_freq = mean(freq),
            se = sd(freq)/sqrt(length(freq)),
            n_mutations = length(freq))

phage_2 <- group_by(phage, gene, description, phage) %>%
  summarise(m_freq = mean(freq),
            se = sd(freq)/sqrt(length(freq)),
            n_mutations = length(freq))

#need to filter out ones seen in nophage from phage

no_phage_genes <- nophage_2[,c(1,2)]
phage_genes <- phage_2[,c(1,2)]

diff_genes <- setdiff(no_phage_genes, phage_genes)
same_genes <- semi_join(no_phage_genes, phage_genes)

phage_unique <- phage_2[!(phage_2$gene %in% same_genes$gene),]
nophage_unique <- nophage_2[!(nophage_2$gene %in% same_genes$gene),]

#NMDS plot

is.numeric(mut$freq)

mut$freq2 <- round(mut$freq, 0)

d_clust <- select(mut, Sample, position, gene, phage, freq2) %>%
  unite(., 'id', c(Sample, phage), sep ='_')  %>%
  spread(., id, freq2)%>%
  arrange(position) %>%
  group_by(gene) %>%
  ungroup() 

d_clust[is.na(d_clust)] <- 0

d_clust$name <- interaction(d_clust$position, d_clust$gene)

d_clust <- d_clust[,-c(1,2)]
d_clust <- d_clust[,c(13,1:12)]
nams <- d_clust[,1]
d_clust <- d_clust[,-1]
row.names(d_clust) <- nams$name
d_clust2 <- t(d_clust)

d_clust2 <- as.matrix(d_clust2)

nmds_1 <- metaMDS(d_clust2, distance = 'euclidean')
stressplot(nmds_1)

plot(nmds_1)
treat=c(rep("No Phage",6),rep("Phage",6))
ordihull(nmds_1,groups=treat,draw="polygon",col="grey90",label=T)

data.scores <- as.data.frame(scores(nmds_1)$sites)
data.scores$site <- rownames(data.scores)  
data.scores$grp<-treat

adon.results<-adonis2(d_clust2 ~ treat)
print(adon.results)

dist_from_00 <- function(x, y){
  return(sqrt((0 - x)^2+(0-y)^2))
}
library(ggvegan)
d_nmds <- fortify(nmds_1) %>%
  janitor::clean_names() %>%
  mutate_if(., is.factor, as.character)


# wrangle sites
d_sample <- filter(d_nmds, score == 'sites') %>%
  separate(., label, c('rep', 'phage'), sep = '_')  %>%
  mutate(dist = dist_from_00(nmds1, nmds2))

d_gene <- filter(d_nmds, score == 'species') %>%
  rename(., gene_name = label) %>%
  mutate(dist = dist_from_00(nmds1, nmds2))

d_sample$label <- c("1", "2", "3", "4", "5", "6", "1", "2", "3", "4", "5", "6")

ggplot() +
  geom_segment(aes(x = 0, y = 0, yend = nmds2, xend = nmds1, group = score), d_gene, arrow = arrow(length = unit(0.01, "npc")), col = 'lightgrey') +
  geom_point(aes(nmds1, nmds2, col = phage), size = 5, data = filter(d_sample)) +
  theme_bw() +
  palettetown::scale_color_poke(pokemon = 'blastoise', spread = 4, labels = c("Susceptible", "Resistant")) +
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  theme(axis.text = element_text(size = 11, colour = 'black'), axis.title = element_text(size = 12, colour = 'black'), strip.background = element_blank(), strip.text = element_text(size = 11, hjust = 0), legend.text = element_text(size = 12, colour = 'black'), legend.title = element_text(size = 12, colour = 'black'),legend.position = "bottom") +
  labs(col = "Phage") +
  ggforce::geom_mark_hull(aes(nmds1, nmds2, group = phage), d_sample) 
mut3 <- group_by(mut, phage, Sample) %>%
  summarise(., dist = sum(freq))

m1 <- lm(dist ~ phage, data = mut3)
m2 <- lm(dist ~ 1, data = mut3)
anova(m1,m2)
