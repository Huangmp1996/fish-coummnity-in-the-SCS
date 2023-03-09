library(dplyr)
library(ggplot2)
library(rfishbase)

## read data
total_clean_from_lq <- read.csv('total_clean_from_lq.csv', encoding = 'UTF-8',
                                header = T, colClasses=c("X.U.FEFF.date" = "character"))
total_clean_from_lq <- mutate(total_clean_from_lq, 
                       site = recode(site,
                                     '美济礁' = 'Meiji',
                                     '永暑礁北沙洲' = 'Yongshu', 
                                     '七连屿' = 'Qilian',
                                     '永乐环礁' = 'Yongle',
                                     '南扉暗沙' = 'Zhongsha', 
                                     '漫步暗沙' = 'Zhongsha')) %>% 
  dplyr::filter(site != 'Zhongsha')

distribution <- read.csv('species_distribution.csv', header = T, encoding = 'UTF-8')
species_name <- select(distribution, X.U.FEFF.species, latin, Order, Family) %>% 
  distinct()
total <- left_join(total_clean_from_lq, species_name, by = c('species' = 'X.U.FEFF.species'))
total <- dplyr::filter(total, count > 0)

## survey duration
video_time <- group_by(total_clean_from_lq, X.U.FEFF.date, camera, video, site) %>% 
  slice_tail() %>% tidyr::separate(timestamp, into = c('start', 'end'), sep = '-')


## community compare
distribution <- mutate(distribution, 
                       site = recode(site,
                                    '美济礁' = 'Meiji',
                                    '永暑礁北沙洲' = 'Yongshu', 
                                    '七连屿' = 'Qilian',
                                    '永乐环礁' = 'Yongle',
                                    '南扉暗沙' = 'Zhongsha', 
                                    '漫步暗沙' = 'Zhongsha')) %>% 
  dplyr::filter(site != 'Zhongsha')
fish_richness <- group_by(distribution, site, Family) %>% 
    summarise(richness = n_distinct(latin)) %>% arrange(desc(richness))
fish_richness$site <- factor(fish_richness$site,
                            levels = c('Meiji', 'Yongshu', 'Qilian',
                                       'Zhongsha', 'Yongle'))
fish_richness_wide <- tidyr::spread(fish_richness, site, richness)
fish_richness_wide <- tibble::column_to_rownames(fish_richness_wide, 'Family')
fish_richness_wide[is.na(fish_richness_wide)] <- 0
# fish_richness_wide <- t(fish_richness_wide)
library(pheatmap)
pheatmap(fish_richness_wide,
         cellwidth = 20, cellheight = 10,
         cluster_cols = F, angle_col = c(45),
         treeheight_row = 0,
         legend_breaks = c(0, 10, 20, 30, 40, 43), 
         legend_labels = c("0", "10", "20", "30", "40", "species number"),
         filename = 'Result/heatmap.pdf')

## community cluster
distribution_wide <- select(distribution, site, latin, count) %>% 
  group_by(site, latin) %>% summarise(count = sum(count)) %>% ungroup()

distribution_wide <- tidyr::spread(distribution_wide, latin, count)
distribution_wide[is.na(distribution_wide)] <- 0
distribution_wide <- tibble::column_to_rownames(distribution_wide, 'site')
distribution_wide <- distribution_wide[c(1,4,2,3),]

library(vegan)
library(ape)
df_dist <- vegdist(distribution_wide, method = 'bray')
df_hc1 <- hclust(df_dist, method="average") %>% as.dendrogram()
par(mar=c(3,0,3,0))
pdf('Result/cluster.pdf', width = 9.5, height = 4)
plot(df_hc1, type = "rectangle", horiz = T)
dev.off()

# species diversity estimation
library(iNEXT)
library(ggpubr)
out <- iNEXT(t(distribution_wide), q = c(0,1,2), datatype ="abundance",
             knots = 40,se = TRUE,conf = 0.95,nboot = 50)
p1 <- ggiNEXT(out,type = 2, color.var = 'Assemblage')+
  theme_bw() + theme(plot.margin = unit(c(0.5,1,1,1),'cm'))
p2 <- ggiNEXT(x=out, type = 3, facet.var = 'Assemblage')+
  theme_bw() + theme(plot.margin = unit(c(0.2,1,1,1),'cm'),
                     strip.background = element_rect(fill = "white"))
p_inext <- ggarrange(p1, p2, ncol = 1, labels = c("A", "B"))
p_inext
ggsave('Result/inext.pdf', width = 28, height = 20, units = 'cm')

# abund_rank
abund_rank <- group_by(distribution, site, latin) %>% 
  summarise(abund = sum(count)) %>% 
  mutate(abund_percentage = round(abund / sum(abund), 5)) %>% 
  arrange(desc(abund), .by_group = T) %>% 
  slice_head(n = 20) %>% left_join(distinct(select(distribution, X.U.FEFF.species, latin)), 
                                   by = c('latin'))
abund_rank$site <- factor(abund_rank$site,
                            levels = c('Meiji', 'Yongshu', 
                                       'Qilian', 'Yongle'))
latin_order <- unique(abund_rank$latin)[order(unique(abund_rank$latin))]

ggplot(abund_rank, aes(x = abund_percentage, y = forcats::fct_rev(latin))) + 
  geom_bar(aes(fill = site), stat = 'identity') + theme_bw() + 
  ylab('species') + xlab('Relative abundance') + 
  scale_y_discrete(position="right")+
  theme(legend.position = 'none',
        axis.text.y = element_text(face = 'italic'),
        panel.grid.minor = element_blank())+
  facet_grid(~ site)+
  scale_fill_manual(values = c('#f8766d', '#f8766d', '#00bfc4', '#00bfc4'))
ggsave('Result/abund_rank.pdf', width = 25, height = 18.5, units = 'cm')
write.csv(abund_rank, 'Result/abund_rank_20.csv', row.names = F)

## trophic group
# fishbase traits download
# source('fishbase_download.R')
# fb_result <- data.frame()
# for (i in species$Latin) {
#   result_i <- fishbase_downloading(i)
#   fb_result <- rbind(fb_result, result_i)
# }
# write.csv(fb_result, 'fb_result.csv', row.names = F)
species_trophic <- read.csv('trophic_group.csv', header = T,
                            encoding="UTF-8", fileEncoding="UTF-8-BOM")
fb_result <- read.csv('fb_result.csv')

fb_result$TL_group <- cut(fb_result$MaxLengthTL, 
                          breaks = c(-Inf, 10, 20, 30, 40, 50, 100, Inf),
                          labels = c("< 10","10-20","20-30","30-40","40-50",
                                     "50-100", "> 100"),
                          right = FALSE)

total_trophic <- left_join(total, species_trophic, by = c('latin')) %>% 
  left_join(fb_result, by = c('latin' = 'Species')) %>% 
  dplyr::filter(!is.na(Primary_functional_group))%>% 
  group_by(site) %>% mutate(abund = count/sum(count))

trophic_abund_table <- group_by(total_trophic, site, Primary_functional_group,
                          Secondary_functional_group) %>% 
  summarise(abund_sum = sum(abund))
trophic_species_table <- group_by(total_trophic, site, Primary_functional_group,
                          Secondary_functional_group) %>% 
  summarise(species_number = n_distinct(latin)) %>% 
  mutate(sp_prot = case_when(
    site == 'Meiji' ~ species_number/225,
    site == 'Yongshu' ~ species_number/182,
    site == 'Qilian' ~ species_number/114,
    site == 'Yongle' ~ species_number/72,
  ))

trophic_abund_table$site <- factor(trophic_abund_table$site,
                          levels = c('Meiji', 'Yongshu', 
                                     'Qilian', 'Yongle'))

trophic_species_table$site <- factor(trophic_species_table$site,
                             levels = c('Meiji', 'Yongshu', 
                                        'Qilian', 'Yongle'))

group_rank <- c("piscivore",
                "generalist carnivore (fish and invertebrates)",
                "invertivore",
                "scaly fish",
                "ectoparasite feeder",
                "omnivore",
                "planktivore",
                "corallivore",
                "browsers",
                "scraper/excavators",
                "grazer/detritivores"
                )
trophic_abund_table$Secondary_functional_group <- factor(trophic_abund_table$Secondary_functional_group,
                                                   levels = group_rank)
trophic_species_table$Secondary_functional_group <- factor(trophic_species_table$Secondary_functional_group,
                                                         levels = group_rank)

trophic_theme <- theme(legend.text = element_text(size = 6), 
                 legend.title = element_text(size = 5), 
                 legend.key.size = unit(0.5, 'cm'),
                 axis.title = element_text(size = 8),
                 axis.title.x = element_blank())
library(ggalluvial)
library(ggpubr)
p_trophic_abund <- ggplot(trophic_abund_table,
                          aes(x = site, weight = 100*abund_sum, fill = Secondary_functional_group))+
  geom_bar(position = "stack", width = 0.75) + labs(y = 'Relative abundance') + theme_bw() +
  trophic_theme
p_trophic_species <- ggplot(trophic_species_table,
                            aes(x = site, weight = sp_prot, fill = Secondary_functional_group))+
  geom_bar(position = "stack", width = 0.75) + labs(y = 'Relative species number') + theme_bw() +
  trophic_theme
ggarrange(p_trophic_species, p_trophic_abund,
          ncol = 2, labels = c("A", "B"),
          common.legend = TRUE, legend = 'bottom')
ggsave('Result/trophic_group_stack.pdf', width = 18, height = 11, units = 'cm')

# debug
# trophic_species <- dplyr::filter(total_trophic, site=='Yongshu') %>% pull(latin) %>% unique()
# distribution_species <- dplyr::filter(distribution, site=='Yongshu') %>% pull(latin) %>% unique()
# setdiff(distribution_species, trophic_species)

## violin plot
schooling_species <- group_by(abund_rank, site) %>% slice_head(n = 5) %>%
  pull(latin) %>% unique()
violin_total <- filter(total, latin %in% schooling_species)
violin_total$site <- factor(violin_total$site,
                            levels = c('Meiji', 'Yongshu', 
                                       'Qilian', 'Yongle'))
my_comparisons <- list(c("Meiji","Yongshu"),
                       c("Qilian", "Yongle"))
library(ggpubr)
p_sch <- ggplot(violin_total, aes(x = site, y = count)) +
  geom_violin(aes(fill = site))+
  stat_summary(fun = median, geom = "point")+
  stat_compare_means(comparisons = my_comparisons, method = 'wilcox.test')+
  facet_wrap(~ latin, ncol = 3, scales = "free_y")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "italic"),
        legend.position = "none")
p_sch
# ggsave('Result/violin_plot_median.pdf', height = 9, width = 10)

## shoaling size distribution
shoaling_mean <- dplyr::filter(total, !is.na(latin)) %>% group_by(site, latin) %>%
  summarise(shoaling_mean = mean(count)) %>%
  group_by(site) %>% 
  arrange(desc(shoaling_mean), .by_group = TRUE) %>% slice_head(n = 30) %>% 
  tidyr::unite('site_species', site, latin, sep = '_', remove = FALSE) %>%
  mutate(site_species = factor(site_species, levels = site_species)) %>%
  left_join(species_trophic, by = 'latin')

shoaling_mean$site <- factor(shoaling_mean$site,
                                     levels = c('Meiji', 'Yongshu', 
                                                'Qilian', 'Yongle'))
ggplot(shoaling_mean,
       aes(x = site_species, y = shoaling_mean, fill = Primary_functional_group)) +
  geom_bar(stat = 'identity') +  
  facet_wrap(~ site, ncol = 4, scales="free_x")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = 'italic'),
        panel.grid = element_blank(),
        legend.position = "bottom")
ggsave('Result/shoaling_mean_distribution.pdf', height = 6, width = 18)

## anova
library(FSA)
nonpar_result <- data.frame()
for (j in unique(abund_rank$latin)) {
  # j = 'Hemitaurichthys polylepis'
  tmp <- dplyr::filter(total, latin %in% j) %>% mutate(site = as.factor(site))
  if (n_distinct(tmp$site) >= 3) {
    dunn <- dunnTest(count ~ site, tmp, method = "bh")
    print(j)
    tmp_result <- as.data.frame(dunn$res)
    tmp_result$species <- rep(j, nrow(tmp_result))
    nonpar_result <- bind_rows(nonpar_result, tmp_result)
  }else if(n_distinct(tmp$site) == 2){
    kt <- kruskal.test(count ~ site, tmp)$p.value
    kt_result <- data.frame(Comparison=paste(unique(tmp$site), collapse = ' - '),
                            species=j,
                            P.unadj=kt)
    nonpar_result <- bind_rows(nonpar_result, kt_result)
  }
}
write.csv(nonpar_result, 'Result/nonpar_result_schooling.csv', row.names = F)
