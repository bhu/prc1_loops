library(data.table)
library(tidyverse)
library(pals)
library(patchwork)
library(gghalves)
library(ggnewscale)
library(ggh4x)
library(survival)
library(readxl)
library(scattermore)
library(shadowtext)

read_csv('../data/hic/metadata.csv') %>%
  filter(biotype == 'primary_tissue') %>%
  mutate(tumorType = factor(tumorType, c('Pediatric high-grade glioma',
                    'Ependymoma',
                    'Medulloblastoma',
                    'Juvenile pilocytic astrocytoma',
                    'Fetal brain',
                    'Adult brain'))) %>%
  arrange(tumorType, subgroup) %>%
  mutate(grp = sprintf('%s (%s)', tumorType, subgroup) %>% fct_inorder(),
         shape = c(15,16,17,18,7,13)[as.numeric(tumorType)]) -> md

clrs <- setNames(kelly(14)[-1],levels(md$grp))
readRDS('../data/pup/GBM.cPRC1.nonBiv.rds') %>%
  filter(between(x, 10, 12) & between(y, 10, 12)) %>%
  group_by(s, r) %>%
  summarise(z = mean(z), .groups = 'drop') %>%
  merge(md, by.x = 's', by.y = 'samp_name') %>%
  filter(r == 'both') %>%
  group_by(grp) %>%
  mutate(med = median(z)) %>%
  ungroup() %>%
  arrange(-med) %>%
  mutate(grp = fct_inorder(grp)) %>%
  ggplot(aes(x = grp, y = z, color = grp, fill = grp)) +
  geom_boxplot(outlier.color = NA, width = .3, position = position_nudge(x = .2)) +
  stat_summary(geom = "crossbar", width=0.15, fatten=0, color="white",
               fun.data = function(x) c(y=median(x), ymin=median(x), ymax=median(x)),
               position = position_nudge(x = .2)) +
  geom_point(position = position_nudge(x = -.2)) +
  scale_y_continuous('Average interaction (O/E) between\ncPRC1 H3K4me3- sites\ndefined in K27M cell lines',
                     breaks = c(2, 4, 6)) +
  facet_wrap(~'Conservation of cPRC1 looping in primary tissues') +
  scale_color_manual(values = clrs) +
  scale_fill_manual(values = clrs) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white', size = 11),
        axis.line.x = element_line(color = 'black'),
        axis.title.x = element_blank(),
        plot.margin = margin(5,5,5,100),
        legend.position = 'none',
        panel.grid.major = element_line(color = 'grey75', linetype = 'dashed'),
        axis.ticks.y = element_blank()) -> p1
p1

readRDS('../data/pup/GBM.cPRC1.nonBiv.rds') %>%
  filter(s %in% c('HSJ050T', 'E2061')) %>%
  filter(between(x, 10, 12) & between(y, 10, 12)) %>%
  filter(r == 'both') %>%
  group_by(s, r) %>%
  summarise(v = round(mean(z), 2), .groups = 'drop') %>%
  merge(readRDS('../data/pup/GBM.cPRC1.nonBiv.rds')) %>%
  mutate(s = ifelse(s == 'HSJ050T', 'Strongest', 'Weakest') %>%
           factor(c('Weakest', 'Strongest'))) %>%
  ggplot(aes(x, y, fill = z)) +
  geom_raster() +
  geom_label(aes(x = Inf, y = Inf, hjust =1, vjust = 1, label = v),
             inherit.aes = F,  alpha = .8, label.size = NA,
             data = ~distinct(., s, v)) +
  facet_wrap(~s) +
  scale_fill_gradientn('O/E', colors = coolwarm(),
                       #breaks = scales::pretty_breaks(3),
                       limits = c(1/7, 7),
                       breaks = c(.15, 1, 7),
                       labels = c('0.15', '1', '7'),
                       oob = scales::squish,
                       trans = 'log2',
                       guide = guide_colorbar(barheight = .5, barwidth = 2.5,
                                              title.position = 'right', title.vjust = 1)) +
  scale_x_continuous(breaks = c(10.5),
                     labels = 'cPRC1 site') +
  scale_y_continuous(breaks = 10.5,
                     labels = 'cPRC1 site') +
  coord_cartesian(expand = F, clip = 'off') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white', size = 11),
        axis.text = element_text(size =11, color = 'black'),
        axis.title = element_blank(),
        legend.position = 'bottom',
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        legend.justification = 'left') -> p2

{wrap_plots(p1,p2) &
    theme(plot.background = element_blank())} %>%
  ggsave('sf10_a.pdf', ., height = 5.42, width = 10.2, bg = 'transparent')



signif.num <- function(x) {
  as.character(symnum(x, corr = FALSE, na = FALSE, legend = FALSE,
                      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                      symbols = c("****", "***", "**", "*", "ns")))
}
tclrs <- c('#e45756', '#4c78a8')





c('H3.3 K27M' = 'h3.3k27m',
  'H3.3 G34R/V' = 'h3.3g34rv',
  #'H3.1/2 K27M' = 'h3.12k27m',
  'H3 WT' = 'h3wt',
  'Fetal hindbrain' = 'human.fetal.hindbrain',
  'Fetal thalamus' = 'human.fetal.thalamus') %>%
  lapply(function(x) {
    c('_umap.tsv', '.single_cells.ssgsea_scores.tsv', '.meta_data.tsv') %>%
      lapply(function(y) {
        sprintf('../data/rna/sc_xeno/%s%s', x, y) %>%
          fread()
      }) %>%
      Reduce(merge, .)
  }) %>%
  bind_rows(.id = 'samp') -> d

clrs <- c('RGC' = '#FFCC00',
          'Glial progenitors' = '#D5D98B',
          'Proliferating OPC' = '#E6F957',
          'OPC' = '#E0DE53',
          'Oligodendrocytes' = '#B4E04E',
          'Astrocytes' = '#00A385',
          'Ependymal' = '#8EE5CF',
          'Neuronal progenitors' = '#FFBDA3',
          'Neurons' = '#135CA0', 
          'Immune' = '#7F7F7F', 
          'Vascular & other' = '#B3B3B3')


cclrs <- kelly()[c(3,7,6,2)]
d %>%
  add_count(samp) %>%
  mutate(samp = factor(samp, c('H3.3 K27M', 'H3 WT',
                               'Fetal hindbrain', 'Fetal thalamus')),
         sz = 10000 / n) %>%
  na.omit() %>%
  dplyr::select(samp, cellID, UMAP_1, UMAP_2, avg_2a, avg_2b) %>%
  pivot_longer(c(avg_2a, avg_2b), names_to = 'clu', values_to = 'v') %>%
  mutate(clu = c('avg_2a' = 'H3K4me3-', 'avg_2b' = 'H3K4me3+')[clu]) %>%
  ggplot(aes(x = UMAP_1, y = UMAP_2, color = v)) +
  geom_scattermore(pointsize = 1.1, data = ~subset(., samp %in% c('H3.3 K27M', 'Fetal thalamus') & clu == 'H3K4me3-')) +
  geom_scattermore(pointsize = 3.3, data = ~subset(., samp %in% c('H3 WT') & clu == 'H3K4me3-')) +
  geom_scattermore(pointsize = 5.5, data = ~subset(., samp %in% c('Fetal hindbrain') & clu == 'H3K4me3-')) +
  facet_nested(samp ~ 'cPRC1 targets' + clu,  scales = 'free',
               strip = strip_nested(background_x = elem_list_rect(fill = c(cclrs[2], scales::hue_pal()(2))))) +
  scale_color_gradientn('Mean expression',
                        colors = coolwarm(21)[11:21],
                        limits = c(0, .06), 
                        oob = scales::squish,
                        breaks = c(0, .06),
                        guide = guide_colorbar(title.position = 'left', order = 1,
                                               title.vjust = 1.2,
                                               barheight = .5, barwidth = 2)) +
  new_scale_color() +
  geom_scattermore(aes(color = v), pointsize = 1.1, 
                   data = ~subset(., samp %in% c('H3.3 K27M', 'Fetal thalamus') & clu == 'H3K4me3+')) +
  geom_scattermore(aes(color = v), pointsize = 3.3,
                   data = ~subset(., samp %in% c('H3 WT') & clu == 'H3K4me3+')) +
  geom_scattermore(aes(color = v), pointsize = 5.5,
                   data = ~subset(., samp %in% c('Fetal hindbrain') & clu == 'H3K4me3+')) +
  scale_color_gradientn('',
                        colors = coolwarm(21)[11:1],
                        limits = c(0, .25), 
                        oob = scales::squish,
                        breaks = c(0, .25),
                        guide = guide_colorbar(title.position = 'left', order = 2,
                                               title.vjust = 1.2,
                                               barheight = .5, barwidth = 2)) +
  labs(x = 'UMAP 1', y = 'UMAP 2') +
  theme(strip.background = element_rect(fill = 'black'),
        strip.text = element_text(size = 11, color = 'white'),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        plot.background = element_blank(),
        axis.text = element_blank(),
        legend.background = element_blank(),
        legend.position = 'bottom',
        legend.text = element_text(size = 11),
        #legend.title = element_text(angle = 270),
        axis.ticks = element_blank()) -> p
ggsave('sf10_b.pdf', p, height = 6.3, width = 3.4, bg = 'transparent')



d %>%
  mutate(samp = factor(samp, c('H3.3 K27M', 'H3 WT',
                               'Fetal hindbrain', 'Fetal thalamus')),
         class = Cell_ontological_class %>% 
           sub('Oligodendrocyte precursors', 'OPC', .) %>%
           factor(names(clrs))) %>%
  na.omit() %>%
  dplyr::select(samp, cellID, class, avg_2a, avg_2b) %>%
  pivot_longer(c(avg_2a, avg_2b), names_to = 'clu', values_to = 'v')  %>%
  mutate(clu = c('avg_2a' = 'H3K4me3-', 'avg_2b' = 'H3K4me3+')[clu]) %>%
  ggplot(aes(x = class, y = v, color = class, fill = class)) +
  geom_boxplot() +
  facet_grid2(clu ~ samp, scales = 'free_y',
              strip = strip_themed(background_y = elem_list_rect(fill = scales::hue_pal()(2)))) +
  scale_fill_manual(values = clrs) +
  scale_color_manual(values = clrs) +
  stat_summary(geom = "crossbar", width=0.4, fatten=0, color="white", 
               fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
  labs(y = 'Geneset-average expression', x = 'Cell type') +
  theme(axis.text = element_text(size = 11, color = 'black'),
        legend.position = 'none',
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(size = 11, color = 'white'),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(color = 'grey80', linetype = 'dashed'),
        plot.background = element_blank()) -> pp
ggsave('sf10_c.pdf', p, height = 6, width = 7, bg = 'transparent')


c('DIPGXIII', 'HSJ019') %>%
  setNames(., .) %>%
  lapply(function(x) {
    dat <- read_xlsx('../data/quant/Survival_cumulative_2022.xlsx', sheet = x)
    fit <- survfit(Surv(time, status) ~ genotype, data = dat)
    
    survdiff(Surv(time, status) ~ genotype, data = dat) %>%
      {pchisq(.$chisq, length(.$n)-1, lower.tail = FALSE)} %>%
      sprintf('p = %.1g', .) -> pv
    
    fortify(fit, surv.connect = T) %>%
      mutate(p = pv)
  }) %>%
  bind_rows(.id = 'line') %>%
  ggplot(aes(x = time, y = surv, group = strata, color = strata,
             ymax = upper, ymin = lower, fill = strata)) +
  geom_step(size = 2) +
  ggfortify:::geom_confint(data = ~subset(., strata != 'KO' | line != 'DIPGXIII'), alpha = .3, color = NA) +
  scale_fill_manual('Genotype', values = tclrs) +
  scale_color_manual('Genotype', values = tclrs) +
  labs(x = 'Days', y = 'Survival probability') +
  facet_grid2(line~'Patient-derived xenograft') +
  geom_label(aes(label = p), x = -Inf, y = -Inf, hjust = -.1, vjust = -.3,
             data = ~distinct(., line, .keep_all = T), color = 'black', fill = 'white') +
  #annotate('label', label = pv, x = -Inf, y = -Inf, hjust = -.1, vjust = -.3) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(.1)) +
  theme(panel.background = element_rect(fill = NA, color = 'black', size = 1),
        plot.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = 'grey80', linetype = 'dashed'),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(size = 11, color = 'white'),
        axis.text = element_text(size = 11, color = 'black'),
        legend.position = c(1,.98),
        legend.justification = c(1,1),
        axis.line = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        legend.key = element_blank()) -> p
ggsave('sf10_d.pdf', p, height = 3.8, width =3.5, bg = 'transparent')





c('bt245','dipg13','hsj019') %>%
  lapply(function(x) {
    fread(sprintf('../data/rna/sc_xeno/%s.meta_data.tsv', x)) %>%
      dplyr::select(1:5) %>%
      `colnames<-`(c('cellID', 'genotype', 'cor_cc', '2a', '2b')) %>%
      mutate(line = x) 
  }) %>%
  rbindlist() %>%
  mutate(cor_cc = factor(cor_cc, c('RGC', 'Glial progenitors', 
                                   'Oligodendrocyte precursors',
                                   'Oligodendrocytes', 'Astrocytes')),
         line = toupper(line)) %>%
  na.omit() %>%
  group_by(cor_cc, genotype, line) %>%
  summarise(num = n(),
            `2a` = mean(`2a`, na.omit = T),
            `2b` = mean(`2b`, na.omit = T)) %>%
  group_by(genotype, line) %>%
  mutate(prop = num / sum(num) * 100) %>%
  ungroup() %>%
  filter(line != 'BT245') %>%
  mutate(line = sub('DIPG13', 'DIPGXIII', line)) %>%
  group_by(line) %>%
  mutate(`2a` = (`2a` - min(`2a`)) / diff(range(`2a`))) %>%
  ungroup() %>%
  ggplot(aes(x = cor_cc, y = genotype, color = `2a`)) +
  geom_point(aes(size = prop), shape = 16) +
  geom_shadowtext(aes(label = round(prop, 1)), color = 'white',
                  position = position_nudge(y = .2)) +
  scale_size('% of glial cells in sample\nmatching specific cell type',
             range = c(1, 15),
             breaks = c(10, 30, 50),
             guide = guide_legend(title.position = 'left',
                                  barheight = 6.5, barwidth = .5)) +
  scale_color_gradientn('H3K4me3- cPRC1\ntarget expression',
                        colors = coolwarm(21)[11:21],
                        breaks = c(0, 1),
                        labels = c('Low', 'High'),
                        guide = guide_colorbar(title.position = 'left',
                                               barheight = 6.5, barwidth = .5)) +
  facet_wrap(~line) +
  labs(y = 'Xenograft genotype', x = 'Cell type') +
  theme(panel.background = element_rect(fill = NA, color = 'black', size = 1),
        plot.background = element_blank(),
        legend.title = element_text(angle = 90),
        legend.position = 'right',
        legend.justification = 'top',
        #panel.grid.major.y = element_line(color = 'grey80', linetype = 'dashed'),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(size = 11, color = 'white'),
        axis.text = element_text(size = 11, color = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        legend.key = element_blank(),
        panel.grid = element_blank()) -> p

ggsave('sf10_e.pdf', p, bg = 'transparent', height = 4, width = 7)
