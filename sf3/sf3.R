library(tidyverse)
library(pals)
library(plotly)
library(cluster)
library(uwot)
library(patchwork)

md <- fread('../data/hic/metadata.csv') %>%
  filter(biotype != 'culture' & !public) %>%
  mutate(tumorType = fct_inorder(tumorType),
         grp = sprintf('%s (%s)', tumorType, subgroup) %>% fct_inorder(),
         shape = c(15,16,17,18,8)[as.numeric(tumorType)])
shps <- deframe(distinct(md[,c('grp', 'shape')]))
clrs <- setNames(kelly(12)[-1],levels(md$grp))


load('../data/hic/eig.rda')
load('../data/hic/tad.rda')

set.seed(42)
d1 <- emat[,md$samp_name] %>%
  na.omit() %>%
  t() %>%
  umap(n_components = 2, metric = 'correlation',) %>%
  as.data.frame() %>%
  cbind(md) %>%
  mutate(kind = 'Compartment score')
d2 <- bias %>% 
  lapply(function(x) is.finite(x$`50000`)) %>%
  bind_cols() %>%
  .[,md$samp_name] %>%
  apply(1, all) %>%
  {bind_cols(lapply(tad, `[[`, 'TADscore'))[.,]} %>%
  .[,md$samp_name] %>%
  na.omit() %>%
  t() %>%
  umap(n_components = 2, metric = 'correlation') %>%
  as.data.frame() %>%
  cbind(md) %>%
  mutate(kind = 'Boundary score')
d3 <- readRDS('../data/hic/hicrep.rds')$`50000`$`5` %>%
  distinct(s1, s2, r) %>%
  {rbind(., .[,c(2,1,3)] %>% `colnames<-`(c('s1','s2','r')),
         tibble(s1 = unique(.$s1), s2 = s1, r = 1))} %>%
  distinct() %>%
  pivot_wider(names_from = 's1', values_from = 'r') %>%
  column_to_rownames('s2') %>%
  .[md$samp_name, md$samp_name] %>%
  {as.dist(1 - .)} %>%
  umap(n_components = 2) %>%
  as.data.frame() %>%
  cbind(md) %>%
  mutate(kind = 'Matrix similarity')

list(d1, d2, d3) %>%
  bind_rows() %>%
  ggplot(aes(x = V1, y = V2, color = grp, shape = grp)) +
  geom_point(size = 5) +
  scale_shape_manual(values = shps) +
  scale_color_manual(values = clrs) +
  labs(x = 'UMAP 1', y = 'UMAP 2') +
  facet_grid('Embedding'~kind) +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white',  size = 11),
        #axis.text = element_text(color = 'black', size = 11),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none',
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = 'grey80', linetype = 'dashed')) -> pplt


list('Compartment score' = emat,
     'Boundary score' = bias %>% 
       lapply(function(x) is.finite(x$`50000`)) %>%
       bind_cols() %>%
       .[,md$samp] %>%
       apply(1, all) %>%
       {bind_cols(lapply(tad, `[[`, 'TADscore'))[.,]}) %>%
  lapply(function(x) {
    x[,md$samp_name] %>%
      na.omit() %>%
      cor() %>%
      {1 - .} %>%
      as.dist() %>%
      silhouette(as.numeric(md$grp), .) %>%
      .[,1:3] %>%
      as.data.frame() %>%
      mutate(across(c(cluster, neighbor), ~factor(levels(md$grp)[.x], levels(md$grp)))) %>%
      cbind(md)
  }) %>% bind_rows(.id = 'kind') -> c1

readRDS('../data/hic/hicrep.rds')$`50000`$`5` %>%
  distinct(s1, s2, r) %>%
  {rbind(., .[,c(2,1,3)] %>% `colnames<-`(c('s1','s2','r')),
         tibble(s1 = unique(.$s1), s2 = s1, r = 1))} %>%
  distinct() %>%
  pivot_wider(names_from = 's1', values_from = 'r') %>%
  column_to_rownames('s2') %>%
  .[md$samp_name, md$samp_name] %>%
  {as.dist(1 - .)} %>%
  silhouette(as.numeric(md$grp), .) %>%
  .[,1:3] %>%
  as.data.frame() %>%
  mutate(across(c(cluster, neighbor), ~factor(levels(md$grp)[.x], levels(md$grp)))) %>%
  cbind(md) %>%
  mutate(kind = 'Matrix similarity', .before = 1) %>%
  rbind(c1) %>%
  group_by(samp_name, kind) %>%
  mutate(w = mean(sil_width),
         i = as.numeric(grp),
         kind = factor(kind,c('Compartment score',
                              'Boundary score',
                              'Matrix similarity'))) %>%
  ungroup() %>%
  arrange(-i, -w) %>% 
  mutate(samp_name = sub('WG_', '', samp_name) %>% fct_inorder()) %>%
  ggplot(aes(y = samp_name, x = sil_width, fill = grp)) +
  geom_vline(xintercept = 0) +
  geom_col() +
  facet_grid('Classification cohesion'~kind, scales = 'free_x') +
  scale_fill_manual(values = clrs) +
  labs(x = 'Silhouette width', y = 'Sample') +
  theme(plot.background = element_blank(),
        panel.background = element_rect(color = 'black', size = 1, fill = NA),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white',size = 11),
        axis.text = element_text(color = 'black', size = 11),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = 'none',
        axis.line = element_blank(),,
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(color = 'grey80', linetype = 'dashed')) -> pplt2

ggsave('sf3.pdf', wrap_plots(pplt, pplt2, nrow = 2) & theme(plot.background = element_blank()),
       height = 8, width = 10, bg = 'transparent')


pplt + theme(legend.position = 'bottom', legend.text = element_text(size = 11)) + 
  guides(color = guide_legend(nrow = 3), shape = guide_legend(nrow = 3)) -> l1
ggsave('sf3_leg.pdf', l1)
