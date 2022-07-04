library(data.table)
library(tidyverse)
library(rtracklayer)
library(pbapply)
library(pals)
library(eulerr)
library(ggh4x)
library(RobustRankAggreg)

pks <- readRDS('../data/pks/ctcf.rds')

u <- GRangesList(ctcf) %>%
  unlist() %>%
  reduce()


lapply(ctcf, function(x) overlapsAny(u, x)) %>%
  bind_cols() %>%
  euler(shape = 'ellipse') %>%
  plot(quantities = T)

pp <- readRDS('../data/pup/int_ctcf.conv.rds')

ann <- pp %>%
  filter(between(x, 10, 12) & between(y, 10, 12)) %>%
  group_by(samp, cond, line) %>%
  summarise(z = round(mean(z), 2)) %>% ungroup()

tclrs <- c('#e45756', '#4c78a8')
pp %>%
  ggplot(aes(x,y, fill = z)) +
  geom_raster() +
  geom_label(aes(x = Inf, y = Inf, hjust =1, vjust = 1, label = z), data = ann, inherit.aes = F,
             alpha = .8, label.size = NA) +
  facet_wrap(~samp) +
  scale_fill_gradientn('O/E',
                       colors = coolwarm(),
                       trans = 'log',
                       limits = c(.8, 1.25), oob= scales::squish) +
  facet_grid2(cond ~ line, strip = strip_themed(background_y = elem_list_rect(fill = tclrs))) +
  coord_cartesian(expand = F) +
  theme(panel.background = element_blank(),
        plot.background = element_blank(),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(size = 11, color = 'white'),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank()) -> p 

readRDS('../data/hic/bart3d.rds') %>%
  split(., .$dir) %>%
  lapply(function(x) {
    split(x, x$line) %>%
      lapply(function(y) {
        setNames(y$irwin_hall_pvalue, y$TR) %>%
          sort() %>%
          names()
      }) %>%
      aggregateRanks() %>%
      as_tibble() %>%
      mutate(idx = 1:n())
  }) %>%
  bind_rows(.id = 'dir') %>%
  mutate(dir = c('Decreased' = 'K27M < KO',
                 'Increased' = 'K27M > KO')[dir],
         prc2 = Name %in% c('SUZ12', 'EZH2')) %>%
  ggplot(aes(x = idx, y = Score, color = prc2, label = Name)) +
  geom_point() +
  geom_label_repel(data = ~subset(., idx <= 5),
                   point.padding = 0,
                   min.segment.length = 0,
                   nudge_x = 500,
                   box.padding = 0.3) +
  scale_color_manual(values = c('grey50', 'red')) +
  scale_y_continuous('P-value',
                     trans = trans_new('nlog10', \(x) -log10(x), \(x) 10^(-x)),
                     breaks = c(1,  .1, .01, .001, .0001),
                     labels = trans_format("log10", math_format(10^.x))) +
  xlab('Transcriptional regulator ranked by association with differential interaction') +
  facet_wrap(~dir, scales = 'free_y') +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        panel.grid = element_blank(),
        legend.position = 'none',
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white', size = 11),
        axis.text = element_text(color = 'black', size = 11),
        panel.grid.major = element_line(color = 'grey80', linetype = 'dashed')) -> p
ggsave('sf2_c.pdf', p, height = 3, width = 9)

md <- fread('../data/hic/metadata.csv') %>%
  filter(biotype != 'culture' & !public) %>%
  mutate(tumorType = fct_inorder(tumorType),
         grp = sprintf('%s (%s)', tumorType, subgroup) %>% fct_inorder(),
         shape = c(15,16,17,18,8)[as.numeric(tumorType)])
shps <- deframe(distinct(md[,c('grp', 'shape')]))
clrs <- setNames(kelly(12)[-1],levels(md$grp))


load('../data/hic/eig.rda')
load('../data/hic/tad.rda')

list('Insulation score' =  bias %>% 
       lapply(function(x) is.finite(x$`50000`)) %>%
       bind_cols() %>%
       .[,md$samp] %>%
       apply(1, all) %>%
       {bind_cols(lapply(tad, `[[`, 'TADscore'))[.,]},
     'Compartment score' = emat) %>%
  lapply(function(m) {
    list(BT245 = c('BT245_K27M','BT245_KO'),
         DIPGXIII = c('DIPGXIII_K27M','DIPGXIII_KO'),
         HSJ019 = c('HSJ019_K27M','HSJ019_KO'),
         G477 = c("G477_K27M", "G477_WT")) %>%
      lapply(function(x) {
        m[[x[1]]] - m[[x[2]]]
      }) %>%
      bind_cols() %>%
      na.omit() %>%
      cor() %>%
      as.data.frame() %>%
      rownames_to_column('s1') %>%
      pivot_longer(-s1, names_to = 's2', values_to = 'r')
  }) %>%
  bind_rows(.id = 'kind') -> pd

mid_rescaler <- function(mid = 0) {
  function(x, to = c(0, 1), from = range(x, na.rm = TRUE)) {
    scales::rescale_mid(x, to, from, mid)
  }
}
pd %>%
  mutate(rr = round(r, 2)) %>%
  mutate(across(c(s1,s2), ~sub('\\..*', '', .x))) %>%
  ggplot(aes(x = s1, y = s2, fill = r)) +
  geom_tile() +
  geom_text(aes(label = rr)) +
  facet_wrap(~kind) +
  scale_fill_distiller(expression("Pearson's r for"~Delta*'score'),
                       palette = 'BrBG', rescaler = mid_rescaler(),
                       breaks = c(0, .5, 1),
                       guide = guide_colorbar(title.position = 'left',
                                              barwidth = .5, barheight = 8)) +
  coord_cartesian(expand = F) +
  labs(x = 'Isogenic contrast (K27M vs K27M-KO or WT vs K27M-OE)') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white',size = 11),
        axis.title.y = element_blank(),
        legend.background = element_blank(),
        legend.title = element_text(angle = 90),
        panel.spacing.x = unit(2, 'line'),
        axis.text = element_text(color = 'black', size = 11)) -> p

