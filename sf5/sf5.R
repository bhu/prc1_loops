library(data.table)
library(tidyverse)
library(ggh4x)
library(ggdist)
library(tximport)
library(rtracklayer)
library(gghalves)
library(readxl)
library(matrixStats)
library(pals)
library(patchwork)

g <- import.gff3('../data/genes/gencode.v36.annotation.gff3.gz')
e2g <- as_tibble(g[g$type == 'gene']) %>%
  dplyr::select(gene_id, gene_name) %>% deframe()
tx2gene <- g[g$type == 'transcript'] %>% as_tibble() %>%
  dplyr::select(transcript_id, gene_id) %>% distinct()
fs <- list.files('../data/rna/gbm_culture', full.names = T) %>%
  setNames(., sub('.sf', '', basename(.)))
txi <- tximport(fs, 'salmon', tx2gene = tx2gene)
pd <- txi$abundance[names(e2g[grep('CBX[24678]', e2g)]),] %>%
  as.data.frame() %>%
  rownames_to_column('id') %>%
  pivot_longer(-id, names_to = 'samp', values_to = 'tpm') %>%
  mutate(name = e2g[id]) %>%
  filter(!grepl('G477', samp)) %>%
  mutate(line = case_when(grepl('HSJ', samp) ~ 'HSJ019',
                          grepl('BT|245', samp) ~ 'BT245',
                          grepl('DIPG|D13', samp) ~ 'DIPGXIII',
                          T ~ NA_character_),
         media = case_when(grepl('SerDif', samp) ~ 'DM',
                           T ~ 'SCM') %>%
           factor(c('SCM','DM')),
         geno = case_when(line == 'HSJ019' & grepl('[cC](8|10)', samp) ~ 'KO',
                          line == 'HSJ019' ~ 'K27M',
                          line == 'DIPGXIII' & grepl('[cC](5|10)', samp) ~ 'KO',
                          line == 'DIPGXIII' ~ 'K27M',
                          line == 'BT245' & grepl('[cC](2|4|5|7)|D2', samp) ~ 'KO',
                          line == 'BT245' ~ 'K27M')) %>%
  filter(media != 'DM')
tclrs <- c('#e45756', '#4c78a8')
ggplot(pd, aes(x = geno, y = tpm)) +
  stat_eye() +
  geom_point(aes(color = geno), alpha = .7, shape = 16) +
  scale_color_manual(values = tclrs) +
  facet_nested(~name + line ) +
  ylab('TPM') +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = 'grey80', linetype = 'dashed'),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(size = 11, color = 'white'),
        axis.text = element_text(color = 'black', size = 11),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        legend.position = 'none') -> p
ggsave('sf5_a.pdf', p, height = 2.8, width = 12, bg = 'transparent')

cgi <- fread('../data/aggr/BT245_KO_H3K27me3.mat.gz', skip = 1, select = 1:6) %>%
  {!grepl('^ENS', .$V4)}

t1k <- list(BT245 = c(KO = "BT245_KO_H3K27me3.mat.gz",
                      K27M = "BT245_K27M_H3K27me3.mat.gz"),
            DIPGXIII = c(KO = "DIPGXIII_KO_H3K27me3.mat.gz",
                         K27M = "DIPGXIII_K27M_H3K27me3.mat.gz"),
            HSJ019 = c(KO = "HSJ019_KO_H3K27me3.mat.gz",
                       K27M = "HSJ019_K27M_H3K27me3.mat.gz")) %>%
  lapply(function(x) {
    rs <- lapply(x, function(y) {
      fread(file.path('../data/aggr', y), skip = 1, drop = 1:6)[,990:1010] %>%
        rowMeans(na.rm = T) %>%
        order(decreasing = T) %>%
        intersect(., which(cgi)) %>%
        head(1000)
    })
    rs$union <- unique(unlist(rs))
    rs$intersect <- intersect(rs$K27M, rs$KO)
    rs
  })

agg <- list.files('../data/aggr', full.names = T) %>%
  tibble(f = .) %>%
  filter(grepl('BT|HSJ|DIPG', f)) %>%
  mutate(s = sub('.mat.gz', '', basename(f))) %>%
  separate(s, c('line','cond','mark')) %>%
  filter(mark %in% c('RING1B','CBX2')) %>%
  split(., .$line) %>%
  lapply(function(x) {
    split(x, x$cond) %>%
      lapply(function(y) {
        split(y$f, y$mark) %>%
          lapply(function(z) {
            fread(z, skip = 1, drop = 1:6) %>%
              {.[t1k[[x$line[1]]]$union,]} %>%
              as.matrix() %>%
              {tibble(mu = colMeans(., na.rm = T),
                      std = colSds(., na.rm = T),
                      num = colSums(is.finite(.)))} %>%
              mutate(idx = 1:n())
          }) %>%
          bind_rows(.id = 'mark')
      }) %>%
      bind_rows(.id = 'cond')
  }) %>%
  bind_rows(.id = 'line')


agg %>%
  mutate(se = std / sqrt(num)) %>%
  ggplot(aes(x = idx, y = mu, color = cond, 
             ymin = mu - se, ymax = mu + se, fill = cond)) +
  geom_ribbon(color = NA, alpha = .25) +
  geom_line() +
  facet_grid2(mark ~ line, 
              #independent = 'y',
              scales = 'free_y') +
  scale_x_continuous(breaks = c(850, 1000, 1150),
                     labels = c('-150kb', 'H3K27me3-\nenriched CGI', '+150kb')) +
  ylab('CPM') +
  coord_cartesian(xlim = c(820, 1180)) +
  scale_color_manual(values = tclrs) +
  scale_fill_manual(values = tclrs) +
  theme(plot.background = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        panel.grid.major = element_line(color = 'grey80', linetype = 'dashed'),
        axis.text = element_text(size = 11, color = 'black'),
        strip.background = element_rect(fill = 'black'),
        legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.title = element_blank(),
        legend.text = element_text(size = 11, color ='black'),
        legend.background = element_blank(),
        legend.key = element_blank(),
        strip.text = element_text(color = 'white', size = 11)) -> p

ggsave('sf5_b.pdf', p,  height = 4, width = 6.92, bg = 'transparent')

ssp <- readRDS('../../sim/gbm.fcs.rds') %>%
  filter(grepl('ring|cbx', smp, ignore.case = T)) %>%
  filter(shf == 10000) %>%
  dplyr::select(smp, fcs)
ssp <- readRDS('../../pp_rna/jab11405.ssp.rds') %>%
  filter(grepl('-RING|-CBX', Sample) & grepl('HSJ', Sample)) %>%
  dplyr::select(smp = Sample, fcs = `FCS(10k)`) %>%
  rbind(ssp)


readRDS('../data/chip/PRC1.ssp.rds') %>%
  separate(samp, c('line', 'geno', 'mark'), '_') %>%
  ggplot(aes(x = geno, y = fcs, fill = geno)) +
  geom_col() +
  facet_grid(mark ~ line, scales = 'free_y') +
  ylab('Confinement') +
  scale_fill_manual(values = tclrs) +
  scale_y_continuous(expand = expansion(c(0, .1))) +
  theme(plot.background = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        #axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        panel.grid.major = element_line(color = 'grey80', linetype = 'dashed'),
        axis.text = element_text(size = 11, color = 'black'),
        strip.background = element_rect(fill = 'black'),
        legend.position = 'none',
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        strip.text = element_text(color = 'white', size = 11)) -> p
ggsave('sf5_c.pdf', p, height = 3.8, width = 3.3)

load('../data/pup/BT245.qtr.rda')

ann <- pups %>%
  filter(between(x, 11, 11) & between(y, 11, 11)) %>%
  group_by(q1, q2, cond) %>%
  summarise(z = mean(z), .groups = 'drop') %>%
  mutate(z = round(z, 2))

tclrs <- c('#e45756', '#4c78a8')
p <- ggplot(pups, aes(x = x, y = y, fill = z)) +
  geom_raster() +
  geom_label(aes(x = Inf, y = Inf, hjust =1, vjust = 1, label = z), data = ann, inherit.aes = F,
             alpha = .5, label.size = NA) +
  geom_label(aes(x = -Inf, y = -Inf, hjust =0, vjust = 0, label = num), 
             alpha = .8, label.size = NA, fill = NA,
             data = nums, inherit.aes = F)+
  scale_fill_gradientn('O/E', colors = coolwarm(), 
                       limits = c(1/1.43,1.43),
                       oob = scales::squish,
                       trans = 'log2', breaks = c(0.7,1,1.4)) +
  facet_nested(q1~cond +q2,strip = strip_nested(
    background_x = elem_list_rect(fill = c(tclrs,
                                           rep('black', 8)))
  )) +
  labs(x = 'H3K27me3-based CGI quartile', y = 'Ring1B-based CGI quartile') +
  coord_cartesian(expand = F) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        legend.background = element_blank(),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(size = 11, color = 'white'),
        axis.text = element_blank(),
        axis.ticks = element_blank())

pd <- readRDS('../data/aggr/BT245.qrt.H3K27ac.rds') %>%
  lapply(function(x) {
    bind_rows(x, .id = 'qrt') %>%
      mutate(mark = 'H3K27ac')
  }) %>%
  bind_rows(.id = 'cond') %>%
  arrange(cond, qrt) %>%
  mutate(grp = paste(cond, qrt)  %>% fct_inorder())

p2 <- pd %>%
  ggplot(aes(x = idx, y = mu, color = cond)) +
  geom_line() +
  facet_grid(mark ~ grp, scales = 'free') +
  labs(y = 'CPM') +
  scale_color_manual(values = tclrs) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +
  facetted_pos_scales(x = list(
    grp == 'K27M 1' ~ scale_x_continuous(breaks = c(1, 100, 200),
                                         labels = c('-20kb', 'CGI', '+20kb'),
                                         expand = expansion(0)),
    T ~ scale_x_continuous(breaks = c(1, 100, 200),
                           labels = c('', '', ''),
                           expand = expansion(0))
  )) +
  theme(legend.position = 'none',
        plot.background = element_blank(),
        panel.background = element_rect(fill = NA,color = 1, size = 1),
        panel.grid = element_blank(),
        strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        axis.text = element_text(color = 'black', size = 11))

{wrap_plots(p, p2, nrow = 2, heights = c(4, 1)) &
    theme(plot.background = element_blank())} %>%
  ggsave('sf5_d.pdf', ., height = 7.5, width = 11, bg = 'transparent')


