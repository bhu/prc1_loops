library(data.table)
library(tidyverse)
library(rtracklayer)
library(patchwork)
library(ggh4x)
library(ggnewscale)
library(ggforce)
library(matrixStats)

reg <- GRanges('chr22', IRanges(26730000, 27930000))

bws <- list.files('../data/sim', pattern = '.bw$', full.names = T) %>%
  setNames(., sub('.bw', '', basename(.))) %>%
  lapply( function(f) {
    import.bw(f, selection = BigWigSelection(reg))
  })
bins <- tile(reg, 200)[[1]]
seqlevels(bins) <- lapply(bws, function(x) seqlevels(x)) %>%
  Reduce(intersect, .) 
xs <- mid(bins)
scs <- lapply(bws, function(x) {
  seqlevels(x) <- seqlevels(bins)
  mcolAsRleList(x, 'score') %>%
    binnedAverage(bins, ., 'score') %>%
    score() %>%
    tibble(score = .) %>%
    mutate(idx = xs)
}) %>% bind_rows(.id = 'samp')

clrs <- c('K27M' = '#e45756', 'K27M-KO' = '#4c78a8',
          setNames(colorRampPalette(c('#e45756', '#4c78a8'))(7)[2:6],
                   c('+0kb', '+5kb', '+10kb', '+50kb', '+100kb')),
          'Simulated' = '#f28e2b', 'ChIP-seq' = '#59a14f')
pp1 <- scs %>%
  mutate(mark = case_when(grepl('CTCF', samp) ~ 'CTCF', T ~ 'H3K27me3') %>%
           factor(c('H3K27me3', 'CTCF')),
         num = case_when(grepl('K27M_H3K27me3', samp) ~ -1,
                         grepl('kb', samp) ~ as.numeric(sub('kb', '', samp)),
                         grepl('KO_H3K27me3', samp) ~ 101,
                         grepl('K27M_CTCF', samp) ~ 102,
                         T ~ 103)) %>%
  arrange(num) %>%
  filter(!(samp %in% paste0(c(2,3,4,6,7,8,9),'0kb'))) %>% 
  mutate(samp = case_when(grepl('kb', samp) ~ paste0('+', samp),
                          T ~ samp),
         cond = case_when(grepl('K27M', samp) ~ 'K27M',
                          grepl('kb', samp) ~ samp,
                          T ~ 'K27M-KO'),
         type = case_when(grepl('kb', samp) ~ 'Simulated',
                          T ~ 'ChIP-seq'),
         samp = fct_inorder(samp)) %>%
  ggplot(aes(x = idx, y = score, color = cond, fill = cond)) +
  geom_line() +
  geom_area(alpha = .5) +
  geom_label(aes(label = cond), x = Inf, y = Inf, hjust = 1, vjust = 1,
             label.size = NA, fill = '#ffffff99', 
             data = ~distinct(., cond, mark, samp)) +
  geom_label(aes(label = type, color = type), x = -Inf, y = Inf, hjust = 0, vjust = 1,
             label.size = NA, fill = '#ffffff99', 
             data = ~distinct(., type, mark, samp), inherit.aes = F) +
  facet_nested(mark + samp ~ ., scales = 'free') +
  scale_y_continuous('CPM', breaks = scales::pretty_breaks(2)) +
  scale_color_manual(values = clrs) +
  scale_fill_manual(values = clrs) +
  scale_x_continuous(sprintf('%s:%s-%skb', seqnames(reg),
                             formatC(as.integer(start(reg)/1e3), big.mark = ','), 
                             formatC(as.integer(end(reg)/1e3), big.mark = ',')),
                     expand = expansion(0)) +
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        axis.line.y = element_line(color = 'black'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white', size = 11))

ggsave('sf1_a.pdf', pp1, height = 7.5, width = 3.5, bg = 'transparent')

p1 <- readRDS('../data/sim/fcs.exp.rds') %>%
  ggplot(aes(x = log10(shf), y=  fcs, color = cond, group = interaction(cond, mark))) +
  #geom_line()  +
  geom_smooth(aes(linetype = mark), fill = NA, span = .3)+
  scale_x_continuous('Shift',
                     expand = expansion(0),
                     breaks = c(2, 3, 4, 5, 6),
                     labels = c('100b','1kb', '10kb', '100kb', '1mb')) +
  coord_cartesian(xlim = c(1.6, 6)) +
  scale_color_manual('Condition', values = clrs[1:2]) +
  scale_linetype('ChIP target') +
  ylab('Fragment cluster score') +
  facet_wrap(~'Real ChIP-seq') +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black'),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white', size = 11),
        axis.text = element_text(color = 'black', size = 11),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.background = element_blank(),
        legend.position = c(1,1),
        legend.key = element_blank(),
        legend.justification = c(1,1),
        panel.grid.major = element_line(color = 'grey80', linetype = 'dashed'))

p2 <- readRDS('../data/sim/fcs.rds') %>%
  mutate(i = sub('kb', '', ext) %>%
           fct_inseq(),
         f = formatC(shf, format = 'd') %>%
           substr(1,1) %>% 
           as.numeric()) %>%
  arrange(i) %>%
  filter(ext %in% c('0kb', '5kb', '10kb' ,'50kb' ,'100kb')) %>%
  mutate(ext = fct_inorder(paste0('+', ext))) %>%
  filter(between(shf, 200, 1e4) | (shf > 1e4 & f < 5)) %>%
  ggplot(aes(x = log10(shf), y=  fcs, color = ext)) +
  #geom_line()  +
  geom_smooth(aes(linetype = ext), fill = NA, span = .3)+
  scale_x_continuous('Shift',
                     expand = expansion(0),
                     breaks = c(3, 4, 5, 6),
                     labels = c('1kb', '10kb', '100kb', '1mb')) +
  coord_cartesian(xlim = c(2.3, 6)) +
  scale_color_manual('Peak extension', values = clrs[3:7]) +
  scale_linetype('Peak extension') +
  ylab('Fragment cluster score') +
  facet_wrap(~'Simulated H3K27me3') +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black'),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white', size = 11),
        axis.text = element_text(color = 'black', size = 11),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.background = element_blank(),
        legend.position = c(1,1),
        legend.key = element_blank(),
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        legend.justification = c(1,1),
        panel.grid.major = element_line(color = 'grey80', linetype = 'dashed'))

pd1 <- readRDS('../data/sim/fcs.exp.rds') %>%
  filter(shf == 1000)

cns1 <- pd1[,c('fcs','mark','cond')] %>% pivot_wider(names_from = 'cond', values_from = 'fcs')

p3 <- ggplot(pd1, aes(x = cond, y = fcs, color = cond)) +
  geom_point(size = 3) +
  scale_color_manual(values = clrs[1:2]) +
  new_scale_color() +
  geom_link(aes(x = 'K27M', xend = 'K27M-KO', y = K27M, yend = `K27M-KO`, color = stat(index)),
            data = cns1, inherit.aes = F, show.legend = F) +
  scale_color_gradient(low = clrs[1], high = clrs[2]) +
  facet_wrap(~mark) +
  labs(x = 'Condition', y = 'FCS @ 1kb\n("Confinement")') +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black'),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white', size = 11),
        axis.text = element_text(color = 'black', size = 11),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.background = element_blank(),
        legend.position = 'none',
        legend.key = element_blank(),
        legend.justification = c(1,1),
        panel.grid.major = element_line(color = 'grey80', linetype = 'dashed'))


pd2 <- readRDS('../data/sim/fcs.rds') %>%
  mutate(i = sub('kb', '', ext) %>%
           fct_inseq(),
         f = formatC(shf, format = 'd') %>%
           substr(1,1) %>% 
           as.numeric()) %>%
  arrange(i) %>%
  filter(ext %in% c('0kb', '5kb', '10kb' ,'50kb' ,'100kb')) %>%
  mutate(ext = fct_inorder(paste0('+', ext))) %>%
  filter(shf == 1000) 

p4 <- ggplot(pd2, aes(x = ext, y = fcs, color = ext)) +
  geom_point(size = 3) +
  scale_color_manual(values = clrs[3:7]) +
  geom_line(group = 1) +
  facet_wrap(~'Simulated H3K27me3') +
  labs(x = 'Peak extension', y = 'FCS @ 1kb\n("Confinement")') +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black'),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white', size = 11),
        axis.text = element_text(color = 'black', size = 11),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.background = element_blank(),
        legend.position = 'none',
        legend.key = element_blank(),
        legend.justification = c(1,1),
        axis.title.y = element_blank(),
        panel.grid.major = element_line(color = 'grey80', linetype = 'dashed'))


{wrap_plots(p1,p2,p3,p4,nrow=2,
            heights = c(2,1)) &
    theme(plot.background = element_blank())} %>%
  ggsave('sf1_bc.pdf', ., height = 7.5, width = 6.5, bg = 'transparent')


list.files('../data/sim', full.names = T, pattern = 'mat.gz') %>%
  setNames(., sub('.mat.gz', '', basename(.))) %>%
  lapply(function(x) {
    fread(x, skip = 1, drop = 1:6) %>%
      as.matrix() %>%
      {tibble(mu = colMeans(., na.rm = T),
              std = colSds(., na.rm = T),
              num = colSums(!is.na(.)))} %>%
      mutate(x = 1:n())
  }) %>%
  bind_rows(.id = 'ext') -> a

fread('../data/sim/chr22.bed') %>%
  mutate(w = V3 - V2) %>%
  pull(w) -> ws

scl <- seq(10,100,10) %>%
  lapply(function(x) {
    tibble(fac = sum(ws) + length(ws) * x * 2 * 1e3,
           ext = sprintf('%skb', x))
  }) %>%
  bind_rows()

a %>%
  mutate(i = sub('kb', '', ext) %>% as.numeric()) %>%
  merge(scl) %>%
  filter(!(ext %in% c('0kb', '5kb'))) %>%
  arrange(i) %>%
  mutate(ext = fct_inorder(ext)) %>%
  filter(between(x, 501, 1500)) %>%
  ggplot(aes(x = x, y = mu)) +
  geom_line(aes(color = ext)) +
  scale_color_manual('Peak extension', values = pals::turbo(10)) +
  facet_wrap(~'Aggregation of simulated profiles') +
  labs(y = 'CPM') +
  scale_x_continuous(breaks = c(501,1000,1500),
                     labels = c('-500kb', 'H3K27me3 peak', '+500kb')) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white', size = 11),
        axis.text = element_text(color = 'black', size = 11),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        legend.background = element_blank(),
        legend.position = 'right',
        legend.key = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_line(color = 'grey80', linetype = 'dashed')) -> p

ggsave('sf1_d.pdf', p, height = 4, width = 9)


readRDS('../data/sim/fcs.dev.rds') %>% 
  separate(Sample, c('reg','day','mark'), '_') %>% 
  separate(mark, c('mark','rep'),'\\.') %>%
  filter(mark == 'H3K27me3') %>% 
  mutate(reg = sub('facial', 'craniofacial ', reg) %>%
           sub('neural', 'neural ', .)) %>%
  filter(grepl('brain', reg)) %>%
  ggplot(aes(x = day, y = `FCS(1k)`)) +
  geom_smooth(aes(group = 1), color = 'firebrick') +
  geom_point() +
  facet_wrap(~reg, nrow = 1) +
  labs(x = 'Time point', y = 'H3K27me3 confinement') +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white', size = 11),
        axis.text = element_text(color = 'black', size = 11),
        plot.background = element_blank(),
        panel.grid.major = element_line(color = 'grey80', linetype = 'dashed'),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) -> p
ggsave('sf1_e.pdf', p, height = 2.8, width = 5, bg = 'transparent')

