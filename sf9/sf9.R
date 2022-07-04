library(data.table)
library(tidyverse)
library(rtracklayer)
library(tximport)
library(eulerr)
library(pals)
library(ggpubr)
library(pals)
library(patchwork)
library(ggplotify)
library(enrichR)
library(gghalves)
library(ggh4x)

setEnrichrSite("Enrichr")
dbs <- listEnrichrDbs()
db.use <- c('WikiPathway_2021_Human', 'BioPlanet_2019')


g <- import.gff3('../data/genes/gencode.v36.annotation.gff3.gz')
e2g <- as_tibble(g[g$type == 'gene']) %>%
  dplyr::select(gene_id, gene_name) %>% deframe()
tx2gene <- g[g$type == 'transcript'] %>% as_tibble() %>%
  dplyr::select(transcript_id, gene_id) %>% distinct()
fs <- list.files('../data/rna/gbm_culture', full.names = T) %>%
  setNames(., sub('.sf', '', basename(.)))
txi <- tximport(fs, 'salmon', tx2gene = tx2gene)



signif.num <- function(x) {
  as.character(symnum(x, corr = FALSE, na = FALSE, legend = FALSE,
                      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                      symbols = c("****", "***", "**", "*", "ns")))
}
tclrs <- c('#e45756', '#4c78a8')

cclrs <- kelly()[c(3,7,6,2)]
pd <- readRDS('../data/pmtcgi/GBM.embed.rds') %>%
  mutate(cl = case_when(clust == '0' ~ '1',
                        clust == '1_1' ~ '4',
                        clust %in% c('1_0_0', '1_0_3_0') ~ '2b',
                        clust %in% c('1_0_1', '1_0_3_1') ~ '2a',
                        clust %in% c('1_0_2_0_0', '1_0_2_1_1') ~ '3b',
                        clust %in% c('1_0_2_0_1', '1_0_2_1_0') ~ '3a',
                        T ~ NA_character_),
         cl2 = c('1' = 'Active', '2a' = 'cPRC1', '2b' = 'cPRC1',
                 '3a' = 'PRC2', '3b' = 'PRC2', '4' = 'Other')[cl] %>%
           factor(c('Active','cPRC1','PRC2','Other'))) %>%
  merge(mutate(fread('../data/pmtcgi/hg38.5kb.bed',
                     col.names = c('chr', 'start', 'end', 'name', 'score', 'strand')),
               idx = 1:n()), by = 'idx') %>%
  na.omit()


d2 <- fread('../data/pmtcgi/GBM.full.csv',
            drop = c(names(pd)[3:13]))

e <- list(BT245 = c( "SD_BT245_p18", 'SD_BT245_p18'),
          DIPGXIII = c("DIPGXIII_C12-P2", "DIPGXIII_C8-P2"),
          HSJ019 = c("HSJ019_par_p10", "HSJ019_par_p11")) %>%
  lapply(function(x) {
    txi$abundance[,x] %>%
      rowMeans(na.rm = T)
  }) %>%
  bind_cols() %>%
  mutate(name = rownames(txi$abundance))

fread('../data/pmtcgi/hg38.5kb.bed') %>%
  mutate(idx = 1:n()) %>%
  select(idx, name = V4) %>%
  merge(e) %>%
  pivot_longer(-c(name, idx), names_to = 'line', values_to = 'Expression') %>%
  mutate(Expression = log2(Expression + 1)) %>%
  merge(pd, by = c('line', 'idx')) %>%
  ggplot(aes(x = cl2, y = Expression)) +
  geom_boxplot(aes(color = cl2, fill = cl2), notch = T) +
  stat_summary(geom = "crossbar", width=0.2, fatten=0, color="white",
               fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
  facet_wrap(~line, nrow = 1) +
  scale_y_log10('Gene expression (TPM)', breaks = c(0.01, .1, 1, 10, 100),
                labels = c('0.01', '0.1', '1', '10', '100')) +
  scale_color_manual(values = cclrs) +
  scale_fill_manual(values = cclrs) +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black'),
        strip.background = element_rect(fill = 'black'),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        panel.grid.major = element_line(color ='grey75', linetype = 'dashed'),
        strip.text = element_text(color = 'white', size = 11),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) -> p

ggsave('sf9_a.pdf', p, height = 5, width = 4)

pd %>%
  split(., .$cl2) %>%
  lapply(function(x) {
    v <- split(x$idx, x$line)
    u <- unlist(v) %>% unique()
    lapply(v, function(y) u %in% y) %>%
      bind_cols() %>%
      euler(shape = 'ellipse') %>%
      plot(quantities = T) %>%
      as.ggplot()
  }) %>%
  wrap_plots(nrow=2) -> p

ggsave('sf9_b.pdf', p, height = 5, width = 5)


readRDS('../data/pmtcgi/GBM.embed.rds') %>%
  mutate(cl = case_when(clust == '0' ~ '1',
                        clust == '1_1' ~ '4',
                        clust %in% c('1_0_0', '1_0_3_0') ~ '2b',
                        clust %in% c('1_0_1', '1_0_3_1') ~ '2a',
                        clust %in% c('1_0_2_0_0', '1_0_2_1_1') ~ '3b',
                        clust %in% c('1_0_2_0_1', '1_0_2_1_0') ~ '3a',
                        T ~ NA_character_)) %>%
  na.omit() %>%
  dplyr::select(line, cl, idx) %>%
  merge(mutate(fread('../data/pmtcgi/hg38.5kb.bed'), idx = 1:n()), by = 'idx') %>%
  filter(grepl('ENSG', V4)) %>%
  dplyr::select(line, cl, name = V4) %>%
  split(., .$line) %>%
  lapply(function(x) {
    tmp <- split(x$name, x$cl)
    tmp$`3` <- c(tmp$`3a`, tmp$`3b`)
    tmp$`2` <- c(tmp$`2a`, tmp$`2b`)
    tmp
  }) -> gl

lapply(gl, `[[`, '2') %>%
  Reduce(intersect, .) %>%
  e2g[.] %>%
  unname() %>%
  unique() %>%
  enrichr(db.use) -> res

res %>%
  rbindlist(idcol = 'db') %>%
  filter(db == 'WikiPathway_2021_Human') %>%
  filter(Adjusted.P.value < .05) %>%
  arrange(P.value) %>%
  mutate(pway = sub(' WP.*', '', Term) %>% fct_inorder(),
         olap = sub('/.*', '', Overlap) %>% as.integer()) %>%
  ggplot(aes(y = pway, x = Odds.Ratio, color = -log10(Adjusted.P.value))) +
  geom_segment(aes(yend = pway, xend = 0)) +
  geom_point(aes(size = olap)) +
  scale_x_continuous(expand = expansion(c(0, .05))) +
  scale_size('Overlap') +
  scale_color_gradientn(expression(-log[10]~'FDR'), colors = viridis(50)) +
  labs(x =  'Odds ratio', y = 'Pathway') +
  facet_wrap(~'Consensus cPRC1 targets: WikiPathways') +
  theme(legend.direction = "vertical",
        legend.box = "horizontal",
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = 11, color = 'black'),
        axis.text = element_text(size = 11, color = 'black'),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(size = 11, color = 'white'),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        plot.background = element_blank(),
        panel.grid.major.x = element_line(color = 'grey80', linetype = 'dashed'),
        panel.grid = element_blank()) -> p
ggsave('sf9_c.pdf', p, height = 2.42, width = 8.3, bg = 'transparent')




readRDS('../data/pmtcgi/GBM.embed.rds') %>%
  mutate(cl = case_when(clust == '0' ~ '1',
                        clust == '1_1' ~ '4',
                        clust %in% c('1_0_0', '1_0_3_0') ~ '2b',
                        clust %in% c('1_0_1', '1_0_3_1') ~ '2a',
                        clust %in% c('1_0_2_0_0', '1_0_2_1_1') ~ '3b',
                        clust %in% c('1_0_2_0_1', '1_0_2_1_0') ~ '3a',
                        T ~ NA_character_),
         cl2 = sub('[ab]$', '', cl),
         cl2 = c('1' = 'Active', '2' = 'cPRC1', '3' = 'PRC2', '4' = 'Other')[cl2] %>%
           factor(c('Active', 'cPRC1', 'PRC2', 'Other'))) %>%
  na.omit() %>% 
  filter(cl %in% c('3a','3b')) %>%
  mutate(cl2 = c('3a' = 'H3K4me3-', '3b' = 'H3K4me3+')[cl]) %>%
  ggplot(aes(x = cl2, y = H3K4me3, color = cl2, fill = cl2)) + 
  geom_hline(yintercept = 0) +
  geom_half_violin(side = 'l', color = NA, alpha = .5, scale = 'width') +
  geom_boxplot(position = position_nudge(x = .2), width = .3) +
  stat_summary(geom = "crossbar", width=0.2, fatten=0, color="white", position = position_nudge(x = .2),
               fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
  ylab('H3K4me3 enrichment') +
  facet_wrap(~'PRC2') +
  coord_cartesian(ylim = c(-5,7.5)) +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        strip.background = element_rect(fill = 'black'),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = 'none',
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_line(color = 'grey70', linetype = 'dashed'),
        legend.text = element_text(size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(color = 'white', size = 11)) -> p

ggsave('sf9_d.pdf', p, height = 2.8, width = 1.5, bg = 'transparent')



pdd <- readRDS('../data/pup/GBM.clust.rds') %>%
  filter(set == 'int') %>%
  filter(reg == 'both') 

ann <- pdd %>%
  group_by(line, cond, cl, scl) %>%
  filter(between(x, 10, 12) & between(y, 10, 12)) %>%
  summarise(z = round(mean(z), 2)) %>%
  ungroup()

ggplot(pdd, aes(x, y, fill = z)) +
  geom_raster() +
  geom_label(aes(x = Inf, y = -Inf, hjust =1, vjust = 1, label = z), 
             data = ann, inherit.aes = F,
             alpha = .8, label.size = NA) +
  facet_nested(line + cond ~ cl + scl,
               strip = strip_nested(
                 text_x = elem_list_text(
                   color = c(rep('white', 4), rep('black', 8))
                 ),
                 background_x = elem_list_rect(
                   fill = c(pals::kelly()[c(3,7,6,2)], rep('white', 8))
                 ),
                 background_y = elem_list_rect(
                   fill = c(rep('black', 3), rep(tclrs, 3))
                 )
               )) +
  scale_y_reverse() +
  coord_cartesian(expand = F) +
  scale_fill_gradientn('O/E', colors = coolwarm(), trans = 'log2',
                       oob = scales::squish, limits = c(1/3.5,3.5)) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'right',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.background = element_blank(),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white', size = 11)) -> p
ggsave('sf9_e.pdf', p, height = 5.8, width = 8.8, bg = 'transparent')


