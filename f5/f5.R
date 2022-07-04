library(data.table)
library(tidyverse)
library(patchwork)
library(tximport)
library(pals)
library(gghalves)
library(scattermore)
library(ggh4x)
library(ggdist)
library(ggtext)
library(readxl)
library(DESeq2)
library(GSVA)
library(rtracklayer)
library(ggpubr)
library(survival)
library(shadowtext)
library(ggnewscale)
source('../data/diverging_map.R')
signif.num <- function(x) {
  as.character(symnum(x, corr = FALSE, na = FALSE, legend = FALSE,
                      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                      symbols = c("****", "***", "**", "*", "ns")))
}
tclrs <- c('#e45756', '#4c78a8')
cmap <- diverging.colormap(seq(0, 1, .01),
                           rgb1 = hex2RGB(tclrs[2]),
                           rgb2 = hex2RGB(tclrs[1]),
                           outColorspace = "sRGB") %>%
  {.[. > 1] <- 1; .} %>%
  {rgb(.[,1], .[,2], .[,3])}

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
p1 <- ggplot(pd, aes(x = V1, y = V2)) +
  geom_scattermore(aes(color = cl2)) +
  labs(x = 'UMAP 1', y = 'UMAP 2') +
  facet_wrap(~'K27M pGBM CGIs & promoters') +
  scale_color_manual('Cluster', values = cclrs,
                     guide = guide_legend(keyheight=0.2,
                                          default.unit="inch")) +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        strip.background = element_rect(fill = 'black'),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.text = element_text(size = 11),
        strip.text = element_text(color = 'white', size = 11)) 

p2 <- pd %>%
  select(V1, V2, H3K27me3, line, idx, CBX2, H3K4me3, H3K27ac) %>%
  pivot_longer(c(H3K27me3, H3K27ac, H3K4me3, CBX2), 
               names_to = 'mark', values_to = 'v') %>%
  ggplot(aes(x = V1, y = V2)) +
  geom_scattermore(aes(color = v)) +
  labs(x = 'UMAP 1', y = 'UMAP 2') +
  facet_wrap(~mark, nrow = 2) +
  scale_color_gradientn('Enrichment', colours = coolwarm(), limits = c(-3, 3),
                        oob = scales::squish, breaks = c(-3, 0, 3),
                        guide = guide_colorbar(barheight = 2, barwidth = .5,
                                               title.position = 'top',
                                               title.vjust = 0,)) +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        strip.background = element_rect(fill = 'black'),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = 'right',
        legend.justification = 'bottom',
        legend.title = element_text(angle = 270),
        axis.title = element_blank(),
        legend.text = element_text(size = 11),
        strip.text = element_text(color = 'white', size = 11)) 


{ wrap_plots(p1,p2,nrow =1) &
    theme(plot.background = element_blank()) } %>%
  ggsave('f5_a.pdf', height = 5.2, width = 5.9, bg = 'transparent')






g <- import.gff3('../data/genes/gencode.v36.annotation.gff3.gz')
e2g <- as_tibble(g[g$type == 'gene']) %>%
  dplyr::select(gene_id, gene_name) %>% deframe()
tx2gene <- g[g$type == 'transcript'] %>% as_tibble() %>%
  dplyr::select(transcript_id, gene_id) %>% distinct()
fs <- list.files('../data/rna/gbm_culture', full.names = T) %>%
  setNames(., sub('.sf', '', basename(.)))
txi <- tximport(fs, 'salmon', tx2gene = tx2gene)


d2 <- fread('../data/pmtcgi/GBM.full.csv',
            drop = c(names(pd)[3:13]))

e <- list(BT245 = c( "SD_BT245_p18", 'SD_BT245_p18' ),
          DIPGXIII = c( "DIPGXIII_C12-P2", "DIPGXIII_C8-P2"),
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
  select(-name) %>%
  merge(d2, all = T) -> d3

m1 <- c('Expression', 'H3K4me3', 'H3K27ac',
        'H33' ,'H3K27M', 'CTCF', 'SMC1',
        'H3K36me3', 'H3K36me2', 'H3K9me3',
        'RING1B', 'H2Aub', 'CBX2', 'SUZ12', 'H3K27me3','H3K27me2', 'H3K27me1')

pdd <- pd %>%
  merge(d3, by = c('line','idx')) %>%
  dplyr::select(all_of(m1), cl2) %>%
  pivot_longer(-cl2, names_to ='mark', values_to = 'v') %>%
  group_by(cl2, mark) %>%
  summarise(v = mean(v, na.rm = T), .groups = 'drop') %>%
  mutate(mark = factor(mark, m1),
         x = as.numeric(cl2),
         y = as.numeric(mark))


pp1 <- ggplot(pdd, aes(x = x, y = y, fill = v)) +
  geom_tile(aes(fill = v), data = ~subset(., mark=='Expression')) +
  scale_fill_gradientn('log2(TPM + 1)', colors = cividis(25), breaks = c(1, 3),
                       limits = range(pdd$v[pdd$mark == 'Expression']),
                       guide = guide_colorbar(barheight = .5, barwidth = 2,
                                              title.position = 'left',
                                              title.vjust = 1, order = 1)) +
  new_scale_fill() +
  geom_tile(aes(fill = v), data = ~subset(., mark !='Expression')) +
  geom_hline(yintercept = 1.5, size = 1, color = 'white') +
  scale_fill_gradientn('Enrichment', colors = coolwarm(25), 
                       limits = c(-3, 3), breaks = c(-3, 0, 3), oob = scales::squish,
                       guide = guide_colorbar(barheight = .5, barwidth = 2,
                                              title.position = 'left',
                                              title.vjust = 1, order = 2)) +
  scale_x_reverse(breaks = unique(pdd$x), labels = unique(pdd$cl2)) +
  scale_y_continuous(breaks = unique(pdd$y), labels = unique(pdd$mark)) +
  #coord_cartesian(expand = F) +
  labs(x = 'Cluster') +
  coord_flip(expand = F) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        strip.background = element_rect(fill = 'black'),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = 'top',
        legend.justification = 'right',
        axis.title = element_blank(),
        legend.text = element_text(size = 11),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        strip.text = element_text(color = 'white', size = 11)) 

ggsave('f5_b.pdf', pp1, height = 3, width = 5)


pp1 <- pd %>% na.omit() %>% 
  ggplot(aes(x = cl2, y = H3K4me3, color = cl2, fill = cl2)) + 
  geom_hline(yintercept = 0) +
  geom_half_violin(side = 'l', color = NA, alpha = .5, scale = 'width') +
  geom_boxplot(position = position_nudge(x = .2), width = .3) +
  stat_summary(geom = "crossbar", width=0.2, fatten=0, color="white", position = position_nudge(x = .2),
               fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
  scale_color_manual(values = cclrs) +
  scale_fill_manual(values = cclrs) +
  ylab('H3K4me3 enrichment') +
  facet_wrap(~'Major clusters') +
  coord_cartesian(ylim = c(-5,7.5)) +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        strip.background = element_rect(fill = 'black'),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = 'none',
        axis.title.x = element_blank(),
        panel.grid.major = element_line(color = 'grey70', linetype = 'dashed'),
        legend.text = element_text(size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(color = 'white', size = 11)) 


pp2 <- pd %>% na.omit() %>% 
  filter(cl %in% c('2a','2b')) %>%
  mutate(cl2 = c('2a' = 'H3K4me3-', '2b' = 'H3K4me3+')[cl]) %>%
  ggplot(aes(x = cl2, y = H3K4me3, color = cl2, fill = cl2)) + 
  geom_hline(yintercept = 0) +
  geom_half_violin(side = 'l', color = NA, alpha = .5, scale = 'width') +
  geom_boxplot(position = position_nudge(x = .2), width = .3) +
  stat_summary(geom = "crossbar", width=0.2, fatten=0, color="white", position = position_nudge(x = .2),
               fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
  #scale_color_manual(values = cclrs) +
  # scale_fill_manual(values = cclrs) +
  ylab('H3K4me3 enrichment') +
  facet_wrap(~'cPRC1') +
  coord_cartesian(ylim = c(-5,7.5)) +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        strip.background = element_rect(fill = cclrs[2]),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = 'none',
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_line(color = 'grey70', linetype = 'dashed'),
        legend.text = element_text(size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(color = 'white', size = 11)) 

{wrap_plots(pp1, pp2, nrow = 1, widths = c(2,1)) &
    theme(plot.background = element_blank())} %>%
  ggsave('f5_c.pdf', ., height = 3.6, width = 5)

pdd <- readRDS('../data/pup/GBM.clust.rds') %>%
  filter(set == 'int') %>%
  filter(reg == 'both') %>%
  filter(line == 'BT245') %>%
  pivot_wider(names_from = 'cond', values_from = 'z') %>%
  mutate(oe = K27M/KO) %>%
  filter(cl == 'cPRC1' | scl == 'All') %>%
  mutate(grp = ifelse(cl == 'cPRC1' & scl != 'All', scl, as.character(cl)) %>%
           factor(c('Active','cPRC1','H3K4me3-','Other','PRC2','H3K4me3+'))) 
  
ann <- pdd %>%
  group_by(grp) %>%
  filter(between(x, 10, 12) & between(y, 10, 12)) %>%
  summarise(z = round(mean(oe), 2)) %>%
  ungroup()

ggplot(pdd, aes(x, y, fill = oe)) +
  geom_raster() +
  geom_label(aes(x = Inf, y = Inf, hjust =1, vjust = 1, label = z),
             data = ann, inherit.aes = F,
             alpha = .8, label.size = NA) +
  facet_wrap2(~grp) +
  coord_cartesian(expand = F) +
  scale_fill_gradientn(expression(frac('K27M O/E','K27M-KO O/E')),
                       colors = cmap,
                       limits =  c(1/1.9,1.9),
                       oob = scales::squish,
                       trans = 'log2',
                       guide = guide_colorbar(label.position = 'bottom',
                                              title.vjust = 1,
                                              barheight = 2),
                       breaks = scales::pretty_breaks(3)) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'bottom',
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.background = element_blank(),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white', size = 11)) -> p 
ggsave('f5_d.pdf', p , height = 3.3, width = 3.6, bg = 'transparent')


signif.num <- function(x) {
  as.character(symnum(x, corr = FALSE, na = FALSE, legend = FALSE,
                      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                      symbols = c("****", "***", "**", "*", "ns")))
}

load('../data/pmtcgi/GBM.5kb.rda')

dif <- cts %>%
  Map(function(x, nm) {
    1e6 * x / lsz[nm]
  }, ., names(.)) %>%
  split(., sub('_[0-9]$','', names(.))) %>%
  lapply(function(x) {
    bind_cols(x) %>%
      rowMeans()
  }) %>% bind_cols() %>%
  mutate(idx = 1:n()) %>%
  pivot_longer(-idx, names_to = 'track', values_to = 'v') %>%
  filter(grepl('EV', track)) %>%
  separate(track, c('line', 'geno', 'cons', 'mark'), '_') %>%
  dplyr::select(-cons) %>%
  pivot_wider(names_from = 'geno', values_from = 'v') %>%
  mutate(LFC = log2(KO / K27M)) %>%
  dplyr::select(-K27M, -KO) %>%
  pivot_wider(names_from = 'mark', values_from = 'LFC')

ma <- cts %>%
  Map(function(x, nm) {
    1e6 * x / lsz[nm]
  }, ., names(.)) %>%
  split(., sub('_[0-9]$','', names(.))) %>%
  lapply(function(x) {
    bind_cols(x) %>%
      rowMeans()
  }) %>% bind_cols() %>%
  mutate(idx = 1:n()) %>%
  pivot_longer(-idx, names_to = 'track', values_to = 'v') %>%
  filter(grepl('EV', track)) %>%
  separate(track, c('line', 'geno', 'cons', 'mark'), '_') %>%
  dplyr::select(-cons) %>%
  pivot_wider(names_from = 'geno', values_from = 'v')


scl <- .4
shf <- .15
tclrs <- c('#e45756', '#4c78a8')


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
  merge(ma, by = c('line', 'idx')) %>%
  filter(is.finite(K27M) & is.finite(KO)) %>% 
  filter(!(mark %in% c('IgG', 'Input'))) %>%
  pivot_longer(-c(line, idx, cl, mark), names_to = 'geno', values_to = 'v') %>%
  mutate(v = log10(v + 1)) %>%
  filter(line == 'BT245') %>%
  filter(mark %in% c('K4me3')) %>%
  arrange(cl) %>%
  mutate(mark = paste0('H3', mark),
         clr = cclrs[as.integer(sub('[ab]$', '', cl))],
         clu = c('1' = 'Active',
                 '2a' = 'cPRC1, H3K4me3-',
                 '2b' = 'cPRC1, H3K4me3+',
                 '3a' = 'PRC2, H3K4me3-',
                 '3b' = 'PRC2, H3K4me3+',
                 '4' = 'Other')[cl],
         clu = sprintf("<span style='color:%s'>%s</span>", clr, clu) %>% fct_inorder()) %>%
  ggplot(aes(x = clu, y = v, fill = geno, color = geno)) +
  geom_half_violin(side = 'l', data = ~subset(., geno == 'K27M'), color = NA,
                   position = position_nudge(x = -shf), alpha = .3,
                   show.legend = F, scale = 'width') +
  stat_pointinterval(data = ~ subset(., geno == 'K27M'), show.legend = F,
                     position = position_nudge(x = -shf)) +
  geom_half_violin(side = 'r', data = ~subset(., geno == 'KO'), color = NA,
                   position = position_nudge(x = shf), alpha = .3, 
                   show.legend = T, scale = 'width') +
  stat_pointinterval(data = ~ subset(., geno == 'KO'), show.legend = T,
                     position = position_nudge(x = shf)) +
  facet_wrap(~'Promoter/CGI H3K4me3', scales = 'free') +
  scale_size(guide = guide_none()) +
  ylab('CPM') +
  facetted_pos_scales(y = list(
    scale_y_continuous(limits = c(0, 2.5), expand = expansion(c(0, .05))),
    scale_y_continuous(limits = c(0, 3), expand = expansion(c(0, .05)))
  )) +
  scale_color_manual(values = tclrs) +
  scale_fill_manual(values = tclrs) +
  theme(plot.background = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white', size = 11),
        axis.text = element_text(color = 'black', size = 11),
        axis.title.x = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.position = c(.65,1),
        legend.justification = c(.6,1),
        axis.text.x = element_markdown(angle = 45, hjust = 1),
        panel.grid.major.y = element_line(color = 'grey80', linetype = 'dashed')) -> p
ggsave('f5_e.pdf', p, height = 3.8, width = 3.2)



cts <- fread('../data/rna/bulk_tumor/Ensembl.ensGene.exon.raw.tsv.gz') %>% 
  column_to_rownames('V1')

md <- read_xlsx('../data/rna/bulk_tumor/info.samples.ext.xlsx') %>%
  filter(Source == 'Nada') %>%
  filter(Nickname %in% colnames(cts)) %>%
  column_to_rownames('Nickname')

cts <- cts[,rownames(md)]  

dds <- DESeqDataSetFromMatrix(cts, md, ~Genotype)
vsd <- vst(dds) %>%
  assay() %>%
  as.data.frame() %>%
  rownames_to_column('id') %>%
  mutate(id = sub('.*:', '', id)) %>%
  group_by(id) %>%
  summarise_all(mean) %>%
  ungroup() %>%
  column_to_rownames('id') %>%
  as.matrix()

g <- import.gff3('../data/genes/gencode.v36.annotation.gff3.gz')
e2g <- as_tibble(g[g$type == 'gene']) %>% dplyr::select(gene_id, gene_name) %>% deframe()

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
  merge(mutate(fread('../data/pmtcgi/hg38.5kb.bed'), idx = 1:n())) %>%
  filter(grepl('^ENSG', V4)) %>%
  split(., .$line) %>%
  lapply(function(x) {
    tmp <- x %>%
      split(., .$cl) %>%
      lapply(function(y) {
        unname(e2g[y$V4])
      })
    tmp$`2` <- unlist(tmp[c('2a', '2b')])
    tmp$`3` <- unlist(tmp[c('3a', '3b')])
    tmp
  }) -> gls

gls <- names(gls[[1]]) %>%
  setNames(., .) %>%
  lapply(function(x) {
    lapply(gls, `[[`, x) %>%
      Reduce(intersect, .)
  })



pd <- gsva(vsd, gls, method = 'ssgsea', 
           min.sz = 15, max.sz = 1e5, verbose = F) %>%
  as.data.frame() %>%
  rownames_to_column('cl') %>%
  pivot_longer(-cl, names_to = 'samp', values_to = 'score') %>%
  merge(rownames_to_column(md,'samp')) %>%
  filter(Genotype %in% c('H3.3K27M', 'H3WT')) %>%
  mutate(grp = sub('.*K27M', 'H3.3 K27M', Genotype) %>%
           sub('H3.3G34RV', 'H3.3 G34R/V', .) %>%
           sub('WT', ' WT', .) %>%
           factor(c('H3.3 K27M', 'H3 WT'))) %>%
  filter(cl %in% c('1', '2a', '2b')) %>%
  mutate(cl1 = c('1' = 'Active', '2a' = 'cPRC1',
                 '2b' = 'cPRC1')[cl],
         cl2 = c('1' = '', '2a' = 'H3K4me3-', '2b' = 'H3K4me3+')[cl])

yr <- diff(range(pd$score))

ann <- pd %>%
  compare_means(score ~ grp, ., group.by = 'cl') %>%
  merge(distinct(pd, cl, cl1, cl2)) %>%
  mutate(symb = as.character(signif.num(p.adj)),
         y = case_when(cl == '1' & group2 == 'H3.3 G34R/V' ~ .8 - yr * .1,
                       cl == '1' & group1 == 'H3.3 G34R/V' ~ .8 - yr * .3,
                       cl == '1'  ~ .8 - yr * .2,
                       cl == '2a' & group2 == 'H3.3 G34R/V' ~ .1 + yr * .1,
                       cl == '2a' & group1 == 'H3.3 G34R/V' ~ .1 + yr * .3,
                       cl == '2a'  ~ .1 + yr * .2,
                       cl == '2b' & group2 == 'H3.3 G34R/V' ~ .3 - yr * .1,
                       cl == '2b' & group1 == 'H3.3 G34R/V' ~ .3 - yr * .3,
                       cl == '2b'  ~ .3 - yr * .2))

p <- ggplot(pd, aes(x = grp, y = score, color = grp, fill = grp)) +
  geom_boxplot(outlier.colour = NA, width = .3,
               position = position_nudge(x = -.25)) +
  stat_summary(geom = "crossbar", width = 0.2, fatten = 0, color = "white", 
               fun.data = function(x) c(y=median(x), 
                                        ymin=median(x),
                                        ymax=median(x)),
               position = position_nudge(x = -.25)) +
  geom_half_point(shape = 16, size = 1) +
  facet_nested(.~cl1 + cl2,
               strip = strip_nested(background_x = elem_list_rect(fill = c(cclrs[1:2], 'black', scales::hue_pal()(2))))) +
  labs(x = 'pHGG subgroup', y = 'ssGSEA score in tumor RNA-seq') +
  scale_color_manual(values = tclrs) +
  scale_fill_manual(values = tclrs) +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = 'grey80', linetype = 'dashed'),
        axis.text = element_text(size =11, color = 'black'),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white', size = 11),
        legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave('f5_f.pdf', p, height = 3.6, width = 3.3)


dat <- read_xlsx('../data/quant/Survival_cumulative_2022.xlsx', sheet = 'BT245')
fit <- survfit(Surv(time, status) ~ genotype, data = dat)

survdiff(Surv(time, status) ~ genotype, data = dat) %>%
  {pchisq(.$chisq, length(.$n)-1, lower.tail = FALSE)} %>%
  sprintf('p = %.1g', .) -> pv

fortify(fit, surv.connect = T) %>%
  ggplot(aes(x = time, y = surv, group = strata, color = strata,
             ymax = upper, ymin = lower, fill = strata)) +
  geom_step(size = 2) +
  ggfortify:::geom_confint(data = ~subset(., strata != 'KO'), alpha = .3, color = NA) +
  scale_fill_manual('Genotype', values = tclrs) +
  scale_color_manual('Genotype', values = tclrs) +
  labs(x = 'Days', y = 'Survival probability') +
  facet_wrap(~'Patient-derived xenograft') +
  annotate('label', label = pv, x = -Inf, y = -Inf, hjust = -.1, vjust = -.3) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(.1)) +
  theme(panel.background = element_rect(fill = NA, color = 'black', size = 1),
        plot.background = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = 'grey80', linetype = 'dashed'),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(size = 11, color = 'white'),
        axis.text = element_text(size = 11, color = 'black'),
        legend.position = c(1, .77),
        legend.justification = c(1, .5),
        axis.line = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        legend.key = element_blank()) -> p1

ggsave('f5_g.pdf', p1,height = 3.1, width = 2.5, bg = 'transparent')


fread('../data/rna/sc_xeno/bt245.meta_data.expr.tsv') %>%
  dplyr::select(cellBC, TOP2A, SOX2, MBP, UMAP_1, UMAP_2, genotype) %>%
  pivot_longer(c(TOP2A, SOX2, MBP), names_to = 'gene', values_to = 'v') %>%
  mutate(gene = factor(gene, c('SOX2', 'TOP2A', 'MBP'))) %>%
  ggplot(aes(x = genotype, y = v, fill = genotype, color = genotype)) +
  geom_violin(alpha = .3, color = NA) +
  stat_summary(geom = 'pointrange', fun.data = mean_cl_boot, size = .6) +
  geom_signif(comparisons = list(c('K27M', 'KO')), map_signif_level  = signif.num,
              color = 'black', tip_length = 0, y = 1.8) +
  facet_wrap(~gene) +
  scale_y_continuous('Expression', expand = expansion(c(0, .15))) +
  xlab('Xenograft sample') +
  scale_color_manual(values = tclrs) +
  scale_fill_manual(values = tclrs) +
  coord_cartesian(ylim = c(0, 2)) +
  theme(plot.background = element_blank(),
        panel.grid.major.y = element_line(color = 'grey80', linetype = 'dashed'),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black'),
        strip.background = element_rect(fill = 'black'),
        legend.position = 'none',
        #axis.title.x = element_blank(),
        axis.text = element_text(size = 11, color = 'black'),
        strip.text = element_text(color = 'white', size = 11)) -> p2

ggsave('f5_h.pdf', p2 + theme(axis.text.x = element_text(angle = 45, hjust = 1)),
       height = 3.3, width = 2.4, bg = 'transparent')

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
  filter(line == 'BT245') %>%
  #mutate(line = sub('DIPG13', 'DIPGXIII', line)) %>%
  group_by(line) %>%
  mutate(`2a` = (`2a` - min(`2a`)) / diff(range(`2a`))) %>%
  ungroup() %>%
  ggplot(aes(x = cor_cc, y = genotype,  color = `2a`)) +
  geom_point(aes(size = prop), shape = 16) +
  geom_shadowtext(aes(label = round(prop, 1)), color = 'white',
                  position = position_nudge(y = .2)) +
  scale_size('% of glial cells in sample\nmatching specific cell type',
             range = c(1, 20),
             breaks = c(10, 20, 30),
             guide = guide_legend(title.position = 'left',
                                  barheight = 6.5, barwidth = .5)) +
  scale_color_gradientn('H3K4me3- cPRC1\ntarget expression',
                        colors = coolwarm(21)[11:21],
                        breaks = c(0, 1),
                        labels = c('Low', 'High'),
                        guide = guide_colorbar(title.position = 'left',
                                               barheight = 6.5, barwidth = .5)) +
  facet_wrap(~'Xenograft composition vs cPRC1 target expression') +
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
ggsave('f5_i.pdf', p, bg = 'transparent', height = 3.8, width = 5.9)


