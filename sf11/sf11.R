library(data.table)
library(tidyverse)
library(rtracklayer)
library(tximport)
library(DESeq2)
library(pals)
library(patchwork)
library(ggsignif)
library(ggtext)
library(ggh4x)
library(RRHO2)

g <- import.gff3('../data/genes/gencode.v36.annotation.gff3.gz')
e2g <- as_tibble(g[g$type == 'gene']) %>%
  dplyr::select(gene_id, gene_name) %>% deframe()
tx2gene <- g[g$type == 'transcript'] %>% as_tibble() %>%
  dplyr::select(transcript_id, gene_id) %>% distinct()


fs <- list.files('../data/rna/gbm_culture', full.names = T) %>%
  setNames(., sub('.sf', '', basename(.)))
txi <- tximport(fs, 'salmon', tx2gene = tx2gene)

md <- tibble(samp = colnames(txi$abundance)) %>%
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
  na.omit() %>%
  column_to_rownames('samp')

txii <- lapply(txi[1:3], `[`, , rownames(md))
txii$countsFromAbundance <- txi$countsFromAbundance


dds <- DESeqDataSetFromTximport(txii, md, ~media) %>% estimateSizeFactors()

tclrs <- c('#e45756', '#4c78a8')
signif.num <- function(x) {
  as.character(symnum(x, corr = FALSE, na = FALSE, legend = FALSE,
                      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                      symbols = c("****", "***", "**", "*", "ns")))
}

res <- md %>% 
  split(., .$line) %>% 
  lapply(function(x) {
    print(x$line[1])
    x %>%
      split(.,.$geno) %>%
      lapply(function(y) {
        print(y$geno[1])
        dds[,rownames(y)] %>%
          DESeq() %>%
          results(c('media', 'DM', 'SCM')) %>%
          as.data.frame() %>% 
          rownames_to_column('gene_id')
      })
  })


dds2 <- DESeqDataSetFromTximport(txii, md, ~geno) %>% estimateSizeFactors()


res2 <- md %>% 
  split(., .$line) %>% 
  lapply(function(x) {
    print(x$line[1])
    x %>%
      split(.,.$media) %>%
      lapply(function(y) {
        print(y$media[1])
        dds2[,rownames(y)] %>%
          DESeq() %>%
          results(c('geno',  'KO', 'K27M')) %>%
          as.data.frame() %>% 
          rownames_to_column('gene_id')
      })
  })

list(media = res2, geno = res) %>%
  lapply(function(x) {
    lapply(x, bind_rows, .id = 'grp') %>%
      bind_rows(.id = 'line')
  }) %>%
  bind_rows(.id = 'kind') -> rres


cclrs <- kelly()[c(3,7,6,2)]
readRDS('../data/pmtcgi/GBM.embed.rds') %>%
  mutate(cl = case_when(clust == '0' ~ '1',
                        clust == '1_1' ~ '4',
                        clust %in% c('1_0_0', '1_0_3_0') ~ '2b',
                        clust %in% c('1_0_1', '1_0_3_1') ~ '2a',
                        clust %in% c('1_0_2_0_0', '1_0_2_1_1') ~ '3b',
                        clust %in% c('1_0_2_0_1', '1_0_2_1_0') ~ '3a',
                        T ~ NA_character_)) %>%
  na.omit() %>%
  merge(mutate(fread('../data/pmtcgi/hg38.5kb.bed'), idx = 1:n()), by = 'idx') %>%
  dplyr::select(gene_id = V4, line, cl) %>%
  merge(rres, by = c('line', 'gene_id')) %>%
  filter(kind == 'media') %>%
  filter(line != 'BT245') %>%
  na.omit() %>%
  arrange(cl) %>%
  mutate(clr = kelly()[c(3,7,6,2)][as.integer(sub('[ab]$', '', cl))],
         clu = c('1' = 'Active',
                 '2a' = 'cPRC1, H3K4me3-',
                 '2b' = 'cPRC1, H3K4me3+',
                 '3a' = 'PRC2, H3K4me3-',
                 '3b' = 'PRC2, H3K4me3+',
                 '4' = 'Other')[cl],
         media = c('DM' = 'Differentiation media', 'SCM'  = 'Stem cell media')[grp] %>%
           factor(c('Stem cell media', 'Differentiation media')),
         clu = sprintf("<span style='color:%s'>%s</span>", clr, clu) %>% fct_inorder()) %>%
  ggplot(aes(x = clu, y= log2FoldChange, fill = clr, color = clr)) +
  scale_x_discrete() +
  annotate('rect', xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf, alpha = .1,
           fill = 'chartreuse4') +
  geom_hline(yintercept = 0) +
  geom_violin(alpha = .3, color = NA, scale = 'width') +
  geom_boxplot(outlier.color = NA, width = .2) +
  stat_summary(geom = "crossbar", width=0.1, fatten=0, color="white",
               fun.data = function(x) c(y=median(x), ymin=median(x), ymax=median(x))) +
  scale_color_identity() +
  scale_fill_identity() +
  facet_nested(line~media) +
  coord_cartesian(ylim = c(-10, 10)) +
  labs(y = expression(log[2]*' '*frac('K27M-KO', 'K27M')*' expression'),
       x = 'Gene class') +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black'),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = 'grey80', linetype = 'dashed'),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(size = 11, color = 'white'),
        axis.text = element_text(color = 'black', size = 11),
        axis.text.x = element_markdown(angle = 45, hjust = 1),
        axis.title.x = element_blank()) -> p
ggsave('sf11_a.pdf', p, height = 8, width = 6)




fss <- list.files('../data/rna/CBXAM', full.names = T) %>%
  setNames(., sub('.sf', '', basename(.)))

mdatt <-  data.frame(s = names(fss)) %>%
  mutate(rep = ifelse(grepl('R2$', s), 2, 1),
         line = case_when(grepl('BT', s) ~ 'BT245',
                          grepl('HSJ', s) ~ 'HSJ019',
                          T ~ 'DIPGXIII'),
         geno = ifelse(grepl('c5|c8', s), 'K27M-KO', 'K27M'),
         cond = case_when(grepl('mix', s) ~ 'mix',
                          grepl('SW2', s) ~ 'SW2',
                          grepl('3866', s) ~ '3866',
                          grepl('4976', s) ~ '4976',
                          T ~ 'ctrl'),
         grp = case_when(geno == 'K27M' & cond == '4976' ~ 'K27M 4976',
                         geno == 'K27M' ~ 'K27M ctrl',
                         geno == 'K27M-KO' ~ 'K27M-KO ctrl',
                         T ~ 'Other') %>%
           factor(c('K27M ctrl', 'K27M-KO ctrl', 'K27M 4976', 'Other'))) %>%
  column_to_rownames('s')


salmonn <- tximport(fss, type = "salmon", tx2gene = tx2gene)
salmonn$length[salmonn$length == 0] <- 1

ddss <- DESeqDataSetFromTximport(txi = salmonn,
                                 colData = mdatt,
                                 design = ~grp)

ress <- unique(mdatt$line) %>%
  setNames(.,.) %>%
  lapply(function(x) {
    dds2 <- DESeq(ddss[, ddss$line == x])
    c('K27M 4976', 'K27M-KO ctrl') %>%
      setNames(.,.) %>%
      lapply(function(y) {
        results(dds2, c('grp', y, 'K27M ctrl')) %>%
          as.data.frame() %>%
          rownames_to_column('gene_id')
      })
  })

cls <- readRDS('../data/pmtcgi/GBM.embed.rds') %>%
  mutate(cl = case_when(clust == '0' ~ '1',
                        clust == '1_1' ~ '4',
                        clust %in% c('1_0_0', '1_0_3_0') ~ '2b',
                        clust %in% c('1_0_1', '1_0_3_1') ~ '2a',
                        clust %in% c('1_0_2_0_0', '1_0_2_1_1') ~ '3b',
                        clust %in% c('1_0_2_0_1', '1_0_2_1_0') ~ '3a',
                        T ~ NA_character_)) %>%
  na.omit() %>%
  merge(mutate(fread('../data/pmtcgi/hg38.5kb.bed'), idx = 1:n()), by = 'idx') %>%
  dplyr::select(gene_id = V4, line, cl) %>%
  filter(grepl('ENSG', gene_id)) %>%
  split(., .$line)

c('DIPGXIII', 'HSJ019') %>%
  lapply(function(s) {
    ress[[s]] %>%
      rbindlist(idcol = 'cmp') %>%
      na.omit() %>%
      dplyr::select(gene_id, cmp, log2FoldChange) %>%
      merge(cls[[s]], by = 'gene_id') %>%
      mutate(line = s)
  }) %>%
  bind_rows() %>%
  arrange(cl) %>%
  mutate(clr = kelly()[c(3,7,6,2)][as.integer(sub('[ab]$', '', cl))],
         clu = c('1' = 'Active',
                 '2a' = 'cPRC1, H3K4me3-',
                 '2b' = 'cPRC1, H3K4me3+',
                 '3a' = 'PRC2, H3K4me3-',
                 '3b' = 'PRC2, H3K4me3+',
                 '4' = 'Other')[cl],
         clu = sprintf("<span style='color:%s'>%s</span>", clr, clu) %>% fct_inorder()) %>%
  filter(cmp == 'K27M 4976') %>%
  ggplot(aes(x = clu, y = log2FoldChange, fill = clr, color = clr)) +
  scale_x_discrete() +
  annotate('rect', xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf, alpha = .1,
           fill = 'chartreuse4') +
  geom_hline(yintercept = 0) +
  geom_violin(alpha = .3, color = NA, scale = 'width') +
  geom_boxplot(outlier.color = NA, width = .2) +
  stat_summary(geom = "crossbar", width=0.1, fatten=0, color="white",
               fun.data = function(x) c(y=median(x), ymin=median(x), ymax=median(x))) +
  scale_color_identity() +
  scale_fill_identity() +
  labs(y = expression(log[2]*' '*frac('K27M + CBX-AM', 'K27M')*' expression'),
       x = 'Gene class') +
  facet_grid(line~'Differentiation media') +
  coord_cartesian(ylim = c(-10, 10)) +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black'),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = 'grey80', linetype = 'dashed'),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(size = 11, color = 'white'),
        axis.text = element_text(color = 'black', size = 11),
        axis.text.x = element_markdown(angle = 45, hjust = 1),
        axis.title.x = element_blank()) -> p
ggsave('sf11_b.pdf', p, height = 8, width = 3)

rr2 <- ress %>%
  Map(function(x, l) {
    c1 <- 'K27M 4976'
    c2 <- 'K27M-KO ctrl'
    d <- x[[c1]] %>%
      dplyr::select(gene_id, p = pvalue, lfc = log2FoldChange) %>%
      merge(x[[c2]][,c('gene_id','pvalue','log2FoldChange')], by = 'gene_id') %>%
      na.omit() %>%
      mutate(s1 = -log10(p) * ifelse(lfc > 0, 1, -1),
             s2 = -log10(pvalue) * ifelse(log2FoldChange > 0, 1, -1))
    
    cls[[l]] %>%
      filter(grepl('[ab]$', cl)) %>%
      mutate(cl = sub('[ab]$', '', cl)) %>%
      rbind(cls[[l]]) %>%
      split(., .$cl) %>%
      lapply(function(y) {
        dd <- d[d$gene_id %in% y$gene_id,]
        ro <- RRHO2_initialize(dd[,c('gene_id', 's1')], dd[,c('gene_id', 's2')],
                               labels = c(c1, c2), log10.ind = T)
      })
  }, ., names(.))

rr2 %>%
  lapply(function(x) {
    x %>%
      lapply(function(y) {
        y$hypermat %>%
          as.data.frame() %>%
          rownames_to_column('x') %>%
          pivot_longer(-x, names_to = 'y', values_to = 'z') %>%
          mutate(x = as.numeric(x),
                 y = as.numeric(sub('^V', '', y)))
      }) %>% 
      bind_rows(.id = 'cl')
  }) %>%
  bind_rows(.id = 'line') %>%
  mutate(cl = c('1' = 'Active','2' = 'cPRC1', '3' = 'PRC2', '4' = 'Other')[cl]) %>%
  filter(cl == 'cPRC1') %>%
  ggplot(aes(x = x, y = y, fill = z)) +
  geom_tile() +
  scale_fill_viridis_c(name = expression(-log[10]~p[adj]), 
                       na.value = '#00000000',
                       breaks = scales::pretty_breaks(3),
                       guide = guide_colorbar(title.vjust = 1,
                                              barheight = .7,
                                              barwidth = 5)) +
  coord_cartesian(expand = F) +
  labs(x = expression('Rank by '*frac('K27M + CBX-AM','K27M')),
       y = expression('Rank by '*frac('K27M-KO','K27M'))) +
  facet_grid(cl~line, scales = 'free') +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        strip.background = element_rect(fill = 'black'),
        legend.position = 'bottom',
        legend.justification = 'right',
        strip.text = element_text(size = 13, color = 'white', face = 'bold'),
        legend.key = element_blank()) -> p
ggsave('sf11_e.pdf', p, height = 3.7, width = 6.9, bg = 'transparent')