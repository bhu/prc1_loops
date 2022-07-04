library(data.table)
library(tidyverse)
library(rtracklayer)
library(InteractionSet)
library(colorspace)
library(RcppCNPy)
library(pals)
library(ggh4x)
library(patchwork)
library(gggenes)
library(DESeq2)
library(tximport)
library(ggsignif)
library(ggtext)
library(gghalves)
library(ggpubr)


fall <- matrix(c(255, 255, 255,
                 255, 255, 204,
                 255, 237, 160,
                 254, 217, 118,
                 254, 178, 76,
                 253, 141, 60,
                 252, 78, 42,
                 227, 26, 28,
                 189, 0, 38,
                 128, 0, 38,
                 0, 0, 0), nrow = 3) %>%
  apply(2, function(a) rgb(a[1], a[2], a[3], maxColorValue = 255))

rot <- function(df, degree) {
  dfr <- df
  degree <- pi * degree / 180
  l <- sqrt(df$x^2 + df$y^2)
  teta <- atan(df$y / df$x)
  dfr$x <- round(l * cos(teta - degree))
  dfr$y <- round(l * sin(teta - degree))
  return(dfr)
}

reg <- GRanges('chr6', IRanges(99e6, 101e6))
reg2 <- GRanges('chr6', IRanges(99563000, 99618000))
xbrks <- seq(99e6, 101e6, 1e6)
xlbs <- paste0(99:101, 'mb')

tclrs <- c('#e45756', '#4c78a8')
bclr <- '#54a24b'

res <- 10000

mats <- c('K27M' = 'BT245_K27M', 'K27M-KO' = 'BT245_KO') %>%
  lapply(function(samp) {
    crd <- sprintf('../data/hic/%s/Full_Mats_Coords/coords_%s_res_%d.npy',
                   samp, seqnames(reg), res) %>%
      npyLoad('integer')
    r <- (which(crd == start(reg) / res) - 100):(which(crd == end(reg) / res) + 100)
    sprintf('../data/hic/%s/Full_Enhanced/%s_res_%d.npy',
            samp, seqnames(reg), res) %>%
      npyLoad() %>% 
      .[r, r] %>%
      as.data.frame() %>% 
      mutate(y = 1:n()) %>% 
      pivot_longer(-y, names_to = 'x', values_to = 'z') %>% 
      mutate(x = sub('^V', '', x) %>% as.integer())
  })

p1 <- mats %>%
  bind_rows(.id = 'samp') %>%
  group_by(samp) %>%
  mutate(across(c(x, y), ~.x/sqrt(2))) %>%
  do(., rot(., 45)) %>%
  ungroup() %>%
  filter(between(y, 0, 50)) %>%
  filter(between(x, 101,300)) %>%
  ggplot(aes(x, y, fill = z)) +
  geom_raster() +
  facet_wrap(~samp) +
  coord_cartesian(expand = F) +
  scale_fill_gradientn('Interaction',
                       colours = fall, na.value=fall[1],
                       trans = 'log', limits = c(.05, 1), 
                       oob = scales::squish, breaks = c(0.1,.5,1),
                       guide = guide_colorbar(barwidth = .5, barheight = 3.8,
                                              title.hjust = .5,
                                              title.vjust = 0.5,
                                              title.position = 'left')) +
  facet_wrap2(~samp,
              strip = strip_themed(background_x = elem_list_rect(fill = tclrs)))+
  scale_y_continuous(breaks = c(25, 50),
                     labels = c('500kb','1mb')) +
  scale_x_continuous(breaks = c(101, 200, 300),
                     labels = xlbs) +
  labs(x = seqnames(reg), y = 'Separation') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        strip.background = element_rect(fill = 'black'),
        axis.text.x = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_text(angle = 90),
        axis.title.x = element_blank(),
        strip.text = element_text(color = 'white', size = 11))

bws <- c('K27M' = 'BT245_K27M_H3K27me3.bw',
         'K27M-KO' = 'BT245_KO_H3K27me3.bw') %>%
  lapply(function(x) {
    import.bw(sprintf('../data/bws/%s', x),
              selection = BigWigSelection(reg))
  })

bins <- tile(reg, 500)[[1]]
seqlevels(bins) <- seqlevels(bws[[1]])
xs <- mid(bins)

rx <- readRDS('../../rx.rds') %>%
  {setNames(.$r, .$samp)} %>%
  {list('K27M' = .['BT245-C24-Rx-C_H3K27me3'] / .['BT245_C24-Rx_input'],
        'K27M-KO' = .['BT245-ko-c2p3-Rx-C_H3K27me3'] / .['BT245-ko-c2p3-Rx_input'])} %>%
  lapply(unname) %>%
  unlist() %>%
  {./mean(.)}

scs <- lapply(bws, function(x) {
  mcolAsRleList(x, 'score') %>%
    binnedAverage(bins, ., 'score') %>%
    score() %>%
    tibble(score = .) %>%
    mutate(idx = xs)
}) %>%
  bind_rows(.id = 'samp') %>%
  na.omit() %>%
  mutate(score = score * rx[samp])

p2 <- scs %>% 
  ggplot(aes(x = idx, y = score, fill = samp, color = samp)) +
  #geom_rect(aes(xmin = xmin, xmax = xmax), ymin = -Inf, ymax = Inf,
  #          data = tibble(samp = names(bws), xmin = start(reg2), xmax = end(reg2)),
  #          fill = '#000000', color = NA, alpha = .15, inherit.aes = F) +
  geom_area(alpha = .75) +
  geom_line() +
  coord_cartesian(xlim = c(start(reg), end(reg)), ylim = c(0, .7)) +
  scale_y_continuous('H3K27me3\n(ChIP-Rx)',
                     breaks = c(0, .5)) +
  scale_x_continuous(expand = expansion(0)) +
  facet_wrap(~samp) +
  scale_fill_manual(values = tclrs) +
  scale_color_manual(values = tclrs) +
  #annotate('label', x = Inf, y = Inf, hjust = 1, vjust = 1, 
  #         label = 'K27M', color = tableau20()[1]) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        legend.position = 'none',
        plot.margin = margin(),
        axis.line.y = element_line(color = 'black'))

p3 <- c('K27M' = 'BT245_K27M_CBX2.bw',
        'K27M-KO' = 'BT245_KO_CBX2.bw') %>%
  lapply(function(x) {
    import.bw(sprintf('../data/bws/%s', x),
              selection = BigWigSelection(reg)) %>%
      mcolAsRleList('score') %>%
      binnedAverage(bins, ., 'score') %>%
      score() %>%
      tibble(score = .) %>%
      mutate(idx = xs)
  }) %>%
  bind_rows(.id = 'samp') %>%
  na.omit() %>%
  ggplot(aes(x = idx, y = score, fill = samp, color = samp)) +
  #geom_rect(aes(xmin = xmin, xmax = xmax), ymin = -Inf, ymax = Inf,
  #          data = tibble(samp = names(bws), xmin = start(reg2), xmax = end(reg2)),
  #          fill = '#000000', color = NA, alpha = .15, inherit.aes = F) +
  geom_area(alpha = .75) +
  geom_line() +
  coord_cartesian(xlim = c(start(reg), end(reg))) +
  scale_y_continuous('CBX2\n(CPM)',
                     breaks = c(0, 2)) +
  scale_x_continuous(expand = expansion(0)) +
  facet_wrap(~samp) +
  scale_fill_manual(values = tclrs) +
  scale_color_manual(values = tclrs) +
  #annotate('label', x = Inf, y = Inf, hjust = 1, vjust = 1, 
  #         label = 'K27M', color = tableau20()[1]) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        legend.position = 'none',
        plot.margin = margin(),
        axis.line.y = element_line(color = 'black'))


p4 <- c('K27M' = 'BT245_K27M_H3K4me3.bw',
        'K27M-KO' = 'BT245_KO_H3K4me3.bw') %>%
  lapply(function(x) {
    import.bw(sprintf('../data/bws/%s', x),
              selection = BigWigSelection(reg)) %>%
      mcolAsRleList('score') %>%
      binnedAverage(bins, ., 'score') %>%
      score() %>%
      tibble(score = .) %>%
      mutate(idx = xs)
  }) %>%
  bind_rows(.id = 'samp') %>%
  na.omit() %>%
  ggplot(aes(x = idx, y = score, fill = samp, color = samp)) +
  #geom_rect(aes(xmin = xmin, xmax = xmax), ymin = -Inf, ymax = Inf,
  #          data = tibble(samp = names(bws), xmin = start(reg2), xmax = end(reg2)),
  #          fill = '#000000', color = NA, alpha = .15, inherit.aes = F) +
  geom_area(alpha = .75) +
  geom_line() +
  coord_cartesian(xlim = c(start(reg), end(reg))) +
  scale_y_continuous('H3K4me3\n(CPM)',
                     breaks = c(0, 2)) +
  scale_x_continuous(expand = expansion(0)) +
  facet_wrap(~samp) +
  scale_fill_manual(values = tclrs) +
  scale_color_manual(values = tclrs) +
  #annotate('label', x = Inf, y = Inf, hjust = 1, vjust = 1, 
  #         label = 'K27M', color = tableau20()[1]) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        legend.position = 'none',
        plot.margin = margin(),
        axis.line.y = element_line(color = 'black'))



gclrs <- c('#f28e2b', '#b07aa1')


txs <- list(mm10 = 'GCF_000001635.26_GRCm38.p6_genomic.gff.gz',
            hg38 = 'GCF_000001405.40_GRCh38.p14_genomic.gff.gz') %>%
  lapply(function(x) {
    tx <- import(file.path('..', 'data', 'genes', x))
    rnm <- sub('_genomic.gff.gz', '_assembly_report.txt', x) %>%
      file.path('..', 'data', 'genes', .) %>%
      read_delim(comment = '#', delim = '\t', col_names = F) %>%
      filter(X2 == 'assembled-molecule') %>%
      dplyr::select(X7, X10) %>%
      deframe()
    tx[!is.na(tx$tag) & seqnames(tx) %in% names(rnm)] %>%
      keepSeqlevels(names(rnm)) %>%
      renameSeqlevels(rnm)
  })


subsetByOverlaps(txs$hg38, reg) %>%
  as_tibble() %>%
  filter(type %in% c('exon',  'mRNA')) %>%
  dplyr::select(strand, start, end, type, gene) %>%
  rbind(tibble(strand = c('+', '+', '-', '-') %>% fct_inorder(),
               start = rep(-1, 4),
               end = rep(-1, 4),
               gene = rep(c('a','b'), each = 2),
               type = rep(c('mRNA', 'exon'), 2))) %>%
  mutate(forward = strand == '+',
         position = case_when(strand == '-' ~ end, T ~ start)) %>%
  {rbind(mutate(., s = 'K27M'), mutate(., s = 'K27M-KO'))} %>%
  ggplot(aes(xmin = start, xmax = end, y = strand, color = strand,
             fill = strand, forward = forward)) +
  geom_rect(aes(xmin = xmin, xmax = xmax), ymin = -Inf, ymax = Inf,
            data = tibble(s = names(bws), xmin = start(reg2), xmax = end(reg2)),
            fill = '#000000', color = NA, alpha = .15, inherit.aes = F) +
  geom_gene_arrow(data = ~subset(., type == 'mRNA')) +
  scale_alpha_manual(values = c('UTR' = .5, 'not' = 1)) +
  scale_x_continuous(seqnames(reg), limits =  c(start(reg), end(reg)),
                     breaks = xbrks,
                     labels = xlbs,
                     expand = expansion(0)) +
  scale_fill_manual(values = gclrs) +
  scale_color_manual(values = gclrs) +
  facet_wrap(~s) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        axis.line.x = element_line(color = 'black'),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = 'none') -> p5



bws2 <- c('K27M' = 'BT245_K27M_H3K4me3.bw',
          'K27M-KO' = 'BT245_KO_H3K4me3.bw') %>%
  lapply(function(x) {
    import.bw(sprintf('../data/bws/%s', x),
              selection = BigWigSelection(reg2))
  })

bins2 <- tile(reg2, 500)[[1]]
seqlevels(bins2) <- seqlevels(bws2[[1]])
xs2 <- mid(bins2)

scs2 <- lapply(bws2, function(x) {
  mcolAsRleList(x, 'score') %>%
    binnedAverage(bins2, ., 'score') %>%
    score() %>%
    tibble(score = .) %>%
    mutate(idx = xs2)
}) %>%
  bind_rows(.id = 'samp') 



p6 <- scs2 %>%
  na.omit() %>%
  ggplot(aes(x = idx, y = score, fill = samp, color = samp)) +
  geom_rect(aes(xmin = xmin, xmax = xmax), ymin = -Inf, ymax = Inf,
            data = tibble(samp = names(bws2), xmin = 99601614, xmax = 99605935),
            fill = '#BCBD22', color = NA, alpha = .2, inherit.aes = F) +
  geom_area(alpha = .75) +
  geom_line() +
  coord_cartesian(xlim = c(start(reg2), end(reg2))) +
  scale_y_continuous('H3K4me3\n(CPM)',
                     breaks = c(0, 7)) +
  scale_x_continuous(expand = expansion(0)) +
  facet_wrap(~samp) +
  scale_fill_manual(values = tclrs) +
  scale_color_manual(values = tclrs) +
  #annotate('label', x = Inf, y = Inf, hjust = 1, vjust = 1, 
  #         label = 'K27M', color = tableau20()[1]) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        legend.position = 'none',
        plot.margin = margin(),
        axis.line.y = element_line(color = 'black'))


ggplot(sg, aes(xmin = start, xmax = end, y = strand, color = strand,
               fill = strand, forward = forward)) +
  #geom_rect(aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = 'gold2',
  #          color = NA, data = distinct(rcts, xmin, xmax), inherit.aes = F, alpha = .4) +
  #geom_linerange(data = ~subset(., type == 'transcript')) 
  geom_gene_arrow(data = ~subset(., type == 'transcript')) +
  #geom_gene_arrow(aes(alpha = utr), color = NA, data = ~subset(., type != 'transcript')) +
  scale_alpha_manual(values = c('UTR' = .5, 'not' = 1)) +
  geom_text(aes(x = end, label = gene_name), hjust = -.1,
            data = ~subset(., type == 'transcript' & strand == '-')) +
  geom_text(aes(x = start, label = gene_name), hjust = 1.1,
            data = ~subset(., type == 'transcript' & strand == '+')) +
  #geom_feature(aes(x = position, y = strand, forward = forward, color = strand),
  #             data = ~subset(., type == 'transcript')) +
  coord_cartesian(xlim = c(start(reg2), end(reg2))) +
  scale_x_continuous(seqnames(reg), 
                     breaks = c(99.58e6, 99.6e6),
                     labels = c('99.58mb', '99.6mb'),
                     expand = expansion(0)) +
  scale_fill_manual(values = gclrs) +
  scale_color_manual(values = gclrs) +
  facet_wrap(~s) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        axis.line.x = element_line(color = 'black'),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.margin = margin(),
        legend.position = 'none') -> p7


subsetByOverlaps(txs$hg38, reg2) %>%
  as_tibble() %>%
  filter(type %in% c('exon',  'mRNA')) %>%
  dplyr::select(strand, start, end, type, gene) %>%
  rbind(tibble(strand = c('+', '+', '-', '-') %>% fct_inorder(),
               start = rep(-1, 4),
               end = rep(-1, 4),
               gene = rep(c('a','b'), each = 2),
               type = rep(c('mRNA', 'exon'), 2))) %>%
  mutate(forward = strand == '+',
         position = case_when(strand == '-' ~ end, T ~ start)) %>%
  {rbind(mutate(., s = 'K27M'), mutate(., s = 'K27M-KO'))} %>%
  ggplot(aes(xmin = start, xmax = end, y = strand, color = strand,
             fill = strand, forward = forward)) +
  geom_gene_arrow(data = ~subset(., type == 'mRNA')) +
  scale_alpha_manual(values = c('UTR' = .5, 'not' = 1)) +
  scale_x_continuous(seqnames(reg), breaks = c(99.58e6, 99.6e6),
                     labels = c('99.58mb', '99.6mb'),
                     expand = expansion(0)) +
  scale_fill_manual(values = gclrs) +
  scale_color_manual(values = gclrs) +
  facet_wrap(~s) +
  coord_cartesian(xlim = c(start(reg2), end(reg2))) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        axis.line.x = element_line(color = 'black'),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = 'none') -> p7

{ wrap_plots(p1,p2,p3,p5,p6,p7,ncol=1, heights = c(1.2,1,1,.7,1,.7)) &
    theme(plot.background = element_blank()) } %>%
  ggsave('f6_a.pdf', ., width = 7.4, height = 5.3, bg = 'transparent')



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

gene_name <- 'PRDM13'
txii$abundance[names(e2g[e2g == gene_name]),] %>% 
  data.frame(v = .) %>% rownames_to_column('samp') %>%
  merge(rownames_to_column(md, 'samp')) %>%
  filter(line == 'BT245') %>%
  filter(media == 'SCM') %>% 
  filter(!grepl('D', samp)) %>%
  mutate(gene_name = gene_name) %>%
  ggplot(aes(x = geno, y = v, color = geno, fill = geno)) + 
  geom_boxplot(position = position_nudge(x = .3), width = .2,
               outlier.color = NA) +
  geom_jitter(width = .1, height = 0) +
  geom_signif(comparisons = list(c('K27M','KO')),
              test = 't.test', color = 'black', tip_length = 0,
              map_signif_level = signif.num) +
  scale_color_manual(values = tclrs) +
  scale_fill_manual(values = tclrs) +
  stat_summary(geom = "crossbar", width=0.15, fatten=0, color="white", 
               fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) },
               position = position_nudge(x = .3)) +
  ylab('Expression (TPM)') +
  scale_y_continuous(expand = expansion(c(0.03, .1))) +
  facet_wrap(~gene_name) +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black'),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = 'grey80', linetype = 'dashed'),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(size = 11, color = 'white'),
        axis.text = element_text(color = 'black', size = 11),
        legend.position = 'none',
        axis.title.x = element_blank()) -> p
ggsave('f6_b.pdf', p, width = 3, height = 2.5, bg = 'transparent')






















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
  dplyr::select(line, cl, idx) %>%
  merge(mutate(fread('../data/pmtcgi/hg38.5kb.bed'), idx = 1:n())) %>%
  dplyr::select(line, cl, gene_id = V4) %>%
  merge(rres, by = c('line', 'gene_id')) %>%
  filter(line == 'BT245') %>%
  filter(kind == 'media') %>%
  na.omit() %>%
  arrange(cl) %>%
  mutate(clr = cclrs[as.integer(sub('[ab]$', '', cl))],
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
  geom_hline(yintercept = 0) +
  geom_violin(alpha = .3, color = NA, scale = 'width') +
  geom_boxplot(outlier.color = NA, width = .2) +
  stat_summary(geom = "crossbar", width=0.1, fatten=0, color="white",
               fun.data = function(x) c(y=median(x), ymin=median(x), ymax=median(x))) +
  scale_color_identity() +
  scale_fill_identity() +
  facet_wrap(~media) +
  coord_cartesian(ylim = c(-6, 6)) +
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
        axis.title.x = element_blank()) -> plt1








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




s <- 'BT245'
ress[[s]] %>%
  rbindlist(idcol = 'cmp') %>%
  na.omit() %>%
  dplyr::select(gene_id, cmp, log2FoldChange) %>%
  merge(cls[[s]], by = 'gene_id') %>%
  arrange(cl) %>%
  mutate(clr = cclrs[as.integer(sub('[ab]$', '', cl))],
         clu = c('1' = 'Active',
                 '2a' = 'cPRC1, H3K4me3-',
                 '2b' = 'cPRC1, H3K4me3+',
                 '3a' = 'PRC2, H3K4me3-',
                 '3b' = 'PRC2, H3K4me3+',
                 '4' = 'Other')[cl],
         clu = sprintf("<span style='color:%s'>%s</span>", clr, clu) %>% fct_inorder()) %>%
  filter(cmp == 'K27M 4976') %>%
  ggplot(aes(x = clu, y = log2FoldChange, fill = clr, color = clr)) +
  geom_hline(yintercept = 0) +
  geom_violin(alpha = .3, color = NA, scale = 'width') +
  geom_boxplot(outlier.color = NA, width = .2) +
  stat_summary(geom = "crossbar", width=0.1, fatten=0, color="white",
               fun.data = function(x) c(y=median(x), ymin=median(x), ymax=median(x))) +
  scale_color_identity() +
  scale_fill_identity() +
  labs(y = expression(log[2]*' '*frac('K27M + CBX-AM', 'K27M')*' expression'),
       x = 'Gene class') +
  facet_wrap(~'Differentiation media') +
  coord_cartesian(ylim = c(-6, 6)) +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black'),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = 'grey80', linetype = 'dashed'),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(size = 11, color = 'white'),
        axis.text = element_text(color = 'black', size = 11),
        axis.text.x = element_markdown(angle = 45, hjust = 1),
        axis.title.x = element_blank()) -> plt2

{wrap_plots(plt1, plt2, nrow = 1, widths = c(2,1)) &
    theme(plot.background = element_blank())} %>%
  ggsave('f6_cd.pdf', .)








readRDS('../data/pmtcgi/GBM.embed.rds') %>%
  mutate(cl = case_when(clust == '0' ~ '1',
                        clust == '1_1' ~ '4',
                        clust %in% c('1_0_0', '1_0_3_0') ~ '2b',
                        clust %in% c('1_0_1', '1_0_3_1') ~ '2a',
                        clust %in% c('1_0_2_0_0', '1_0_2_1_1') ~ '3b',
                        clust %in% c('1_0_2_0_1', '1_0_2_1_0') ~ '3a',
                        T ~ NA_character_))  %>%
  dplyr::select(line, cl, idx) %>%
  na.omit() %>%
  merge(mutate(fread('../data/pmtcgi/hg38.5kb.bed'), idx = 1:n())) %>%
  mutate(kind = ifelse(grepl('ENSG', V4), 'TSS', 'CGI')) %>%
  dplyr::select(line, chr = V1, start = V2, end = V3, cl, kind) %>%
  {rbind(., mutate(filter(., grepl('[ab]$', cl)), cl = sub('[ab]$', '', cl)))} %>%
  {rbind(., mutate(., kind = 'both'))} -> d1

readRDS('../data/pmtcgi/GBM.embed.rds') %>%
  mutate(cl = case_when(clust == '0' ~ '1',
                        clust == '1_1' ~ '4',
                        clust %in% c('1_0_0', '1_0_3_0') ~ '2b',
                        clust %in% c('1_0_1', '1_0_3_1') ~ '2a',
                        clust %in% c('1_0_2_0_0', '1_0_2_1_1') ~ '3b',
                        clust %in% c('1_0_2_0_1', '1_0_2_1_0') ~ '3a',
                        T ~ NA_character_))  %>%
  dplyr::select(line, cl, idx) %>%
  na.omit() %>%
  {rbind(., mutate(filter(., grepl('[ab]$', cl)), cl = sub('[ab]$', '', cl)))} %>%
  dplyr::count(cl, idx) %>%
  filter(n == 3) %>%
  merge(mutate(fread('../data/pmtcgi/hg38.5kb.bed'), idx = 1:n())) %>%
  mutate(kind = ifelse(grepl('ENSG', V4), 'TSS', 'CGI'),
         line = 'int') %>%
  dplyr::select(line, chr = V1, start = V2, end = V3, cl, kind) %>%
  {rbind(., mutate(., kind = 'both'))} -> d2


regs <- readRDS('../data/pmtcgi/GBM.embed.rds') %>%
  mutate(cl = case_when(clust == '0' ~ '1',
                        clust == '1_1' ~ '4',
                        clust %in% c('1_0_0', '1_0_3_0') ~ '2b',
                        clust %in% c('1_0_1', '1_0_3_1') ~ '2a',
                        clust %in% c('1_0_2_0_0', '1_0_2_1_1') ~ '3b',
                        clust %in% c('1_0_2_0_1', '1_0_2_1_0') ~ '3a',
                        T ~ NA_character_))  %>%
  dplyr::select(line, cl, idx) %>%
  na.omit() %>%
  merge(mutate(fread('../data/pmtcgi/hg38.5kb.bed'), idx = 1:n())) %>%
  mutate(kind = ifelse(grepl('ENSG', V4), 'TSS', 'CGI')) %>%
  {rbind(., mutate(., kind = 'both'))} %>%
  mutate(line = 'all',
         cl = 'all') %>%
  dplyr::select(line, chr = V1, start = V2, end = V3, cl, kind) %>%
  rbind(d1, d2) %>%
  filter(!(line %in% c('DIPGXIII', 'HSJ019'))) %>%
  mutate(grp = paste(line, cl, kind)) %>%
  split(., .$grp)



list.files('../data/aggr', full.names = T) %>%
  grep('DMSO|UNC', ., value = T) %>%
  setNames(., sub('.mat.gz', '', basename(.))) %>%
  lapply(function(x) {
    r <- fread(x, skip = 1, select = 1:3, col.names = c('chr','start','end')) %>%
      mutate(idx = 1:n())
    m <- fread(x, skip = 1, drop = 1:6) %>% as.matrix()
    lapply(regs, function(y) {
      merge(r, y, by = c('chr', 'start', 'end')) %>%
        pull(idx) %>%
        m[.,] %>%
        {tibble(mu = colMeans(., na.rm = T),
                med = colMedians(., na.rm = T),
                std = colSds(., na.rm = T),
                num = colSums(is.finite(.)))} %>%
        mutate(idx = 1:n()) 
    }) %>% bind_rows(.id = 'grp')
  }) %>% bind_rows(.id = 'f') -> agg


clrs2 <- c('#ff9da7', '#76b7b2')

agg %>%
  mutate(cond = ifelse(grepl('DMSO', f), 'DMSO', 'UNC4976'),
         mark = case_when(grepl('CBX2', f) ~ 'CBX2',
                          grepl('RING1B', f) ~ 'RING1B',
                          T ~ 'H3K27me3')) %>%
  filter(mark != 'H3K27me3') %>%
  separate(grp, c('line', 'cl', 'kind')) %>%
  filter(line == 'BT245' & kind == 'CGI' &
           !grepl('a|b', cl)) %>%
  arrange(cl) %>%
  mutate(cl = c('1' = 'Active', '2' = 'cPRC1',
                '3' = 'PRC2', '4' = 'Other')[cl] %>%
           fct_inorder(),
         cond = sub('UNC4976', 'CBX-AM', cond) %>%
           factor(c('DMSO','CBX-AM')),
         se = std / sqrt(num)) %>%
  ggplot(aes(x = idx, y = mu, color = cond, fill = cond,
             ymin = mu - se, ymax = mu + se)) +
  geom_ribbon(color = NA, alpha = .5) +
  geom_line() +
  facet_nested(cl ~ mark, scales = 'free', independent = 'y',
               strip = strip_nested(background_y = elem_list_rect(
                 fill = pals::kelly()[c(3,7,6,2)]
               ))) +
  scale_x_continuous(limits = c(950, 1050),
                     breaks = c(950, 1000, 1050),
                     labels = c('-50kb', 'CGI', '+50kb')) +
  ylab('CPM') +
  scale_color_manual(values = clrs2) +
  scale_fill_manual(values = clrs2) +
  facetted_pos_scales(y = list(
    mark == 'CBX2' ~ scale_y_continuous(limits = c(0, .8), breaks = c(0, .3, .6)),
    mark == 'RING1B' ~ scale_y_continuous(limits = c(0.03, .11), breaks = c(.03, .06, .09))
  )) +
  theme(axis.text = element_text(size = 11, color = 'black'),
        plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(size = 11, color = 'white'),
        axis.title.x = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.text = element_text(size = 11),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) -> p

ggsave('f6_e.pdf', p, width = 3.6, height = 4.8, bg = 'transparent')






uclr <- '#76b7b2'
gene_name <- 'DLX1'
txii$abundance[names(e2g[e2g == gene_name]),] %>% 
  data.frame(v = .) %>% rownames_to_column('samp') %>%
  merge(rownames_to_column(md, 'samp')) %>%
  filter(line == 'BT245') %>%
  filter(media == 'SCM') %>% 
  filter(!grepl('D', samp)) %>%
  mutate(gene_name = gene_name) -> dat1
salmonn$abundance[names(e2g[e2g == gene_name]),] %>%
  data.frame(v = .) %>% rownames_to_column('samp') %>%
  merge(rownames_to_column(mdatt, 'samp')) %>%
  mutate(media = 'DM', gene_name = gene_name) %>%
  filter(line == 'BT245') %>%
  dplyr::select(samp, v, line, media, geno = grp, gene_name) %>%
  rbind(dat1) %>%
  arrange(v) %>%
  mutate(media = c('DM' = 'Differentiation media', 'SCM'  = 'Stem cell media')[media] %>%
           factor(c('Stem cell media', 'Differentiation media')),
         geno = sub(' ctrl', '', geno) %>%
           sub('K27M-KO', 'KO', .) %>%
           sub('K27M 4976', 'K27M + CBX-AM', .) %>%
           factor(c('K27M','KO','K27M + CBX-AM'))) %>%
  ggplot(aes(x = geno, y = v, color = geno, fill = geno)) + 
  geom_boxplot(position = position_nudge(x = .3), width = .2,
               outlier.color = NA) +
  geom_jitter(width = .1, height = 0) +
  #geom_signif(comparisons = list(c('K27M','KO')),
  #            test = 't.test', color = 'black', tip_length = 0,
  #            map_signif_level = signif.num) +
  scale_color_manual(values = c(tclrs, uclr)) +
  scale_fill_manual(values = c(tclrs, uclr)) +
  stat_summary(geom = "crossbar", width=0.15, fatten=0, color="white", 
               fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) },
               position = position_nudge(x = .3)) +
  ylab('DLX1 expression (TPM)') +
  scale_y_continuous(expand = expansion(c(0.03, .1))) +
  facet_nested(.~'De-repression of H3K4me3- cPRC1 target'+media, scales = 'free_x') +
  coord_cartesian(ylim = c(0, 5)) +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black'),
        panel.grid = element_blank(),
        panel.grid.major = element_line(color = 'grey80', linetype = 'dashed'),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(size = 11, color = 'white'),
        axis.text = element_text(color = 'black', size = 11),
        legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())  -> p

ggsave('f6_f.pdf', p, width = 3.5, height = 2.8, bg = 'transparent')












pd <- read_csv('../data/quant/marker.IF.csv') %>% 
  filter(line == 'BT245') %>%
  mutate(cond = factor(cond, c('KO', 'H3K27M+CBX-AM', 'H3K27M')))

ann <- compare_means(signal ~ cond, pd) %>%
  mutate(y.position = case_when(group1 == 'KO' & group2 == 'H3K27M' ~ .75,
                                group1 == 'KO' & group2 == 'H3K27M+CBX-AM' ~ .2,
                                group1 == 'H3K27M+CBX-AM' & group2 == 'H3K27M' ~ .85),
         cond = group1)

ggplot(pd, aes(x = cond, y = signal, fill = cond, color = cond)) + 
  geom_boxplot(outlier.color = NA, width = .3, position = position_nudge(x = .2)) +
  stat_summary(geom = "crossbar", width=0.15, fatten=0, color="white",
               fun.data = function(x) c(y=median(x), ymin=median(x), ymax=median(x)),
               position = position_nudge(x = .2)) +
  geom_half_point(side = 'l', width = .2) +
  stat_pvalue_manual(data = ann[ann$group2 == 'H3K27M',], 
                     label = 'p.signif', coord.flip = T, tip.length = 0,
                     vjust = .2) +
  stat_pvalue_manual(data = ann[ann$group2 != 'H3K27M',],
                     label = 'p.signif', coord.flip = T,
                     vjust = 1.2,
                     tip.length = 0) +
  coord_flip(clip = 'off') +
  scale_fill_manual(values = c(tclrs[2], uclr, tclrs[1])) +
  scale_color_manual(values = c(tclrs[2], uclr, tclrs[1])) +
  scale_y_continuous(breaks = c(0, .5)) +
  scale_x_discrete(expand = expansion(add = c(.3, .5))) +
  ylab('Normalized SOX10 signal') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(color = 'black'),
        legend.position = 'none',
        axis.text = element_text(size = 11, color = 'black'),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.x = element_line(color = 'grey80', linetype = 'dashed')) -> p
p
ggsave('f6_g.pdf', p, height = 3.7, width = 2, bg = 'transparent')