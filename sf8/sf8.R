library(data.table)
library(tidyverse)
library(tximport)
library(rtracklayer)
library(DESeq2)
library(ggpubr)
library(ggalluvial)
library(pals)
library(enrichR)
library(ggh4x)
library(babelgene)
library(shadowtext)
library(broom)
library(patchwork)
library(scattermore)
library(eulerr)
library(ggdist)
library(InteractionSet)

setEnrichrSite("Enrichr")
dbs <- listEnrichrDbs()
db.use <- c('WikiPathway_2021_Human', 'BioPlanet_2019')

cclrs <- kelly()[c(3,7,6,2)]

cls <- list()


cls$hiPSC <- readRDS('../data/pmtcgi/hiPSC.embed.rds') %>%
  mutate(cl = case_when(clust == '1_0' ~ '1',
                        clust == '0' ~ '4',
                        clust == '1_1_1' ~ '2',
                        clust == '1_1_0' ~ '3',
                        T ~ NA_character_),
         cl2 = c('1' = 'Active', '2' = 'cPRC1', '3' = 'PRC2', '4' = 'Other')[cl] %>%
           factor(c('Active', 'cPRC1', 'PRC2', 'Other'))) %>%
  merge(mutate(fread('../data/pmtcgi/hg38.5kb.bed',
                     col.names = c('chr', 'start', 'end', 'name', 'score', 'strand')),
               idx = 1:n()), by = 'idx')

cls$mESC <- readRDS('../data/pmtcgi/mESC.embed.rds')  %>%
  mutate(cl = case_when(clust == '1_0' ~ '1',
                        clust %in% '0' ~ '4',
                        clust == '1_1_1_1' ~ '3',
                        clust == '1_1_0' ~ '2x',
                        clust == '1_1_1_0' ~ '2y',
                        T ~ NA_character_),
         cl2 = c('1' = 'Active', '2' = 'cPRC1', '3' = 'PRC2',
                 '4' = 'Other')[sub('[xy]$', '', cl)] %>%
           factor(c('Active', 'cPRC1', 'PRC2', 'Other'))) %>%
  merge(mutate(fread('../data/pmtcgi/mm10.5kb.bed',
                     col.names = c('chr', 'start', 'end', 'name', 'score', 'strand')),
               idx = 1:n()), by = 'idx')


txs <- list('hiPSC' = '../data/genes/gencode.v36.annotation.gff3.gz',
            'mESC' = '../data/genes/gencode.vM25.annotation.gff3.gz') %>%
  lapply(import)

fss <- list()

fss$hiPSC <- list.files('../data/rna/ncrm1', full.names = T) %>%
  grep('iPSC-(cE4|cN38)', ., value = T) %>%
  setNames(., sub('.sf', '', basename(.)))

fss$mESC <- list.files('../data/rna/mesc', full.names = T) %>%
  setNames(., sub('.sf', '', basename(.)))


Map(function(tx, fs, cl) {
  md <- tibble(s = names(fs)) %>%
    mutate(cond = ifelse(grepl('cE4|A1_ESC', s), 'WT', 'KO')) %>%
    filter(cond == 'WT') %>%
    column_to_rownames('s')
  e2g <- as_tibble(tx[tx$type == 'gene']) %>% 
    dplyr::select(name = gene_id, gene_name, gene_type)
  tx[tx$type == 'transcript'] %>%
    as_tibble() %>%
    dplyr::select(transcript_id, gene_id) %>%
    distinct() %>%
    tximport(fs[rownames(md)], 'salmon', tx2gene = .) %>%
    .$abundance %>%
    rowMeans(na.rm = T) %>% 
    data.frame(v = .) %>% 
    rownames_to_column('name') %>% 
    merge(cl[,c('cl2','name', 'chr', 'start', 'end')], by = 'name') %>%
    merge(e2g, by = 'name')
}, txs, fss, cls) %>%
  bind_rows(.id = 'gn') %>%
  na.omit() -> ee

ann <- split(ee, ee$gn) %>% lapply(function(x) compare_means(v ~ cl2, x)) %>%
  bind_rows(.id = 'gn') %>% filter(group1 == 'cPRC1' | group2 == 'cPRC1') %>%
  mutate(y = rep(c(5,3,4), 2))

ee %>%
  ggplot(aes(x = cl2, y = v)) +
  geom_violin(aes(fill = cl2), color = NA, alpha = .5) +
  geom_boxplot(aes(color = cl2, fill = cl2), notch = T, width = .3) + 
  stat_summary(geom = "crossbar", width=0.1, fatten=0, color="white", 
               fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
  stat_pvalue_manual(data = ann, y = 'y', label = 'p.signif',
                     tip.length = 0) +
  facet_wrap(~gn, nrow = 1) +
  scale_y_log10('Gene expression (TPM)', breaks = c(0.01, 1, 100, 10000),
                labels = c('0.01', '1', '100', '10000')) +
  scale_color_manual(values = cclrs) +
  scale_fill_manual(values = cclrs) +
  coord_cartesian(ylim = c(1e-3, 1e5)) +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black'),
        strip.background = element_rect(fill = 'black'),
        panel.grid = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        panel.grid.major = element_line(color ='grey75', linetype = 'dashed'),
        strip.text = element_text(color = 'white', size = 11),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))->p

ggsave('sf8_a.pdf', p, height = 4, width = 4.9, bg = 'transparent')




Map(function(cl, tx) {
  cl %>% na.omit() %>% filter(cl2 == 'cPRC1') %>%
    pull(name) %>% {tx[tx$type == 'gene' & tx$gene_id %in% .]} %>% 
    .$gene_name %>%
    enrichr(db.use) 
}, cls, txs) ->  rres 


lapply(rres, `[[`, 'WikiPathway_2021_Human') %>% 
  bind_rows(.id = 'gn') %>% 
  group_by(gn) %>%
  slice_min(Adjusted.P.value, n = 5, with_ties = F) %>%
  ungroup() %>%
  arrange(gn, -P.value) %>%
  mutate(pway = sub(' WP.*', '', Term) %>% paste0(gn, ':', .) %>% fct_inorder(),
         olap = sub('/.*', '', Overlap) %>% as.integer()) %>%
  ggplot(aes(y = pway, x = Odds.Ratio, color = -log10(Adjusted.P.value))) +
  geom_segment(aes(yend = pway, xend = 0)) +
  geom_point(aes(size = olap)) +
  scale_x_continuous(expand = expansion(c(0, .05))) +
  scale_size('Overlap') +
  scale_color_gradientn(expression(-log[10]~'FDR'), colors = pals::viridis(50)) +
  labs(x =  'Odds ratio', y = 'Pathway') +
  facet_grid2(gn ~ 'cPRC1 targets: WikiPathways', scales = 'free_y') +
  scale_y_discrete(labels = function(x) sub('.*:', '', x)) +
  theme(#legend.direction = "vertical",
    #legend.box = "horizontal",
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 11, color = 'black'),
    axis.text = element_text(size = 11, color = 'black'),
    strip.background = element_rect(fill = 'black'),
    strip.text = element_text(size = 11, color = 'white'),
    panel.background = element_rect(fill = NA, color = 'black', size = 1),
    plot.background = element_blank(),
    panel.grid.major.x = element_line(color = 'grey80', linetype = 'dashed'),
    panel.grid = element_blank()) -> plt3



h2m <- orthologs(sub('\\..*', '', cls$hiPSC$name), species = 'mouse')


ortho <- cls$hiPSC %>%
  dplyr::select(human_ensembl = name, human_clust = cl2) %>%
  mutate(human_ensembl = sub('\\..*', '', human_ensembl)) %>%
  merge(., orthologs(.$human_ensembl, species = 'mouse')) %>%
  merge(mutate(dplyr::select(cls$mESC, ensembl = name, mouse_clust = cl2), 
               ensembl = sub('\\..*', '', ensembl)), by = 'ensembl') 

ortho %>%
  dplyr::count(mouse_clust, human_clust) %>%
  na.omit() %>%
  to_lodes_form(axes = 1:2) %>%
  mutate(allu = case_when(alluvium == 1 ~ 'Active',
                          alluvium == 6 ~ 'cPRC1',
                          alluvium == 11 ~ 'PRC2',
                          alluvium == 16 ~ 'Other',
                          T ~ NA_character_) %>% 
           factor(c('Active','cPRC1','PRC2','Other'))) %>%
  ggplot(aes(x = x, stratum = stratum, alluvium = alluvium, y = n)) +
  geom_alluvium(aes(fill = allu)) +
  geom_stratum(aes(fill = stratum)) +
  geom_shadowtext(stat = 'stratum', aes(label = stratum), color = 'white') +
  scale_x_discrete(limits = c('mouse_clust', 'human_clust'),
                   labels = c('mESC', 'hiPSC'),
                   expand = expansion(.2)) +
  scale_y_continuous(expand = expansion(0),
                     limits = c(0, 2e4)) +
  scale_fill_manual(values = cclrs, na.translate = T) +
  labs(x = 'Promoter classification', y = 'Number of orthologous genes') +
  facet_wrap(~'Classification conservation') +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black'),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = 'black'),
        legend.position = 'none',
        panel.grid.major = element_line(color = 'grey80', linetype = 'dashed'),
        strip.text = element_text(size = 11, color = 'white'),
        axis.text = element_text(size = 11, color = 'black')) -> plt1

plt2 <- levels(ortho$human_clust) %>%
  lapply(function(x) {
    fisher.test(ortho$human_clust == x,
                ortho$mouse_clust == x) %>%
      tidy() %>%
      mutate(cl = x)
  }) %>%
  bind_rows() %>%
  mutate(cl = fct_inorder(cl)) %>%
  ggplot(aes(x = cl, y = estimate, ymin = conf.low, ymax = conf.high, color = cl)) +
  geom_linerange(size = 1) +
  geom_point(size = 3) +
  ylab('Overlap odds ratio for\northologous gene classification') +
  facet_wrap(~'mESC vs hiPSC promoter class') +
  scale_color_manual(values = cclrs) +
  labs(x = 'Promoter chromatin classification') +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black'),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = 'black'),
        legend.position = 'none',
        panel.grid.major = element_line(color = 'grey80', linetype = 'dashed'),
        strip.text = element_text(size = 11, color = 'white'),
        axis.text = element_text(size = 11, color = 'black'))

{wrap_plots(plt3, plt1, plt2,nrow = 1, widths = c(.95,.85,1)) &
    theme(plot.background = element_blank())} %>%
  ggsave('sf8_bcd.pdf', ., height = 4.2, width = 12.2, bg = 'transparent')

                                       
load('../data/pmtcgi/hiPSC.cts.rda')
rbind(tibble(chip = c('hNPC_WT_CBX2',
                      'hNPC_WT_H3K4me3',
                      'hNPC_WT_H3K27me3',
                      'hNPC_WT_H3K36me2',
                      'hNPC_WT_H3K36me3'),
             inp = 'hNPC_WT_Input',
             type = 'NPC'),
      tibble(chip = c('hiPSC_WT_CBX2',
                      'hiPSC_WT_H3K4me3',
                      'hiPSC_WT_H3K27me3',
                      'hiPSC_WT_H3K36me2',
                      'hiPSC_WT_H3K36me3'),
             inp = 'hiPSC_WT_Input',
             type = 'iPSC')) %>%
  mutate(mark = sub('-R2$', '', chip) %>% sub('.*[-_]', '', .)) %>%
  split(., .$type) %>%
  lapply(function(x) {
    split(x, x$mark) %>%
      lapply(function(y) {
        facs <- c(libsz[y$chip], libsz[y$inp]) %>%
          {min(.) / .}
        log2(((cts[[y$chip]] * facs[1]) + 1) / ((cts[[y$inp]] * facs[2]) + 1))
      }) %>%
      bind_cols() %>%
      mutate(idx = 1:n(), .before = 1)
  }) %>%
  bind_rows(.id = 'type') -> oo


clrs <- c('#00000000', paste0(pals::brewer.greens(60)[-c(1:11)], 'FF'))
ooo <- oo[apply(oo[,-c(1:2)], 1, sd) > 0,] %>%
  add_count(idx) %>%
  filter(n == 2) %>%
  dplyr::select(-n) 

tclrs <- c('#e45756', '#4c78a8')
yl <- c(-5, 7.5)
yb <- c(-4, 0, 4)
xl <- c(-3, 6.5)
xb <- c(-3, 0, 3, 6)
ggplot(ooo, aes(x = CBX2, y = H3K27me3)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_scattermore(aes(color = type), pointsize = .5, alpha = .5) +
  annotate('rect', xmin = 2, xmax = 6.5, ymin = 4, ymax = 7.4, fill =NA,
           color= 'green', size = 2) +
  scale_y_continuous(breaks = yb) +
  scale_x_continuous(breaks = xb) +
  scale_color_manual(values = tclrs)+
  coord_cartesian(expand = F, xlim = xl, ylim = yl) +
  #facet_wrap(~'PRC enrichment\nin iPSCs & NPCs') +
  theme(plot.background = element_blank(),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(size = 11, color = 'white'),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_line(color = 'grey80', linetype = 'dashed'),
        legend.position = c(1, 0),
        axis.text = element_text(size = 11, color = 'black'),
        legend.text = element_text(size = 11, color = 'black'),
        legend.justification = c(1, 0),
        legend.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1)) -> p1

ggplot(ooo, aes(x = H3K27me3, fill = type)) + 
  geom_density(alpha = .5) +
  coord_flip(xlim = yl) +
  scale_x_continuous(breaks = yb, expand = expansion(0)) +
  scale_y_continuous('Density', breaks = c(0, .3), expand = expansion(c(0, .05))) +
  scale_fill_manual(values = tclrs) +
  theme(plot.background = element_blank(),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(size = 11, color = 'white'),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_line(color = 'grey80', linetype = 'dashed'),
        #legend.position = c(1, 0),
        axis.text = element_text(size = 11, color = 'black'),
        legend.text = element_text(size = 11, color = 'black'),
        legend.justification = c(1, 0),
        legend.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        axis.title.y = element_blank(),
        legend.position = 'none',
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) -> p2

ggplot(ooo, aes(x = CBX2, fill = type)) + 
  geom_density(alpha = .5) +
  coord_cartesian(xlim = xl) +
  scale_x_continuous(breaks = xb, expand = expansion(0)) +
  scale_y_continuous('Density', breaks = c(0, .5), expand = expansion(c(0, .05))) +
  facet_wrap(~'PRC targets\nin iPSCs & NPCs') +
  scale_fill_manual(values = tclrs) +
  theme(plot.background = element_blank(),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(size = 11, color = 'white'),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_line(color = 'grey80', linetype = 'dashed'),
        legend.position = c(1, 0),
        axis.text = element_text(size = 11, color = 'black'),
        legend.text = element_text(size = 11, color = 'black'),
        legend.justification = c(1, 0),
        legend.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) -> p3

{wrap_plots(p3, plot_spacer(),  p1, p2,heights = c(1,3), widths = c(3, 1)) &
    theme(plot.background = element_blank())} %>%
  ggsave('sf8_e.joint.pdf', ., height = 3.8, width = 3.5, bg = 'transparent')


filter(oo, CBX2 > 2 & H3K27me3 > 4) %>%
  {split(.$idx, .$type)} -> idxs

u <- unique(unlist(idxs))

pdf('sf8_e.euler.pdf', height = 1, width = 1.9)
lapply(idxs, function(x) u %in% x) %>%
  bind_cols() %>%
  euler() %>%
  plot(fills = tclrs, quantities = T)
dev.off()



readRDS('../data/pmtcgi/GBM.embed.rds') %>%
  mutate(cl = case_when(clust == '0' ~ '1',
                        clust == '1_1' ~ '4',
                        clust %in% c('1_0_0', '1_0_3_0') ~ '2b',
                        clust %in% c('1_0_1', '1_0_3_1') ~ '2a',
                        clust %in% c('1_0_2_0_0', '1_0_2_1_1') ~ '3b',
                        clust %in% c('1_0_2_0_1', '1_0_2_1_0') ~ '3a',
                        T ~ NA_character_),
         cl2 = sub('[ab]$', '', cl)) %>%
  na.omit() %>%
  filter(cl2 == '2') %>%
  dplyr::select(line, idx) %>%
  dplyr::count(idx) %>%
  filter(n == 3) %>% pull(idx) -> c2

pd <- fread('../data/pmtcgi/hg38.5kb.bed') %>%
  mutate(idx = 1:n(),
         kind = ifelse(grepl('ENSG', V4), 'Promoter', 'CGI')) %>%
  {c(split(.$idx, .$kind), list(All = .$idx))} %>%
  lapply(function(x) {
    lapply(idxs, function(y) {
      a <- intersect(x, c2)
      b <- intersect(x, y)
      
      d <- c(sum(a %in% b), sum(!(a %in% b)), sum(!(b %in% a))) %>%
        {c(., length(x) - sum(.))} 
      
      matrix(d, 2, 2) %>%
        fisher.test() %>%
        tidy() %>%
        cbind(data.frame(t(d)))
    }) %>% bind_rows(.id = 'type')
  }) %>% bind_rows(.id = 'kind') %>%
  filter(kind == 'Promoter') 

ggplot(pd, aes(x = type, y = estimate, ymin = conf.low, ymax = conf.high,
               fill = type, color = type)) +
  geom_col(aes(y = X1), alpha = .3, color = NA) +
  geom_linerange(size = 1) +
  geom_point(size = 5) +
  scale_y_continuous('Odds ratio of overlap',
                     expand = expansion(c(0, .05)),
                     sec.axis = dup_axis(name = 'Size of overlap')) +
  facet_wrap(~'Intersect with glioma\ncPRC1 targets') +
  scale_color_manual(values = tclrs) +
  scale_fill_manual(values = tclrs) +
  labs(x = 'Cell type for defining PRC targets') +
  theme(plot.background = element_blank(),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(size = 11, color = 'white'),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_line(color = 'grey80', linetype = 'dashed'),
        legend.position = 'none',
        axis.text = element_text(size = 11, color = 'black'),
        axis.title.y.right = element_text(color = 'grey50'),
        axis.text.y.right = element_text(color = 'grey50'),
        panel.background = element_rect(fill = NA, color = 'black', size = 1)) -> pp1

signif.num <- function(x) {
  as.character(symnum(x, corr = FALSE, na = FALSE, legend = FALSE,
                      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                      symbols = c("****", "***", "**", "*", "ns")))
}


list('PRC\ntargets' = filter(ooo, CBX2 > 2 & H3K27me3 > 4),
     'All\npromoters' = ooo) %>%
  bind_rows(.id = 'kind') %>%
  filter(idx %in% grep('ENSG', fread('../data/pmtcgi/hg38.5kb.bed')$V4)) %>%
  ggplot(aes(x = type, y = H3K4me3, fill = type, color = type)) +
  geom_hline(yintercept = 0) +
  geom_violin(color = NA, alpha = .4) +
  stat_pointinterval() +
  facet_wrap(~kind, nrow = 1) +
  ylab('H3K4me3') +
  scale_y_continuous(expand = expansion(c(0, .1))) +
  geom_signif(comparisons = list(c('iPSC', 'NPC')), color = 'black',
              map_signif_level = signif.num) +
  scale_fill_manual(values = tclrs) +
  scale_color_manual(values = tclrs) +
  theme(plot.background = element_blank(),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(size = 11, color = 'white'),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        legend.key = element_blank(),
        panel.grid.major = element_line(color = 'grey80', linetype = 'dashed'),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text = element_text(size = 11, color = 'black'),
        panel.background = element_rect(fill = NA, color = 'black', size = 1)) -> pp2
{wrap_plots(pp1, pp2, nrow = 1) &
    theme(plot.background = element_blank())} %>%
  ggsave('sf8_fg.pdf', ., height = 3.56, width = 5.8, bg = 'transparent')







load('../data/pmtcgi/hg38.chromosight.rda')
mat <- mat[, c('hiPSC_WT','hNPC_WT')]

i <- bins %>%
  {list(.[,1:3], .[,4:6])} %>%
  lapply(function(x) {
    x %>%
      `colnames<-`(c('chr','start','end')) %>%
      makeGRangesFromDataFrame()
  }) %>%
  {GInteractions(.[[1]], .[[2]], mode = 'strict')} 

cls <- readRDS('../data/pmtcgi/hiPSC.embed.rds') %>%
  mutate(cl = case_when(clust == '1_0' ~ '1',
                        clust == '0' ~ '4',
                        clust == '1_1_1' ~ '2',
                        clust == '1_1_0' ~ '3',
                        T ~ NA_character_),
         cl2 = c('1' = 'Active', '2' = 'cPRC1', '3' = 'PRC2', '4' = 'Other')[cl] %>%
           factor(c('Active', 'cPRC1', 'PRC2', 'Other'))) %>%
  merge(mutate(fread('../data/pmtcgi/hg38.5kb.bed',
                     col.names = c('chr', 'start', 'end', 'name', 'score', 'strand')),
               idx = 1:n()), by = 'idx') %>%
  na.omit() %>%
  mutate(start = start + 2499,
         end = end - 2499) %>%
  split(., .$cl2) %>%
  lapply(makeGRangesFromDataFrame, keep.extra.columns = T)

lapply(cls, function(x) {
  p <- overlapsAny(anchors(i)$first, x) & 
    overlapsAny(anchors(i)$second, x)
  
  rr <- fread('../data/pmtcgi/hg38.5kb.bed', select = 1:3,
              col.names = c('chr','start','end'))
  
  mat[p,] %>% 
    cbind(bins[p,]) %>%
    na.omit() %>%
    mutate(dif = hiPSC_WT  - hNPC_WT,
           c1 = sprintf('%s:%d-%d', chrom1, start1, end1),
           c2 = sprintf('%s:%d-%d', chrom2, start2, end2)) %>%
    dplyr::select(-chrom1,-start1,-end1,-chrom2,-start2,-end2) %>%
    pivot_longer(c(c1,c2), names_to = 'c', values_to = 'i') %>%
    separate(i,c('seqnames','start','end'),'[:-]') %>%
    dplyr::select(-c) %>%
    distinct() %>%
    merge(as_tibble(x)[,c('seqnames','start','end','name', 'idx')]) -> dd
  
  list(All = dd,
       Promoter = dd[grepl('ENS', dd$name),],
       CGI = dd[!grepl('ENS', dd$name),]) %>%
    lapply(function(y) {
      y %>%
        group_by(name, idx) %>%
        summarise(mu = mean(dif)) %>%
        ungroup()
    }) %>%
    bind_rows(.id = 'kind')
}) %>%
  bind_rows(.id = 'clu') -> ii

oo %>%
  pivot_wider(names_from = type, values_from = c('CBX2','H3K27me3', 'H3K36me2', 'H3K36me3','H3K4me3')) %>% 
  arrange(idx) %>%
  merge(ii, by = 'idx') %>%
  filter(clu == 'cPRC1') %>%
  mutate(CBX2_Diff = CBX2_iPSC - CBX2_NPC,
         H3K4me3_Diff = H3K4me3_iPSC - H3K4me3_NPC) %>%
  dplyr::select(idx, kind, mu, contains('CBX2'), contains('H3K4')) %>%
  mutate(qmu = ntile(mu, 10) %>% as.character() %>% fct_inseq()) %>%
  pivot_longer(-c(idx, kind, mu, qmu), names_sep = '_', names_to = c('mark', 'type')) %>%
  mutate(grp = ifelse(mu > 0, 'iPSC > NPC', 'iPSC < NPC')) %>%
  filter(kind == 'All') -> pdat

pdat %>%
  filter(type == 'Diff') %>%
  mutate(grp = factor(grp, c('iPSC > NPC', 'iPSC < NPC'))) %>%
  ggplot(aes(x = grp, y = value, group = interaction(grp, type), fill = grp, color = grp)) +
  geom_violin(color = NA, alpha = .3) +
  geom_boxplot(notch = T, width = .2) +
  stat_summary(geom = "crossbar", width=0.1, fatten=0, color="white", position = position_dodge(.75),
               fun.data = function(x) c(y=median(x), ymin=median(x), ymax=median(x))) +
  facet_nested(.~'Determinants of differential cPRC1 loops'+mark, scales = 'free_x') +
  geom_signif(comparisons = list(c('iPSC > NPC', 'iPSC < NPC')), 
              map_signif_level = signif.num, color = 'black') +
  scale_fill_manual(values = tclrs) +
  scale_color_manual(values = tclrs) +
  scale_y_continuous(expand = expansion(c(0, .1))) +
  labs(y = 'Differential ChIP-seq enrichment (iPSC - NPC)', x = 'Average loop score with other cPRC1 targets') +
  theme(legend.position = 'none',
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(size = 11, color = 'white'),
        axis.text = element_text(size = 11, color=  'black'),
        panel.grid = element_blank(),
        legend.key = element_blank(),
        axis.text.x = element_text(angle= 45, hjust = 1),
        legend.background = element_blank(),
        legend.text = element_text(size = 11),
        legend.title = element_blank(),
        plot.background = element_blank(),
        panel.grid.major = element_line(color = 'grey80', linetype = 'dashed'),
        panel.background = element_rect(fill = NA, color = 'black', size = 1)) -> pplt
ggsave('sf8_h.pdf', pplt, height = 4, width = 4, bg = 'transparent')


