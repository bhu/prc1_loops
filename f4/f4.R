library(data.table)
library(tidyverse)
library(rtracklayer)
library(patchwork)
library(tximport)
library(pals)
library(ggnewscale)
library(scattermore)
library(ggh4x)
library(DESeq2)
library(InteractionSet)
library(broom)

cclrs <- kelly()[c(3,7,6,2)]
pd2 <- readRDS('../data/pmtcgi/mESC.embed.rds') %>%
  mutate(cl = case_when(clust == '1_0' ~ '1',
                        clust %in% '0' ~ '4',
                        clust == '1_1_1_1' ~ '3',
                        clust %in% c('1_1_0', '1_1_1_0') ~ '2',
                        T ~ NA_character_),
         cl2 = c('1' = 'Active', '2' = 'cPRC1', '3' = 'PRC2', '4' = 'Other')[cl] %>%
           factor(c('Active', 'cPRC1', 'PRC2', 'Other'))) %>%
  na.omit() 

p3 <- ggplot(pd2, aes(x = V1, y = V2)) +
  geom_scattermore(aes(color = cl2),  pixels = c(350,350)) +
  labs(x = 'UMAP 1', y = 'UMAP 2') +
  facet_wrap(~'mESC CGIs & promoters') +
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
        legend.position = c(0,0),
        legend.justification = c(0,0),
        legend.text = element_text(size = 11),
        strip.text = element_text(color = 'white', size = 11)) 

p4 <- pd2 %>%
  select(V1, V2, H3K27me3, idx, CBX2, H3K4me3, H3K27ac) %>%
  pivot_longer(c(H3K27me3, H3K27ac, H3K4me3, CBX2), 
               names_to = 'mark', values_to = 'v') %>%
  ggplot(aes(x = V1, y = V2)) +
  geom_scattermore(aes(color = v), pixels = c(350,350)) +
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



readRDS('../data/pmtcgi/hiPSC.embed.rds') %>%
  mutate(cl = case_when(clust == '1_0' ~ '1',
                        clust == '0' ~ '4',
                        clust == '1_1_1' ~ '2',
                        clust == '1_1_0' ~ '3',
                        T ~ NA_character_),
         cl2 = c('1' = 'Active', '2' = 'cPRC1', '3' = 'PRC2', '4' = 'Other')[cl] %>%
           factor(c('Active', 'cPRC1', 'PRC2', 'Other'))) %>%
  na.omit() -> pd3

p5 <- ggplot(pd3, aes(x = V1, y = V2)) +
  geom_scattermore(aes(color = cl2),  pixels = c(350,350)) +
  labs(x = 'UMAP 1', y = 'UMAP 2') +
  facet_wrap(~'hiPSC CGIs & promoters') +
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
        legend.position = c(0,0),
        legend.justification = c(0,0),
        legend.text = element_text(size = 11),
        strip.text = element_text(color = 'white', size = 11)) 

p6 <- pd3 %>%
  select(V1, V2, H3K27me3, idx, CBX2, H3K4me3, H3K27ac) %>%
  pivot_longer(c(H3K27me3, H3K27ac, H3K4me3, CBX2), 
               names_to = 'mark', values_to = 'v') %>%
  ggplot(aes(x = V1, y = V2)) +
  geom_scattermore(aes(color = v), pixels = c(350,350)) +
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

{ wrap_plots(p3,p4,p5,p6,nrow =2) &
    theme(plot.background = element_blank()) } %>%
ggsave('f4_a.pdf', height = 5.2, width = 5.9, bg = 'transparent')


tx2 <- import.gff('../data/genes/gencode.vM25.annotation.gff3.gz')
tx2gene2 <- as_tibble(tx2) %>% filter(type == 'transcript') %>% 
  distinct(transcript_id, gene_id) %>% dplyr::select(2:1)
txi2 <- list.files('../data/rna/mesc', pattern = 'sf', full.names = T) %>%
  setNames(., sub('.sf', '', basename(.))) %>%
  tximport(type = "salmon", tx2gene = tx2gene2) 


m2 <- c('Expression', "H3K4me1", 'H3K4me3', 'H3K27ac','H3K9ac','MLL2',
        'H3' , 'RPB1', 'Pol2Ser5P', 'CTCF', 
        'H3K36me3', 'H3K36me2', 'H3K9me3',
        'RING1B', 'H2Aub', 'CBX2','CBX7','PHC1', 'RYBP',
        'SUZ12','EZH2','JARID2','AEBP2','PCL2', 'EPOP',
        'H3K27me3','H3K27me2', 'H3K27me1')

list(ENCODE = c('H3K4me1', 'H3K9ac', 'H3K9me3'),
     Mas2018 = c('CTCF', 'H3', 'MLL2', 'Pol2Ser5P', 'RPB1', 'H3K27ac', 'H3K4me3'),
     Healy2019 = c('EPOP', 'RYBP', 'SUZ12', 'PCL2', 'JARID2', 'H2Aub', 'CBX7', 'AEBP2', 'H3K27me1', 'H3K27me2'),
     Chen2022 = c('H3K36me2', 'H3K36me3', 'H3K27me3'),
     `This study` = c('CBX2', 'RING1B'),
     Kundu2017 = c('EZH2', 'PHC1')) %>%
  lapply(function(x) tibble(mark = x)) %>%
  bind_rows(.id = 'proj') %>%
  mutate(symb = c('*', '°', '^', '"', '`', '')[factor(proj)],
         nm = paste0(mark, symb)) %>%
  dplyr::select(mark, nm) %>% deframe() -> rnmmrks


pdd2 <- txi2$abundance[,c("A1_ESC_rep1", "A1_ESC_rep2")] %>%
  rowMeans(na.rm = T)  %>%
  data.frame(Expression = .) %>%
  rownames_to_column('V4') %>%
  merge(mutate(fread('../data/pmtcgi/mm10.5kb.bed'), idx = 1:n())) %>%
  mutate(Expression = log2(Expression + 1)) %>%
  select(idx, Expression) %>%
  as.data.table() %>%
  merge(pd2, by = 'idx') %>%
  {.[,setdiff(names(.), c('cl', 'clust', 'idx', 'line', 'V1','V2')),with=F]} %>%
  pivot_longer(-cl2, names_to ='mark', values_to = 'v') %>%
  group_by(cl2, mark) %>%
  summarise(v = mean(v, na.rm = T), .groups = 'drop') %>%
  mutate(mark = factor(mark, m2),
         x = as.numeric(cl2),
         y = as.numeric(mark)) %>%
  arrange(mark) %>%
  mutate(mark = ifelse(mark == 'Expression', 'Expression', rnmmrks[as.character(mark)]) %>% fct_inorder())

pp2 <- ggplot(pdd2, aes(x = x, y = y, fill = v)) +
  geom_tile(aes(fill = v), data = ~subset(., mark=='Expression')) +
  scale_fill_gradientn('log2(TPM + 1)', colors = cividis(25), breaks = c(1, 3),
                       limits = range(pdd2$v[pdd2$mark == 'Expression']),
                       guide = guide_colorbar(barheight = .5, barwidth = 2,
                                              title.position = 'left',
                                              title.vjust = 1, order = 1)) +
  new_scale_fill() +
  geom_tile(aes(fill = v), data = ~subset(., mark !='Expression')) +
  scale_fill_gradientn('Enrichment', colors = coolwarm(25), 
                       limits = c(-3, 3), breaks = c(-3, 0, 3), oob = scales::squish,
                       guide = guide_colorbar(barheight = .5, barwidth = 2,
                                              title.position = 'left',
                                              title.vjust = 1, order = 2)) +
  geom_hline(yintercept = 1.5, size = 1, color = 'white') +
  scale_x_reverse(breaks = unique(pdd2$x), labels = unique(pdd2$cl2)) +
  scale_y_continuous(breaks = unique(pdd2$y), labels = unique(pdd2$mark)) +
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


tx <- import.gff('../data/genes/gencode.v36.annotation.gff3.gz')
tx2gene <- as_tibble(tx) %>% filter(type == 'transcript') %>% 
  distinct(transcript_id, gene_id) %>% dplyr::select(2:1)
txi <- list.files('../data/rna/ncrm1', pattern = 'sf',full.names = T) %>%
  setNames(., sub('.sf', '', basename(.))) %>%
  tximport(type = "salmon", tx2gene = tx2gene)

m3 <- c('Expression', 'H3K4me3', 'H3K27ac','CTCF', 
        'H3K36me3', 'H3K36me2', 'H3K9me3',
        'RING1B', 'H2Aub', 'CBX2',
        'SUZ12', 'H3K27me3')



pdd3 <- txi$abundance[,c("NCRM1-iPSC-cE4-p23" , "NCRM1-iPSC-cE4-p25")] %>%
  rowMeans(na.rm = T)  %>%
  data.frame(Expression = .) %>%
  rownames_to_column('V4') %>%
  merge(mutate(fread('../data/pmtcgi/hg38.5kb.bed'), idx = 1:n())) %>%
  mutate(Expression = log2(Expression + 1)) %>%
  select(idx, Expression) %>%
  as.data.table() %>%
  merge(pd3, by = 'idx') %>%
  {.[,setdiff(names(.), c('cl', 'clust', 'idx', 'line', 'V1','V2', 'type')),with=F]} %>%
  pivot_longer(-cl2, names_to ='mark', values_to = 'v') %>%
  group_by(cl2, mark) %>%
  summarise(v = mean(v, na.rm = T), .groups = 'drop') %>%
  mutate(mark = factor(sub('H2AK119ub', 'H2Aub', mark), m3),
         x = as.numeric(cl2),
         y = as.numeric(mark)) 

pp3 <- ggplot(pdd3, aes(x = x, y = y, fill = v)) +
  geom_tile(aes(fill = v), data = ~subset(., mark=='Expression')) +
  scale_fill_gradientn('log2(TPM + 1)', colors = cividis(25), breaks = c(1, 3),
                       limits = range(pdd3$v[pdd3$mark == 'Expression']),
                       guide = guide_colorbar(barheight = .5, barwidth = 2,
                                              title.position = 'left',
                                              title.vjust = 1, order = 1)) +
  new_scale_fill() +
  geom_tile(aes(fill = v), data = ~subset(., mark !='Expression')) +
  scale_fill_gradientn('Enrichment', colors = coolwarm(25), 
                       limits = c(-3, 3), breaks = c(-3, 0, 3), oob = scales::squish,
                       guide = guide_colorbar(barheight = .5, barwidth = 2,
                                              title.position = 'left',
                                              title.vjust = 1, order = 2)) +
  geom_hline(yintercept = 1.5, size = 1, color = 'white') +
  scale_x_reverse(breaks = unique(pdd3$x), labels = unique(pdd3$cl2)) +
  scale_y_continuous(breaks = unique(pdd3$y), labels = unique(pdd3$mark)) +
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

{wrap_plots(pp2,pp3,ncol=1) &
    theme(plot.background = element_blank())} %>%
  ggsave('f4_b.pdf', ., height = 5, width = 5)


source('../data/diverging_map.R')
tclrs <- c('#e45756', '#4c78a8')
cmap <- diverging.colormap(seq(0, 1, .01),
                           rgb1 = hex2RGB(tclrs[2]),
                           rgb2 = hex2RGB(tclrs[1]),
                           outColorspace = "sRGB") %>%
  {.[. > 1] <- 1; .} %>%
  {rgb(.[,1], .[,2], .[,3])}
readRDS('../../reprod3/data/pup/PSC.clust.rds') -> pups


pd2 <- pups$mESC %>%
  pivot_wider(names_from = 'cond', values_from = 'z') %>%
  mutate(z = WT/`NSD1-KO`)

ann2 <- pd2 %>%
  filter(between(x, 11, 11) & between(y, 11, 11)) %>%
  group_by( clu) %>%
  summarise(z = mean(z,na.rm = T) %>% round(2)) %>%
  ungroup()

ggplot(pd2, aes(x, y, fill = z)) +
  geom_raster() +
  geom_label(aes(x = Inf, y = -Inf, hjust = 1, vjust = 1, label = z),
             data = ann2, inherit.aes = F,
             alpha = .8, label.size = NA) +
  scale_y_reverse() +
  coord_cartesian(expand = F) +
  scale_fill_gradientn(expression(frac('WT O/E','NSD1-KO O/E')),
                       colors = cmap, 
                       limits = c(1/1.5,1.5),
                       breaks = c(0.8, 1, 1.4),
                       oob = scales::squish,
                       trans = 'log2',
                       guide = guide_colorbar(label.position = 'bottom', 
                                              title.vjust = 1,
                                              barheight = 2)) +
  facet_nested(.~'mESC' + clu) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = c(1,0),
        legend.justification = c(1,1),
        plot.margin = margin(5,5,50,0),
        axis.text = element_blank(),
        legend.direction = 'horizontal',
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.background = element_blank(),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white')) -> plt2




pd3 <- pups$hiPSC %>%
  pivot_wider(names_from = 'cond', values_from = 'z') %>%
  mutate(z = WT/`NSD1+/-`)
ann3 <- pd3 %>%
  filter(between(x, 11, 11) & between(y, 11, 11)) %>%
  group_by( clu) %>%
  summarise(z = mean(z,na.rm = T) %>% round(2)) %>%
  ungroup()
ggplot(pd3, aes(x, y, fill = z)) +
  geom_raster() +
  geom_label(aes(x = Inf, y = -Inf, hjust = 1, vjust = 1, label = z),
             data = ann3, inherit.aes = F,
             alpha = .8, label.size = NA) +
  scale_y_reverse() +
  coord_cartesian(expand = F) +
  scale_fill_gradientn(expression(frac('WT O/E','NSD1-KO O/E')),
                       colors = cmap, 
                       limits = c(1/1.7,1.7),
                       breaks = c(.6, 1, 1.7),
                       oob = scales::squish,
                       trans = 'log2',
                       guide = guide_colorbar(label.position = 'bottom', 
                                              title.vjust = 1,
                                              barheight = 2)) +
  facet_nested(.~'hiPSC' +clu) +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = c(1,0),
        legend.justification = c(1,1),
        plot.margin = margin(5,5,50,0),
        axis.text = element_blank(),
        legend.direction = 'horizontal',
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.background = element_blank(),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white')) -> plt3

{wrap_plots(plt2,plt3,nrow=1) &
    theme(plot.background = element_blank())} %>%
  ggsave('f4_c.pdf', ., height = 2.6, width =9.9, bg = 'transparent')


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

deg <- Map(function(tx, fs, cl) {
  md <- tibble(s = names(fs)) %>%
    mutate(cond = ifelse(grepl('cE4|A1_ESC', s), 'WT', 'KO')) %>%
    column_to_rownames('s')
  e2g <- as_tibble(tx[tx$type == 'gene']) %>% 
    dplyr::select(name = gene_id, gene_name, gene_type)
  tx[tx$type == 'transcript'] %>%
    as_tibble() %>%
    dplyr::select(transcript_id, gene_id) %>%
    distinct() %>%
    tximport(fs, 'salmon', tx2gene = .) %>%
    DESeqDataSetFromTximport(md, ~cond) %>%
    DESeq() %>%
    results(c('cond', 'KO', 'WT')) %>%
    as.data.frame() %>% 
    rownames_to_column('name') %>%
    merge(cl[,c('cl2','name', 'chr', 'start', 'end')], by = 'name') %>%
    mutate(start = start + 2499,
           end = end - 2499) %>%
    merge(e2g, by = 'name')
}, txs, fss, cls)

list(hiPSC = list(rda = '../data/pmtcgi/hg38.chromosight.rda',
                  WT = 'hiPSC_WT',
                  KO = 'hiPSC_NSD1het'),
     mESC = list(rda = '../data/pmtcgi/mm10.chromosight.rda',
                 WT = 'PA',
                 KO = 'NSD1KO')) %>%
  Map(function(i, de) {
    load(i$rda)
    mat <- as_tibble(mat)[, c(i$WT, i$KO)]
    names(mat) <- c('WT', 'KO')
    
    i <- bins %>%
      {list(.[,1:3], .[,4:6])} %>%
      lapply(function(x) {
        x %>%
          `colnames<-`(c('chr','start','end')) %>%
          makeGRangesFromDataFrame()
      }) %>%
      {GInteractions(.[[1]], .[[2]], mode = 'strict')}
    
    de %>%
      split(., .$cl2) %>% 
      lapply(makeGRangesFromDataFrame, keep.extra.columns = T) %>%
      lapply(function(x) {
        p <- overlapsAny(anchors(i)$first, x) & 
          overlapsAny(anchors(i)$second, x)
        
        mat[p,] %>% 
          cbind(bins[p,]) %>%
          na.omit() %>%
          mutate(dif = KO  - WT,
                 c1 = sprintf('%s:%d-%d', chrom1, start1, end1),
                 c2 = sprintf('%s:%d-%d', chrom2, start2, end2)) %>%
          dplyr::select(-chrom1,-start1,-end1,-chrom2,-start2,-end2) %>%
          pivot_longer(c(c1,c2), names_to = 'c', values_to = 'i') %>%
          separate(i,c('seqnames','start','end'),'[:-]') %>%
          dplyr::select(-c) %>%
          distinct() %>%
          merge(as_tibble(x)[,c('name', 'seqnames', 'start', 'end')]) %>%
          group_by(name) %>%
          summarise(v = mean(dif)) %>%
          merge(as_tibble(x), by = 'name') 
      }) %>%
      bind_rows(.id = 'clu')
  }, ., deg) %>%
  bind_rows(.id = 'gn') -> ores

list(hiPSC = list(rda = '../data/pmtcgi/hg38.chromosight.rda',
                  WT = 'hiPSC_WT',
                  KO = 'hiPSC_NSD1het'),
     mESC = list(rda = '../data/pmtcgi/mm10.chromosight.rda',
                 WT = 'PA',
                 KO = 'NSD1KO')) %>%
  Map(function(i, de) {
    load(i$rda)
    mat <- as_tibble(mat)[, c(i$WT, i$KO)]
    names(mat) <- c('WT', 'KO')
    
    i <- bins %>%
      {list(.[,1:3], .[,4:6])} %>%
      lapply(function(x) {
        x %>%
          `colnames<-`(c('chr','start','end')) %>%
          makeGRangesFromDataFrame()
      }) %>%
      {GInteractions(.[[1]], .[[2]], mode = 'strict')}
    
    de %>%
      split(., .$cl2) %>% 
      lapply(makeGRangesFromDataFrame, keep.extra.columns = T) %>%
      lapply(function(x) {
        p <- overlapsAny(anchors(i)$first, x) & 
          overlapsAny(anchors(i)$second, x)
        
        mat[p,] %>% 
          cbind(bins[p,]) %>%
          na.omit() %>%
          mutate(dif = KO  - WT,
                 c1 = sprintf('%s:%d-%d', chrom1, start1, end1),
                 c2 = sprintf('%s:%d-%d', chrom2, start2, end2)) %>%
          dplyr::select(-chrom1,-start1,-end1,-chrom2,-start2,-end2) %>%
          pivot_longer(c(c1,c2), names_to = 'c', values_to = 'i') %>%
          separate(i,c('seqnames','start','end'),'[:-]') %>%
          dplyr::select(-c) %>%
          distinct() %>%
          merge(as_tibble(x)[,c('name', 'seqnames', 'start', 'end')]) %>%
          group_by(name) %>%
          summarise(v = mean(dif)) %>%
          merge(as_tibble(x), by = 'name') %>%
          dplyr::select(v, log2FoldChange) %>%
          na.omit() -> ddat
        
        c('pearson', 'spearman') %>%
          setNames(., .) %>%
          lapply(function(z) {
            cor.test(ddat$v, ddat$log2FoldChange, method = z) %>%
              tidy() 
          }) %>%
          bind_rows(.id = 'method')
      }) %>%
      bind_rows(.id = 'clu')
  }, ., deg) %>%
  bind_rows(.id = 'gn') -> cores

signif.num <- function(x) {
  as.character(symnum(x, corr = FALSE, na = FALSE, legend = FALSE,
                      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                      symbols = c("****", "***", "**", "*", "ns")))
}

cores %>%
  filter(method == 'pearson') %>% 
  mutate(clu = factor(clu, rev(levels(cls[[1]]$cl2))),
         symb = signif.num(p.value)) %>%
  ggplot(aes(x = gn, y = clu, fill = estimate)) + geom_tile() +
  scale_fill_gradientn('Pearson\'s r',
                       colors = pals::coolwarm(), limits = c(-.25, .25)) +
  facet_wrap(~'DiffExp vs DiffInt') +
  labs(x = 'Cell type', y = 'Promoter class') +
  coord_cartesian(expand = F) +
  theme(axis.text = element_text(color = 'black', size = 11),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(size = 11, color = 'white'),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        panel.grid.major = element_line(color = 'grey80', linetype = 'dashed'),
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.ticks = element_blank()) -> plt1

ores %>% 
  mutate(clu = factor(clu, levels(cls[[1]]$cl2))) %>%
  ggplot(aes(y = v, x = log2FoldChange, color = clu, fill = clu)) +
  geom_scattermore(pointsize = .9, data = ~subset(., clu != 'cPRC1'), alpha = .5) +
  geom_point(size = 1, shape = 16, data = ~subset(., clu == 'cPRC1'), alpha = .5) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_smooth(method = 'rlm') + 
  facet_wrap(~gn) +
  geom_point(size = 3, shape = 21, fill = NA, stroke = 2, show.legend = F,
             data = ~subset(., gene_name %in% c('EVX2','FGF5','HAND2','HOXB5','HOXB6','HOXB7','PAX9'))) +
  scale_color_manual('Promoter type', values = cclrs) +
  scale_fill_manual('Promoter type', values = cclrs) +
  labs(y = 'Differential intra-class P-P interaction\nNSD1+/- OR NSD1-/- over WT',
       x = expression(log[2]~frac('NSD1+/- or NSD1-/-', 'WT PSC')~'expression')) +
  coord_cartesian(xlim = c(-6, 6), ylim = c(-.3, .3)) +
  scale_x_continuous(breaks = c(-4, 0, 4)) +
  theme(plot.background = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(size = 11, color = 'white'),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        panel.grid.major = element_line(color = 'grey80', linetype = 'dashed'),
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.ticks = element_blank()) -> plt2

{wrap_plots(plt1, plt2,nrow=1, widths = c(1, 3)) &
    theme(plot.background = element_blank())} %>%
  ggsave('f4_de.pdf', ., height = 3.76, width = 9.5,  bg = 'transparent')





ddat <- fread('../data/quant/NCRM1.ddPCR.tsv') %>% 
  pivot_longer(-Gene, names_to = 'cond') %>% 
  mutate(cond = ifelse(grepl('WT', cond), 'WT', 'NSD1+/-'),
         Gene = factor(Gene)) 

rcts <- distinct(ddat, Gene) %>%
  mutate(idx = as.numeric(Gene),
         xmin = idx - .5, 
         xmax = idx + .5,
         ymin = 1e-10, ymax = Inf,
         alt = idx %% 2 == 0)


ggplot(ddat, aes(x = as.numeric(Gene), y = value, group = cond, color = cond)) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = alt),
            show.legend = F, data = rcts, inherit.aes = F) +
  scale_fill_manual(values = c('#ffffff00', '#00000022')) +
  stat_summary(geom = 'pointrange', fun.data = mean_se, position = position_dodge(.3)) +
  scale_y_log10(breaks = c(1e-4, 1e-3, 1e-2), 
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  coord_flip(ylim = c(1e-4, 1e-2)) +
  scale_x_continuous(breaks = rcts$idx, labels = rcts$Gene, expand = expansion(0)) +
  scale_color_manual(values = tclrs) +
  ylab('ddPCR expression') +
  facet_wrap(~'DEG validation')+
  xlab('Gene') +
  theme(plot.background = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(size = 11, color = 'white'),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        panel.grid.major = element_line(color = 'grey80', linetype = 'dashed'),
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank(),
        legend.position = 'bottom') -> pltt3

ggsave('f4_f.pdf', pltt3, height = 3.97, width = 2.12, bg = 'transparent')
