library(data.table)
library(tidyverse)
library(rtracklayer)
library(scales)
library(pals)
library(cowplot)
library(patchwork)
library(ggnewscale)
library(ggforce)
library(ggh4x)
library(ggpubr)
library(RcppCNPy)
library(gggenes)

tclrs <- c('#e45756', '#4c78a8')
bclr <- '#54a24b'

cgis <- list.files(file.path('..', 'data', 'cgi'), full.names = T) %>%
  setNames(., sub('.bed', '', basename(.))) %>%
  lapply(import)

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

gclrs <- c('#f28e2b', '#b07aa1')

list(BT245 = list(gn = 'hg38',
                  nm = 'This study\npHGG (K27M)\nK27M removal',
                  reg = GRanges('chr1', IRanges(78004000,78084000)),
                  bws = list('K27M' = 'BT245_K27M_H3K27me3.bw',
                             'K27M-KO' = 'BT245_KO_H3K27me3.bw')),
     G477 = list(gn = 'hg38',
                 nm = 'This study\npHGG (WT H3)\nK27M addition',
                 reg = GRanges('chr1', IRanges(181270000,181350000)),
                 bws = list('K27M-OE' = 'G477_K27M_H3K27me3.bw',
                            'WT H3' = 'G477_WT_H3K27me3.bw')),
     PFA = list(gn = 'hg38', nm = 'This study\nPFA ependymoma\nOxygenation',
                reg = GRanges('chr1', IRanges(152070000, 152150000)),
                bws = list('EZHIP high' = 'hypoxia_H3K27me3.bw',
                           'EZHIP low' = 'normoxia_H3K27me3.bw')),
     Diff = list(gn = 'hg38', 
                 nm = 'This study\nhiPSC\nNeural induction',
                 reg = GRanges('chr1', IRanges(26020000,26100000)),
                 bws = list('iPSC' = 'hiPSC_WT_H3K27me3.bw',
                            'NPC' = 'hNPC_WT_H3K27me3.bw')),
     mESC_Bonev = list(gn = 'mm10', 
                       nm = 'Bonev 2017\nmESC\nNeural induction',
                       reg = GRanges('chr1', IRanges(189300000,189380000)),
                       bws = list('mESC (serum)' = 'mESC_H3K27me3.bw',
                                  'mNPC' = 'mNPC_H3K27me3.bw')),
     mESC_Joshi = list(gn = 'mm10', 
                       nm = 'McLaughlin 2019\nmESC\nNaïve->primed',
                       reg = GRanges('chr1', IRanges(54883000,54963000)),
                       bws = list('mESC (serum)' = 'mESC_serum_H3K27me3.bw',
                                  'mESC (2i)' = 'mESC_2i_H3K27me3.bw'))) %>%
  Map(function(i, proj) {
    print(proj)
    nms <- strsplit(i$nm, '\n')[[1]]
    i$reg <- resize(i$reg, 2e5 + 1, fix = 'center')
    bws <- lapply(i$bws, function(f) {
      import.bw(file.path('..', 'data', 'bws', f), selection = BigWigSelection(i$reg))
    })
    bins <- tile(i$reg, 200)[[1]]
    seqlevels(bins) <- lapply(bws, function(x) seqlevels(x)) %>%
      Reduce(intersect, .) 
    xs <- mid(bins)
    scs <- lapply(bws, function(x) {
      seqlevels(x) <- seqlevels(bins)
      mcolAsRleList(x, 'score') %>%
        binnedAverage(bins, ., 'score', na.rm = T) %>%
        score() %>%
        tibble(score = .) %>%
        mutate(idx = xs)
    }) %>% bind_rows(.id = 'samp') %>%
      mutate(samp = fct_inorder(samp),
             proj = i$nm,
             t1 = nms[1], t2 = nms[2], t3 = nms[3])
    cgi <- cgis[[i$gn]] %>%
      subsetByOverlaps(i$reg) %>%
      resize(1, fix = 'center') %>%
      as_tibble() %>%
      dplyr::select(1:2) %>%
      {lapply(bws, function(x) .)} %>%
      bind_rows(.id = 'samp') %>%
      mutate(samp = fct_inorder(samp),
             proj = i$nm,
             t1 = nms[1], t2 = nms[2], t3 = nms[3])
    
    pg <-  subsetByOverlaps(txs[[i$gn]], i$reg) %>%
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
      ggplot(aes(xmin = start, xmax = end, y = strand, color = strand,
                 fill = strand, forward = forward)) +
      geom_gene_arrow(color = NA, data = ~subset(., type == 'mRNA')) +
      scale_y_discrete(breaks = c('-', '+'), labels = c('-', '+')) +
      scale_x_continuous(sprintf('%s:%s-%skb', seqnames(i$reg),
                                 formatC(as.integer(start(i$reg)/1e3), big.mark = ','), 
                                 formatC(as.integer(end(i$reg)/1e3), big.mark = ',')),
                         expand = expansion(0)) +
      coord_cartesian(xlim = range(scs$idx)) +
      ylab('Genes') +
      scale_fill_manual(values = gclrs) +
      scale_color_manual(values = gclrs) +
      theme(plot.background = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(size = 11, color = 'black'),
            axis.line = element_blank(),
            axis.ticks = element_blank(),
            legend.position = 'none',
            axis.title.y = element_blank(),
            axis.text.x = element_blank(),
            panel.grid = element_blank())
    
    p <- ggplot(scs, aes(x = idx, y = score)) +
      geom_vline(aes(xintercept = start), data = cgi, linetype = 'dashed',
                 inherit.aes = F, color = bclr, size = 1) +
      geom_line(aes(color = samp), size = .5) +
      geom_area(aes(fill = samp), alpha = .75) +
      
      geom_label(aes(label = samp, color = samp), x = -Inf, y = Inf, hjust = 0, vjust = 1,
                 data = distinct(scs, samp), label.size = NA, fill = '#ffffff99') +
      scale_x_continuous(sprintf('%s:%s-%skb', seqnames(i$reg),
                                 formatC(as.integer(start(i$reg)/1e3), big.mark = ','), 
                                 formatC(as.integer(end(i$reg)/1e3), big.mark = ',')),
                         expand = expansion(0)) +
      scale_color_manual(values = tclrs) +
      scale_fill_manual(values = tclrs) +
      facet_nested(samp ~ t1 + t2 + t3, scales = 'free',
                   strip = strip_nested(background_x = elem_list_rect(fill = c('grey0', 'grey20', 'grey40')))) +
      scale_y_continuous(expand = expansion(c(0, .2)), breaks = pretty_breaks(2)) +
      ylab('H3K27me3 CPM') +
      coord_cartesian(clip = 'off') +
      theme(legend.position = 'none',
            panel.grid = element_blank(),
            plot.background = element_blank(),
            panel.background = element_blank(),
            axis.text = element_text(color = 'black', size = 11),
            axis.line.y = element_line(color = 'black'),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            strip.background = element_rect(fill = 'black'),
            strip.text = element_text(color = 'white', size = 11),
            strip.background.y = element_blank(),
            strip.text.y = element_blank())
    
    p <- p + theme(axis.title.y = element_blank())
    
    plot_grid(p + theme(axis.title.x = element_blank(),
                        plot.margin = margin(1, 1, -1, 1)),
              pg + theme(plot.margin = margin(-1, 1, 1, 1)),
              nrow = 2, rel_heights = c(5, 1), align = 'v', axis = 'lr')
  }, ., names(.)) -> ps

plot_grid(plotlist = ps, nrow = 1) %>%
  ggsave('f1_a.pdf', ., width = 10.8, height = 3)

ps <- readRDS('../data/chip/H3K27me3.ssp.rds') %>%
  split(., .$proj) %>%
  lapply(function(x) {
    pd <-  x %>%
      mutate(mode = factor(mode, c('narrow', 'broad'))) %>%
      arrange(mode) %>%
      mutate(cond = fct_inorder(as.character(cond)),
             x = as.numeric(cond))
    
    cns <- pd %>% 
      dplyr::select(rep, fcs, cond) %>%
      group_by(rep, cond) %>%
      summarise(fcs = mean(fcs), .groups = 'drop') %>%
      pivot_wider(names_from = 'cond', values_from = 'fcs') %>%
      dplyr::rename(y = 2, yend = 3) %>%
      mutate(x = 1, xend = 2)
    
    
    ggplot(pd, aes(x = x, y = fcs, group = rep, shape = rep)) +
      geom_link(aes(x = x, xend = xend, y = y, yend = yend, color = stat(index)),
                data = cns, inherit.aes = F, show.legend = F) +
      scale_color_gradient(low = tclrs[1], high = tclrs[2]) +
      new_scale_color() +
      geom_point(aes(color = cond), size = 3) +
      scale_color_manual(values = tclrs, guide = guide_none()) + 
      scale_x_continuous(breaks = 1:2, labels = levels(pd$cond), 
                         expand = expansion(c(.25))) +
      coord_cartesian(clip  = 'off') +
      facet_wrap(~proj) +
      ylab('H3K27me3 confinement') +
      theme(plot.background = element_blank(),
            panel.background = element_blank(),
            panel.grid = element_blank(),
            axis.title.x = element_blank(),
            strip.background = element_rect(fill = 'black'),
            strip.text = element_text(color = 'white', size = 11),
            axis.text = element_text(size = 11, color = 'black'),
            panel.grid.major = element_line(color = 'grey75', linetype = 'dashed'),
            axis.ticks.y = element_blank(),
            #axis.title.y = element_blank(),
            legend.title = element_blank(),
            legend.background = element_rect(colour = NA, fill = '#ffffff00'),
            legend.key = element_blank(),
            legend.position = 'none',
            legend.justification = c(.1,-.01),
            legend.text = element_text(size = 11),
            axis.line.x = element_line(color = 'black'))
  })

ps[c('pHGG (K27M)', 'pHGG (WT H3)', 'Ependymoma', 'hDiff', 'mDiff', 'mESC')] %>%
  {wrap_plots(., nrow = 1) & theme(plot.background = element_blank())} %>%
  ggsave('f1_b.pdf', ., height = 2.5, width = 10.9)


t1k <- list(BT245 = c(K27M = 'BT245_K27M_H3K27me3.mat.gz',
                       KO = 'BT245_KO_H3K27me3.mat.gz'),
             DIPGXIII = c(K27M = 'DIPGXIII_K27M_H3K27me3.mat.gz',
                       KO = 'DIPGXIII_KO_H3K27me3.mat.gz'),
             HSJ019 = c(K27M = 'HSJ019_K27M_H3K27me3.mat.gz',
                       KO = 'HSJ019_KO_H3K27me3.mat.gz'),
             G477 = c(K27M = 'G477_K27M_H3K27me3.mat.gz',
                      'WT H3' = 'G477_WT_H3K27me3.mat.gz'),
             PFA = c('EZHIP high' = 'hypoxia_H3K27me3.mat.gz',
                     'EZHIP low' = 'normoxia_H3K27me3.mat.gz'),
             hDiff = c(iPSC = 'hiPSC_WT_H3K27me3.mat.gz',
                       NPC = 'hNPC_WT_H3K27me3.mat.gz'),
             mDiff = c(mESC = 'mESC_H3K27me3.mat.gz',
                       mNPC = 'mNPC_H3K27me3.mat.gz'),
             mESC = c(serum = 'mESC_serum_H3K27me3.mat.gz',
                      '2i' = 'mESC_2i_H3K27me3.mat.gz'),
             NSD1KO = c(WT = 'mESC_Par_H3K27me3.mat.gz',
                        NSD1KO = 'mESC_NSD1KO_H3K27me3.mat.gz'),
             NSD1het = c(WT = 'hiPSC_WT_H3K27me3.mat.gz',
                         'NSD1+/-' = 'hiPSC_NSD1het_H3K27me3.mat.gz'),
             BAP1KO = c(WT = 'mESC_WT_H3K27me3.mat.gz',
                        BAP1KO = 'mESC_BAP1KO_H3K27me3.mat.gz'),
             GCB = c(DKO = 'GCB_H1DKO_H3K27me3.mat.gz',
                     WT = 'GCB_WT_H3K27me3.mat.gz')) %>%
  pblapply(function(x) {
    i <- lapply(x, function(y) {
      p <- file.path('..', 'data', 'aggr', y)
      cgi <- fread(p, skip = 1, select = 1:6)$V4 %>%
        grepl('^EN', .) %>%
        {which(!.)}
      fread(p, skip = 1, drop = 1:6) %>%
        .[,950:1050] %>%
        rowMeans(na.rm = T) %>%
        order(decreasing = T) %>%
        intersect(cgi) %>%
        head(1000)
    }) %>%
      unlist() %>%
      unique()
    
    fread(file.path('..', 'data', 'aggr', x[1]), skip = 1, select = 1:3) %>%
      .[i,] 
  })


aggr <- list(BT245 = c(K27M = 'BT245_K27M_H3K27me3.mat.gz',
                       KO = 'BT245_KO_H3K27me3.mat.gz'),
             DIPGXIII = c(K27M = 'DIPGXIII_K27M_H3K27me3.mat.gz',
                          KO = 'DIPGXIII_KO_H3K27me3.mat.gz'),
             HSJ019 = c(K27M = 'HSJ019_K27M_H3K27me3.mat.gz',
                        KO = 'HSJ019_KO_H3K27me3.mat.gz'),
             G477 = c(K27M = 'G477_K27M_H3K27me3.mat.gz',
                      'WT H3' = 'G477_WT_H3K27me3.mat.gz'),
             PFA = c('EZHIP high' = 'hypoxia_H3K27me3.mat.gz',
                     'EZHIP low' = 'normoxia_H3K27me3.mat.gz'),
             hDiff = c(iPSC = 'hiPSC_WT_H3K27me3.mat.gz',
                       NPC = 'hNPC_WT_H3K27me3.mat.gz'),
             mDiff = c(mESC = 'mESC_H3K27me3.mat.gz',
                       mNPC = 'mNPC_H3K27me3.mat.gz'),
             mESC = c(serum = 'mESC_serum_H3K27me3.mat.gz',
                      '2i' = 'mESC_2i_H3K27me3.mat.gz'),
             NSD1KO = c(WT = 'mESC_Par_H3K27me3.mat.gz',
                        NSD1KO = 'mESC_NSD1KO_H3K27me3.mat.gz'),
             NSD1het = c(WT = 'hiPSC_WT_H3K27me3.mat.gz',
                         'NSD1+/-' = 'hiPSC_NSD1het_H3K27me3.mat.gz'),
             BAP1KO = c(WT = 'mESC_WT_H3K27me3.mat.gz',
                        BAP1KO = 'mESC_BAP1KO_H3K27me3.mat.gz'),
             GCB = c(DKO = 'GCB_H1DKO_H3K27me3.mat.gz',
                     WT = 'GCB_WT_H3K27me3.mat.gz')) %>%
  lapply(function(x) {
    i <- lapply(x, function(y) {
      p <- file.path('..', 'data', 'aggr', y)
      cgi <- fread(p, skip = 1, select = 1:6)$V4 %>%
        grepl('^EN', .) %>%
        {which(!.)}
      fread(p, skip = 1, drop = 1:6) %>%
        .[,950:1050] %>%
        rowMeans(na.rm = T) %>%
        order(decreasing = T) %>%
        intersect(cgi) %>%
        head(1000)
    }) %>%
      unlist() %>%
      unique()
    lapply(x, function(y) {
      p <- file.path('..', 'data', 'aggr', y)
      fread(p, skip = 1, drop = 1:6) %>%
        .[i,] %>%
        as.matrix() %>%
        {.[. > quantile(c(.), .99999, na.rm = T)] <- NA; .} %>%
        {tibble(mu = colMeans(., na.rm = T),
                std = colSds(., na.rm = T),
                num = colSums(is.finite(.)))} %>%
        mutate(idx = 1:n())
    }) %>%
      bind_rows(.id = 'cond') %>%
      mutate(cond = fct_inorder(cond),
             i = as.numeric(cond))
  })



aggr %>% 
  bind_rows(.id = 'proj') %>%
  filter(proj %in% c('BT245', 'G477', 'PFA','hDiff','mDiff','mESC')) %>%
  mutate(i = as.character(i),
         se = std / sqrt(num),
         proj = fct_inorder(proj)) %>%
  ggplot(aes(x = idx, y = mu, color = i, fill = i, ymin = mu - se, ymax = mu + se)) +
  geom_ribbon(color = NA, alpha = .5) +
  geom_line() +
  facet_wrap(~proj, scales = 'free', nrow = 1) +
  facetted_pos_scales(x = list(
    proj == 'BT245' ~ scale_x_continuous(breaks = c(900, 1000, 1100),
                                         labels = c('-100kb', 'CGI', '+100kb')),
    T ~ scale_x_continuous(breaks = c(900, 1000, 1100), labels = c('','',''))
  )) +
  scale_color_manual(values = tclrs) +
  scale_fill_manual(values = tclrs) +
  coord_cartesian(xlim = c(900, 1100)) +
  ylab('H3K27me3\nCPM') +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white', size = 11),
        axis.text = element_text(size = 11, color = 'black'),
        panel.grid.major = element_line(color = 'grey75', linetype = 'dashed'),
        axis.ticks.y = element_blank(),
        legend.position = 'none',
        axis.line.x = element_line(color = 'black')) -> p

ggsave('f1_c.pdf', p, height = 1.5, width = 10.9)

mid_rescaler <- function(mid = 0) {
  function(x, to = c(0, 1), from = range(x, na.rm = TRUE)) {
    scales::rescale_mid(x, to, from, mid)
  }
}

ps <- readRDS('../data/pup/H3K27me3_CGI.rds') %>%
  separate(samp, c('proj', 'cond'), '_') %>%
  filter(cond != 'NSD1het') %>%
  mutate(cond = case_when(proj %in% c('hiPSC', 'hNPC') ~ proj,
                          T ~ cond),
         proj = case_when(proj %in% c('hiPSC', 'hNPC') ~ 'hDiff',
                          proj != 'mESC' ~ proj,
                          cond %in% c('2i', 'serum') ~ 'mESC',
                          cond %in% c('Par','NSD1KO') ~ 'NSD1KO',
                          cond %in% c('WT', 'BAP1KO') ~ 'BAP1KO')) %>% 
  split(., .$proj) %>%
  lapply(function(x) {
    dd <- x %>% mutate(cond = fct_inorder(cond))
    ann <- dd %>%
      filter(between(x, 11, 11) & between(y, 11, 11)) %>%
      group_by(cond, proj) %>%
      summarise(z = sprintf('%.2f', mean(z, na.rm = T))) %>%
      ungroup() 
    ll <- quantile(dd$z, .9999) %>% unname()
    p <- ggplot(dd, aes(x = x, y = y, fill = z)) +
      geom_raster() +
      geom_label(aes(x = Inf, y = Inf, hjust =1, vjust = 1, label = z), data = ann, inherit.aes = F,
                 alpha = .8, label.size = NA) +
      scale_fill_gradientn('O/E', colors = coolwarm(),
                           guide = guide_colorbar(barheight = .5,
                                                  barwidth = 3.5,
                                                  title.vjust = 1),
                           #breaks = scales::pretty_breaks(3), 
                           limits = c(1/ll, ll),
                           breaks = c(ceiling(10/ll) / 10, 1, floor(ll * 10)/10),
                           oob = scales::squish,
                           trans = 'log2') +
      facet_nested(.~ proj + cond,
                   strip = strip_nested(background_x = elem_list_rect(fill = c('black', tclrs)))) +
      coord_cartesian(expand = F, clip = 'off') +
      theme(strip.background = element_rect(fill = 'black'),
            strip.text = element_text(color = 'white', size = 11),
            panel.background = element_blank(),
            plot.background = element_blank(),
            axis.title = element_blank(),
            axis.text = element_blank(),
            legend.position = 'bottom',
            legend.text = element_text(size = 11),
            legend.background = element_blank(),
            legend.key = element_blank(),
            axis.ticks = element_blank()) 
  })

{wrap_plots(ps[c('BT245','G477','PFA','hDiff','mDiff','mESC')],nrow=1) &
    theme(plot.background = element_blank())} %>%
  ggsave('f1_d.pdf', ., height = 2.8, width=10.9, bg='transparent')

readRDS('../data/pup/H3K27me3_CGI_chromosight.rds') %>%
  Map(function(dd, nm) {
    print(nm)
    ann <- compare_means(score ~ cond, dd) %>%
      mutate(y.position = max(dd$score) + (diff(range(dd$score))) * .05)
    p <- ggplot(dd, aes(x = cond, y = score)) +
      geom_violin(aes(fill = cond), alpha = .25,color = NA) +
      geom_boxplot(aes(color = cond, fill = cond), notch = T, width = .3) +
      stat_summary(geom = "crossbar", width=0.15, fatten=0, color="white", 
                   fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
      stat_pvalue_manual(data = ann, inherit.aes = F, label = 'p.signif',
                         coord.flip = F, tip.length = 0) +
      facet_wrap(~ proj, scales = 'free') +
      ylab('CGI interactivity') +
      scale_y_continuous(expand = expansion(c(.01,.08))) +
      #coord_flip() +
      scale_color_manual(values = tclrs) +
      scale_fill_manual(values = tclrs) +
      theme(legend.position = 'none',
            axis.text = element_text(color = 'black', size = 11),
            legend.title = element_blank(),
            legend.key = element_blank(),
            legend.background = element_blank(),
            strip.background = element_rect(fill = 'black'),
            strip.text = element_text(color = 'white', size = 11),
            panel.grid = element_blank(),
            panel.background = element_rect(fill = NA, color = 'black', size = 1),
            plot.background = element_blank(),
            axis.text.x = element_text(angle = 30, hjust= 1),
            panel.grid.major.y = element_line(color = 'grey85', linetype = 'dashed'),
            axis.title.x = element_blank()) 
    p
  }, ., names(.))

{wrap_plots(ps[c('BT245','G477','PFA','hDiff','mDiff','mESC')], nrow = 1) &
    theme(plot.background = element_blank())} %>%
  ggsave('f1_d_p.pdf', ., width = 10.9, height = 3.4)


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

res <- 10000
l <- c(174e6, 178e6)
reg <- GRanges('chr2', IRanges(l[1], l[2]))
crds <- seq(l[1], l[2], res)

mats <- c('K27M' = 'HSJ019_K27M', 'K27M-KO' = 'HSJ019_KO') %>%
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
  filter(between(y, 1, 120)) %>%
  filter(between(x, 101, 500)) %>%
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
  scale_y_continuous(breaks = c(50, 100),
                     labels = c('1mb', '2mb')) +
  scale_x_continuous(breaks = c(200,400)) +
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

bws <- c('K27M' = 'HSJ019_K27M_H3K27me3.bw',
         'K27M-KO' = 'HSJ019_KO_H3K27me3.bw') %>%
  lapply(function(x) {
    import.bw(sprintf('../data/bws/%s', x),
              selection = BigWigSelection(reg))
  })

bins <- tile(reg, 500)[[1]]
seqlevels(bins) <- seqlevels(bws[[1]])
xs <- mid(bins)

rx <- readRDS('../../rx.rds') %>%
  {setNames(.$r, .$samp)} %>%
  {list('K27M' = .['HSJ019-C7-2-Rx-C_H3K27me3'] / .['HSJ019-C7-2-Rx_Input'],
        'K27M-KO' = .['ASH73_HSJ019-ko-c10p29-Rx_H3K27me3-C'] / .['ASH51_HSJ019-ko-c10p29-Rx_Input'])} %>%
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
  ggplot(aes(x = idx, y = score, fill = samp)) +
  geom_line(aes(color = samp), size = .5) +
  geom_area(aes(fill = samp), alpha = .75) +
  coord_cartesian(xlim = l, ylim = c(0, .3)) +
  scale_y_continuous('H3K27me3\n(ChIP-Rx)',
                     breaks = c(0, .3)) +
  facet_wrap(~samp) +
  scale_fill_manual(values = tclrs) +
  scale_color_manual(values = tclrs) +
  scale_x_continuous(expand = expansion(0)) +
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
        axis.line.y = element_line(color = 'black'))

gclrs <- c('#f28e2b', '#b07aa1')

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
  geom_gene_arrow(data = ~subset(., type == 'mRNA')) +
  scale_alpha_manual(values = c('UTR' = .5, 'not' = 1)) +
  scale_x_continuous(seqnames(reg), limits = l,
                     breaks = c(175e6, 177e6),
                     labels = c('175mb', '177mb'),
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
        legend.position = 'none') -> p3

{ wrap_plots(p1,p2,p3,ncol=1) &
    theme(plot.background = element_blank()) } %>%
  ggsave('f1_f.pdf', ., height = 4, width = 7)
