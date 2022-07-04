library(data.table)
library(tidyverse)
library(pals)
library(rtracklayer)
library(ggh4x)
library(patchwork)
library(readxl)
library(igraph)
library(tidygraph)
library(ggraph)
library(scattermore)

tclrs <- c('#e45756', '#4c78a8')
bclr <- '#54a24b'

cgis <- list.files('../data/cgi', full.names = T) %>%
  setNames(., sub('.bed', '', basename(.))) %>%
  lapply(import)

list(DIPGXIII = list(gn = 'hg38', 
                   nm = 'K27M glioma (DIPGXIII, pontine)', 
                   reg = GRanges('chr1',  IRanges(9040000,9340000)),
                   cnds = c('K27M','K27M-KO'),
                   bws = list(H3K27me3 = c('DIPGXIII_K27M_H3K27me3.bw',
                                           'DIPGXIII_KO_H3K27me3.bw'),
                              RING1B = c('DIPGXIII_K27M_RING1B.bw',
                                         'DIPGXIII_KO_RING1B.bw'),
                              CBX2 = c('DIPGXIII_K27M_CBX2.bw',
                                       'DIPGXIII_KO_CBX2.bw'),
                              H2Aub = c('DIPGXIII_K27M_H2AK119ub.bw',
                                        'DIPGXIII_KO_H2AK119ub.bw'))),
     HSJ019 = list(gn = 'hg38',
                   nm = 'K27M glioma (HSJ019, thalamic)',
                   reg = GRanges('chr1',  IRanges(114400000,117400000)),
                   cnds = c('K27M','K27M-KO'),
                   bws = list(H3K27me3 = c('HSJ019_K27M_H3K27me3.bw',
                                           'HSJ019_KO_H3K27me3.bw'),
                              RING1B = c('HSJ019_K27M_RING1B.bw',
                                         'HSJ019_KO_RING1B.bw'),
                              CBX2 = c('HSJ019_K27M_CBX2.bw',
                                       'HSJ019_KO_CBX2.bw'),
                              H2Aub = c(NA, NA))),
     G477 = list(gn = 'hg38',
                 nm = 'WT H3 glioma (G477, cortical)', 
                 reg = GRanges('chr1', IRanges(47138000,47308000)),
                 cnds = c('K27M-OE','WT H3'),
                 bws = list(H3K27me3 = c('G477_K27M_H3K27me3.bw',
                                         'G477_WT_H3K27me3.bw'),
                            RING1B = c('G477_K27M_RING1B.bw',
                                       'G477_WT_RING1B.bw'),
                            CBX2 = c('G477_K27M_CBX2.bw',
                                     'G477_WT_CBX2.bw'),
                            H2Aub = c(NA, NA))))  %>%
  Map(function(i, proj) {
    print(proj)
    bws <- lapply(i$bws, function(m) {
      setNames(m, i$cnds) %>%
        lapply(function(f) {
          if (!is.na(f)) {
            print(f)
            import.bw(file.path('../data/bws', f), selection = BigWigSelection(i$reg))
          } else {
            NA
          }
        })
    })
    bins <- tile(i$reg, 400)[[1]]
    seqlevels(bins) <- lapply(bws, function(x) {
      lapply(x[sapply(x, is, 'GenomicRanges')], seqlevels)
    }) %>%
      unlist(recursive = F) %>%
      Reduce(intersect, .) 
    xs <- mid(bins)
    scs <- lapply(bws, function(x) {
      lapply(x, function(y) {
        if (!is(y, 'GenomicRanges')) {
          tibble(score = 0, idx = xs)
        } else {
          seqlevels(y) <- seqlevels(bins)
          mcolAsRleList(y, 'score') %>%
            binnedAverage(bins, ., 'score') %>%
            score() %>%
            tibble(score = .) %>%
            mutate(idx = xs)
        }
      }) %>% bind_rows(.id = 'samp')
    }) %>% bind_rows(.id = 'mark') %>%
      mutate(samp = factor(samp, i$cnds),
             mark = factor(mark, names(i$bws)),
             proj = i$nm)
    cgi <- cgis[[i$gn]] %>%
      subsetByOverlaps(i$reg) %>%
      as_tibble() %>%
      dplyr::select(1:3) %>%
      {lapply(bws, function(x) lapply(x, function(y) .) %>% bind_rows(.id = 'samp'))} %>%
      bind_rows(.id = 'mark') %>%
      mutate(samp = factor(samp, i$cnds),
             mark = factor(mark, names(i$bws)),
             proj = i$nm)
    
    if (!('yl' %in% names(i))) {
      i$yl <- scs %>% group_by(mark, samp) %>% summarise(v = max(score)) %>%
        arrange(mark, samp) %>% pull(v)
    }
    
    p <- ggplot(scs, aes(x = idx, y = score)) +
      #geom_rect(aes(xmin = start, xmax = end), data = cgi, ymin = -Inf, ymax = Inf, 
      #          inherit.aes = F, fill = bclr, alpha = .5) +
      geom_line(aes(color = samp), size = .5) +
      geom_area(aes(fill = samp), alpha = .5) +
      geom_label(aes(label = samp, color = samp), x = -Inf, y = Inf, hjust = 0, vjust = 1,
                 data = distinct(scs, samp, mark), label.size = NA, fill = '#ffffff99') +
      scale_x_continuous(sprintf('%s:%s-%skb', seqnames(i$reg),
                                 formatC(as.integer(start(i$reg)/1e3), big.mark = ','), 
                                 formatC(as.integer(end(i$reg)/1e3), big.mark = ',')),
                         expand = expansion(0)) +
      facet_nested(mark + samp ~ proj, scales = 'free')  +
      scale_color_manual(values = tclrs) +
      scale_fill_manual(values = tclrs) +
      facetted_pos_scales(y = list(
        scale_y_continuous(expand = expansion(c(0, .2)),
                           breaks = scales::pretty_breaks(2),
                           limits = c(0, i$yl[1])),
        scale_y_continuous(expand = expansion(c(0, .2)),
                           breaks = scales::pretty_breaks(2),
                           limits = c(0, i$yl[2])),
        scale_y_continuous(expand = expansion(c(0, .2)),
                           breaks = scales::pretty_breaks(2),
                           limits = c(0, i$yl[3])),
        scale_y_continuous(expand = expansion(c(0, .2)),
                           breaks = scales::pretty_breaks(2),
                           limits = c(0, i$yl[4])),
        scale_y_continuous(expand = expansion(c(0, .2)),
                           breaks = scales::pretty_breaks(2),
                           limits = c(0, i$yl[5])),
        scale_y_continuous(expand = expansion(c(0, .2)),
                           breaks = scales::pretty_breaks(2),
                           limits = c(0, i$yl[6])),
        scale_y_continuous(expand = expansion(c(0, .2)),
                           breaks = scales::pretty_breaks(2),
                           limits = c(0, i$yl[7])),
        scale_y_continuous(expand = expansion(c(0, .2)),
                           breaks = scales::pretty_breaks(2),
                           limits = c(0, i$yl[8]))
      )) +
      ylab('CPM') +
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
            strip.text = element_text(color = 'white', size = 11))
    #if (proj != 'DIPG13') {
    #  p <- p + theme(axis.title.y = element_blank())
    #}
    p
  }, ., names(.)) -> ps

{wrap_plots(ps, ncol = 1) &
    theme(plot.background = element_blank()) } %>%
  ggsave('sf6_adf.pdf', ., width = 3.7, height = 13, bg = 'transparent')


i <- fread('../data/aggr/lfc/BT245_K27M_H3K27me3.mat.gz', skip = 1, select = 1:6) %>%
  {which(!grepl('ENS', .$V4))}


list(DIPGXIII = list(H3K27me3 = c(K27M = 'DIPGXIII_K27M_H3K27me3.mat.gz',
                                  KO = 'DIPGXIII_KO_H3K27me3.mat.gz'),
                     RING1B = c(K27M = 'DIPGXIII_K27M_RING1B.mat.gz',
                                KO = 'DIPGXIII_KO_RING1B.mat.gz'),
                     CBX2 = c(K27M = 'DIPGXIII_K27M_CBX2.mat.gz',
                              KO = 'DIPGXIII_KO_CBX2.mat.gz')),
     HSJ019 = list(H3K27me3 = c(K27M = 'HSJ019_K27M_H3K27me3.mat.gz',
                                KO = 'HSJ019_KO_H3K27me3.mat.gz'),
                   RING1B = c(K27M = 'HSJ019_K27M_RING1B.mat.gz',
                              KO = 'HSJ019_KO_RING1B.mat.gz'),
                   CBX2 = c(K27M = 'HSJ019_K27M_CBX2.mat.gz',
                            KO = 'HSJ019_KO_CBX2.mat.gz')),
     G477 = list(H3K27me3 = c(K27M = 'G477_K27M_H3K27me3.mat.gz',
                              KO = 'G477_WT_H3K27me3.mat.gz'),
                 RING1B = c(K27M = 'G477_K27M_RING1B.mat.gz',
                            KO = 'G477_WT_RING1B.mat.gz'),
                 CBX2 = c(K27M = 'G477_K27M_CBX2.mat.gz',
                          KO = 'G477_WT_CBX2.mat.gz'))) %>%
  lapply(function(x) {
    lapply(x, function(y) {
      lapply(y, function(z) {
        fread(file.path('../data/aggr/lfc', z), skip = 1, drop = 1:6) %>%
          as.matrix() %>%
          .[i,999:1001] %>%
          rowMeans(na.rm = T) %>%
          tibble(v = .) %>%
          mutate(idx = i) 
      }) %>% bind_rows(.id = 'cond')
    }) %>% bind_rows(.id = 'mark')
  }) %>% bind_rows(.id = 'line') -> agg

readRDS('../data/chip/rx.rds') %>%
  separate(samp, c('line', 'cond', NA), '_') %>%
  mutate(rx = log2(rx),
         cond = ifelse(cond == 'WT','KO', cond)) %>%
  split(.$line) %>%
  lapply(function(x) setNames(x$rx, x$cond)) -> rx


split(agg, agg$line) %>%
  lapply(function(x) {
    x %>%
      dplyr::select(-line) %>%
      pivot_wider(names_from = 'cond', values_from = 'v') %>%
      mutate(dif = case_when(mark == 'H3K27me3' ~ (K27M + rx[[x$line[1]]]['K27M']) -
                               (KO + rx[[x$line[1]]]['KO']),
                             T ~ K27M - KO)) %>%
      dplyr::select(-K27M, -KO) %>%
      pivot_wider(names_from = 'mark', values_from = 'dif') %>%
      ggplot(aes(x = H3K27me3, y = RING1B, color = CBX2)) +
      geom_hline(yintercept = 0, color = 'grey50') +
      geom_vline(xintercept = 0, color = 'grey50') +
      geom_scattermore(alpha = .5, pointsize = 4) +
      #annotate('label', x = Inf, y = Inf, hjust = 1.05, vjust = 1.2,
      #         label = sprintf('Pearson\'s r = %.2f', cor(dd[[m1]] - dd[[m4]], dd[[m2]] - dd[[m5]]))) +
      scale_color_gradientn(expression(Delta*'CBX2'), colors = rev(brewer.rdylbu(20)), 
                            guide = guide_colorbar(barheight = 3),
                            limits = c(-3,3), oob = scales::squish,
                            breaks = scales::pretty_breaks(3)) +
      labs(x = expression(atop(Delta*'H3K27me3,'~log[2]~'K27M/KO','ChIP-Rx normalized')), 
           y = expression(Delta*'Ring1B,'~log[2]~'K27M/KO')) +
      facet_grid(~'Differential enrichment in CGIs') +
      theme(plot.background = element_blank(),
            panel.background = element_rect(color= 'black', size = 1, fill = NA),
            panel.grid = element_blank(),
            panel.grid.major = element_line(color = 'grey75', linetype = 'dashed'),
            legend.position = c(0, 1),
            legend.justification = c(0,1),
            legend.background = element_blank(),
            strip.text = element_text(color = 'white', size = 11),
            strip.background = element_rect(fill = 'black'),
            axis.text = element_text(size = 11, color = 'black')) 
  }) -> ps

{wrap_plots(ps, ncol = 1) &
    theme(plot.background = element_blank())} %>%
  ggsave('sf6_beg.pdf', ., height = 9.5, width = 3.7, bg = 'transparent')

list(H3K27me3 = c(K27M = 'DIPGXIII_K27M_H3K27me3.mat.gz',
                  KO = 'DIPGXIII_KO_H3K27me3.mat.gz'),
     RING1B = c(K27M = 'DIPGXIII_K27M_RING1B.mat.gz',
                KO = 'DIPGXIII_KO_RING1B.mat.gz'),
     CBX2 = c(K27M = 'DIPGXIII_K27M_CBX2.mat.gz',
              KO = 'DIPGXIII_KO_CBX2.mat.gz'),
     H2Aub = c(K27M = 'DIPGXIII_K27M_H2Aub.mat.gz',
               KO = 'DIPGXIII_KO_H2Aub.mat.gz')) %>%
  lapply(function(x) {
    lapply(x, function(y) {
      fread(file.path('../data/aggr/lfc', y), skip = 1, drop = 1:6) %>%
        as.matrix() %>%
        .[i,999:1001] %>%
        rowMeans(na.rm = T) %>%
        tibble(v = .) %>%
        mutate(idx = i) 
    }) %>% bind_rows(.id = 'cond')
  }) %>% bind_rows(.id = 'mark') -> agg2

mclrs <- setNames(tableau20()[seq(5,11,2)],c('H2AK119ub', 'CBX2', 'Ring1B', 'H3K27me3'))
names(mclrs) <- paste0('\u0394', names(mclrs))
gcor <- agg2 %>%
  pivot_wider(names_from = 'cond', values_from = 'v') %>%
  mutate(dif = K27M - KO) %>%
  dplyr::select(-K27M, -KO) %>%
  pivot_wider(names_from = 'mark', values_from = 'dif') %>%
  dplyr::select(-idx) %>%
  {names(.)<-paste0('\u0394', names(.)); .} %>%
  cor() %>%
  as.data.frame() %>%
  rownames_to_column('x') %>%
  pivot_longer(-x, names_to = 'y', values_to = 'r') %>%
  graph_from_data_frame(directed = F)

plt2 <- ggraph(gcor, layout = 'linear', circular = T) +
  geom_edge_link(aes(edge_alpha = abs(r), edge_width = abs(r), color = r)) +
  guides(edge_alpha = "none", edge_width = "none") +
  scale_edge_colour_gradientn('Pearson\'s r', limits = c(-1, 1), colors = coolwarm(),
                              breaks = c(-1,0,1),
                              guide = guide_edge_colorbar(barheight = .5,
                                                          barwidth = 3,
                                                          title.vjust = 1,
                                                          title.hjust = 9,
                                                          title.position = 'top')) +
  geom_node_point(aes(color = name), size = 5) +
  geom_node_label(aes(label = name, color = name), repel = TRUE) +
  facet_grid(.~'Correlation of \u0394enrichment') +
  #coord_cartesian(clip = 'off') +
  scale_color_manual(values = mclrs,guide = guide_none()) +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        legend.background = element_blank(),
        legend.position = c(0.05,.5),
        legend.justification = c(0,.5),
        legend.direction = 'horizontal',
        strip.text = element_text(color = 'white', size = 11),
        strip.background = element_rect(fill = 'black')) 

ggsave('sf6_c.pdf', plt2, height = 2.5, width = 3, bg = 'transparent', device = cairo_pdf)


list(H3K27me3 = c(K27M = 'BT245_K27M_H3K27me3.mat.gz',
                  KO = 'BT245_KO_H3K27me3.mat.gz'),
     RING1B = c(K27M = 'BT245_K27M_RING1B.mat.gz',
                KO = 'BT245_KO_RING1B.mat.gz'),
     CBX2 = c(K27M = 'BT245_K27M_CBX2.mat.gz',
              KO = 'BT245_KO_CBX2.mat.gz'),
     H2Aub = c(K27M = 'BT245_K27M_H2Aub.mat.gz',
               KO = 'BT245_KO_H2Aub.mat.gz')) %>%
  lapply(function(x) {
    lapply(x, function(y) {
      fread(file.path('../data/aggr/lfc', y), skip = 1, drop = 1:6) %>%
        as.matrix() %>%
        .[i,999:1001] %>%
        rowMeans(na.rm = T) %>%
        tibble(v = .) %>%
        mutate(idx = i) 
    }) %>% bind_rows(.id = 'cond')
  }) %>% bind_rows(.id = 'mark') -> agg3

list(DIPGXIII = agg2, BT245 = agg3) %>%
  lapply(function(x) {
    x %>%
      pivot_wider(names_from = 'cond', values_from = 'v') %>%
      split(., .$mark) %>%
      lapply(function(y) {
        y %>%
          dplyr::select(-mark, -idx) %>%
          na.omit() %>%
          {cor(.[[1]], .[[2]])} %>%
          tibble(r = .)
      }) %>%
      bind_rows(.id = 'mark') 
  }) %>%
  bind_rows(.id = 'line') %>%
  ggplot(aes(x = mark, y = r)) +
  geom_point(size = 5) +
  ylab(expression('Spearman\'s'~rho)) +
  facet_grid(line ~ .) +
  scale_y_continuous('Inter-condition correlation across genome-wide 10kb windows',
                     expand = expansion(.1)) +
  scale_x_discrete(labels = function(x) paste0('\u0394', x)) +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black'),
        axis.text = element_text(color = 'black', size = 11),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white', size = 11),
        axis.title.x = element_blank(),
        panel.grid.major = element_line(color = 'grey75', linetype = 'dashed'),
        axis.text.x = element_text(angle = 45, hjust = 1)) -> p

ggsave('sf6_h.pdf',p,height=5,width=3, device = cairo_pdf, bg = 'transparent')


