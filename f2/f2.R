library(data.table)
library(tidyverse)
library(rtracklayer)
library(pals)
library(ggh4x)
library(patchwork)
library(InteractionSet)
library(gghalves)
library(scattermore)
library(ggpubr)
library(readxl)
library(igraph)
library(tidygraph)
library(ggraph)
library(DiffBind)
library(tximport)
library(ggnewscale)
library(gghalves)
library(eulerr)
library(RcppCNPy)
library(pals)
library(scattermore)

cgis <- list.files('../data/cgi', full.names = T) %>%
  setNames(., sub('.bed', '', basename(.))) %>%
  lapply(import)



getWhisks <- function(x) {
  x <- as.numeric(x)
  qs <- quantile(x, c(0.25, 0.75), na.rm = T)
  data.frame(lower = qs[1], upper = qs[2], middle = median(x, na.rm = T),
             ymin = min(x[x >= (qs[1] - 1.5 * diff(qs))], na.rm = T),
             ymax = max(x[x <= (qs[2] + 1.5 * diff(qs))], na.rm = T)) %>%
    mutate(notchlower = middle - 1.58 * diff(qs)/sqrt(length(x)),
           notchupper = middle + 1.58 * diff(qs)/sqrt(length(x)))
}


pks <- readRDS('../data/pks/PRC.rds')
clrs <- c('#f28e2b', '#59a14f', '#b07aa1', '#9c755f')

list(BT245 = list(fs = c(K27M = 'BT245_K27M',
                         `K27M-KO` = 'BT245_KO'),
                  nm = 'K27M glioma',
                  gn = 'hg38'),
     hiPSC = list(fs = c(`WT` = 'hiPSC_WT',
                         `NSD1+/-` = 'hiPSC_NSD1het'),
                  nm = 'hiPSC',
                  gn = 'hg38'),
     mESC = list(fs = c(`WT` = 'mESC_Par',
                        `NSD1-KO` = 'mESC_NSD1KO'),
                 nm = 'mESC',
                 gn = 'mm10')) %>%
  Map(function(i, proj) {
    print(proj)
    dd <- lapply(i$fs, function(x) fread(sprintf('../data/prc_csight/%s.tsv', x))) %>%
      {merge(.[[1]], .[[2]], 
             by = c('chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'))} %>%
      mutate(d = score.x - score.y) %>%
      filter(is.finite(d))
    
    ii <- list(dd[,1:3],  dd[,4:6]) %>%
      lapply(function(y) {
        y %>%
          `colnames<-`(c('chr', 'start', 'end')) %>%
          makeGRangesFromDataFrame()
      }) %>%
      {GInteractions(.[[1]], .[[2]], mode = 'strict')}
    
    u <- GRangesList(pks[[proj]]) %>% unlist() %>% reduce() %>% subsetByOverlaps(cgis[[i$gn]])
    p <- lapply(pks[[proj]], function(x) overlapsAny(u, x))
    
    vv <- list(All_three = p$CBX2 & p$RING1B & p$H3K27me3,
               Ring1B_only = p$RING1B & (!p$CBX2) &  (!p$H3K27me3),
               #CBX2_only = p$CBX2 & (!p$RING1B) &  (!p$H3K27me3),
               H3K27me3_only = p$H3K27me3 & (!p$CBX2) &  (!p$RING1B)) %>%
      lapply(function(x) {
        u[x] %>%
          lapply(anchors(ii), overlapsAny, .) %>%
          bind_cols() %>%
          rowSums() %>%
          `==`(2) %>%
          dd$d[.] %>%
          tibble(v = .)
      }) %>%
      bind_rows(.id = 'grp') %>%
      filter(grepl('only|three', grp)) %>%
      mutate(grp = gsub('and', '&', grp) %>%
               gsub('_', ' ', .))
    
    pd <- group_by(vv, grp) %>%
      do(getWhisks(.$v)) %>%
      ungroup() %>%
      mutate(nm = i$nm)
    
    
    ann <- compare_means(v ~ grp, vv)
    
    r <- diff(range(c(pd$ymin, pd$ymax)))
    
    annn <- ann %>%
      filter(group1 == 'All three') %>%
      mutate(y.position = max(pd$upper) + r / 10 * (1:n()),
             nm = i$nm)
    
    yl <- sprintf('%s - %s', names(i$fs)[1], names(i$fs)[2])
    #if (proj == 'mESC') {
    yl <- paste0('\u0394Pairwise interaction\n', yl)
    #}
    ggplot(pd, aes(x = grp, y = middle, group = grp, fill = grp, color = grp)) +
      geom_hline(yintercept = 0, color = 'black') +
      geom_boxplot(aes(ymin = ymin, ymax = ymax, fill = grp,
                       lower = lower, upper = upper,
                       notchupper = notchupper,
                       notchlower = notchlower,
                       middle = middle, color = grp),
                   notch = T, stat = "identity", width = .6,
                   show.legend = F) +
      geom_crossbar(aes(ymin = middle, ymax = middle),
                    color = "white", width = 0.4, fatten = 0) +
      stat_pvalue_manual(annn, inherit.aes = F, label = 'p.signif') +
      facet_wrap(~nm)+
      scale_y_continuous(position = 'right') +
      ylab(yl) +
      scale_fill_manual(values = clrs[c(2,3,1)]) +
      scale_color_manual(values = clrs[c(2,3,1)]) +
      theme(plot.background = element_blank(),
            legend.position = 'none',
            axis.text = element_text(color = 'black', size = 11),
            panel.background = element_blank(),
            strip.background = element_rect(fill = 'black'),
            panel.grid = element_blank(),
            axis.ticks = element_blank(),
            panel.grid.major = element_line(color = 'grey80', linetype = 'dashed'),
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text = element_text(color = 'white', size = 11))
  }, ., names(.)) -> pp
ggsave('f2_a.csight.pdf',pp$BT245,height = 5.2, width = 2.5, bg = 'transparent', device = cairo_pdf)

pdf('f2_a.euler.pdf', height = 2, width = 2)
u <- GRangesList(pks$BT245) %>% unlist() %>% reduce() %>% subsetByOverlaps(cgis$hg38)
p <- lapply(pks$BT245, function(x) overlapsAny(u, x)) %>% bind_cols()
euler(p, shape = 'ellipse') -> fit
names(fit$fitted.values) %>% tibble(nm = .) %>% 
  mutate(clr = case_when(nm == 'CBX2' ~ clrs[4], nm == 'RING1B' ~ clrs[1],
                         nm == 'H3K27me3' ~ clrs[3], 
                         nm == 'CBX2&RING1B&H3K27me3' ~ clrs[2],
                         T ~ '#dddddd')) %>% 
  pull(clr) %>%
  plot(fit, quantities = T, fills = .)
dev.off()



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

reg <- GRanges('chr1', IRanges(93.5e6, 95.5e6))
xbrks <- c(94e6, 95e6)
xlbs <- paste0(94:95, 'mb')

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

mats %>%
  bind_rows(.id = 'samp') %>%
  group_by(samp) %>%
  mutate(across(c(x, y), ~.x/sqrt(2))) %>%
  do(., rot(., 45)) %>%
  ungroup() %>%
  filter(between(y, 0, 50)) %>%
  filter(between(x, 101,301)) %>%
  ggplot(aes(x, y, fill = z)) +
  geom_raster() +
  coord_cartesian(expand = 0) +
  scale_fill_gradientn('Interaction',
                       colours = fall, na.value=fall[1],
                       trans = 'log', limits = c(.05, 1), 
                       oob = scales::squish, breaks = c(0.1,.5,1),
                       guide = guide_colorbar(barwidth = .5, barheight = 3.8,
                                              title.hjust = .5,
                                              title.vjust = 0.5,
                                              title.position = 'left')) +
  facet_nested('Hi-C' + samp~'pHGG (K27M)')+
  geom_label(aes(label = samp, color = samp), x = -Inf, y = Inf, hjust = 0, vjust = 1,
             data = ~distinct(., samp), label.size = NA, fill = '#ffffff99', show.legend = F) +
  scale_color_manual(values = tclrs) +
  scale_y_continuous(breaks = c(50),
                     labels = c('1mb')) +
  scale_x_continuous(breaks = c(451,551,651,751),
                     labels = c('110', '120', '130', '140mb')) +
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
        strip.text = element_text(color = 'white', size = 11)) -> p1

gn <- 'hg38'
cnds = c('K27M','K27M-KO')
nm = ''

bws = list(H3K27me3 = c('BT245_K27M_H3K27me3.bw',
                        'BT245_KO_H3K27me3.bw'),
           RING1B = c('BT245_K27M_RING1B.bw',
                      'BT245_KO_RING1B.bw'),
           CBX2 = c('BT245_K27M_CBX2.bw',
                    'BT245_KO_CBX2.bw'),
           H2Aub = c('BT245_K27M_H2AK119ub.bw',
                     'BT245_KO_H2AK119ub.bw')) %>%
  lapply(function(m) {
    setNames(m, cnds) %>%
      lapply(function(f) {
        print(f)
        import.bw(file.path('../data/bws', f), selection = BigWigSelection(reg))
      })
  })


bins <- tile(reg, 500)[[1]]
seqlevels(bins) <- lapply(unlist(bws, recursive = F), function(x) seqlevels(x)) %>%
  Reduce(intersect, .) 
xs <- mid(bins)
scs <- lapply(bws, function(x) {
  lapply(x, function(y) {
    seqlevels(y) <- seqlevels(bins)
    mcolAsRleList(y, 'score') %>%
      binnedAverage(bins, ., 'score') %>%
      score() %>%
      tibble(score = .) %>%
      mutate(idx = xs)
  }) %>% bind_rows(.id = 'samp')
}) %>% bind_rows(.id = 'mark') %>%
  mutate(samp = factor(samp, cnds),
         mark = factor(mark, names(bws)),
         proj = nm)

cgi <- cgis[[gn]] %>%
  subsetByOverlaps(reg) %>%
  as_tibble() %>%
  dplyr::select(1:3) %>%
  {lapply(bws, function(x) lapply(x, function(y) .) %>% bind_rows(.id = 'samp'))} %>%
  bind_rows(.id = 'mark') %>%
  mutate(samp = factor(samp, cnds),
         mark = factor(mark, names(bws)),
         proj = nm)


yl <- scs %>% group_by(mark, samp) %>% summarise(v = max(score)) %>%
  arrange(mark, samp) %>% pull(v)



scs %>%
  arrange(mark, samp) %>%
  mutate(grp = paste(mark, samp) %>% fct_inorder()) %>%
  split(., .$grp) %>%
  lapply(function(x) {
    m <- max(x$score)
    l <- if(m >= 10) {
      floor(m)
    } else if (m >= 1) {
      floor(m * 10) / 10
    } else if (m >= .1) {
      floor(m * 100) / 100
    } else if (m >= .01) {
      floor(m * 1000) / 1000
    } else if (m >= .01) {
      floor(m * 10000) / 10000
    } else {
      floor(m * 100000) / 100000
    }
    scale_y_continuous(limits = c(0, m),
                       breaks = c(0, l),
                       labels = c('', l),
                       expand = expansion(c(0, .2)))
  }) %>% unname() -> yls

p2 <- ggplot(scs, aes(x = idx, y = score)) +
  geom_line(aes(color = samp), size = .5) +
  geom_area(aes(fill = samp), alpha = .75) +
  geom_label(aes(label = samp, color = samp), x = -Inf, y = Inf, hjust = 0, vjust = 1,
             data = distinct(scs, samp, mark), label.size = NA, fill = '#ffffff99') +
  scale_x_continuous(seqnames(reg), breaks = xbrks,
                     labels = xlbs,
                     limits = c(start(reg), end(reg)),
                     expand = expansion(0)) +
  facet_nested(mark + samp ~ proj, scales = 'free')  +
  scale_color_manual(values = tclrs) +
  scale_fill_manual(values = tclrs) +
  facetted_pos_scales(y = yls) +
  ylab('CPM') +
  coord_cartesian(clip = 'off') +
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = 11, color = 'black'),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_line(color = 'black'),
        strip.text.x = element_blank(),
        strip.background.x = element_blank(),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white', size = 11))



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
  ggplot(aes(xmin = start, xmax = end, y = strand, color = strand,
             fill = strand, forward = forward)) +
  #geom_segment(aes(x = start, xend = end, y = strand, yend = strand),
  #             data = ~subset(., type == 'mRNA')) +
  #geom_gene_arrow(data = ~subset(., type == 'mRNA')) + 
  # geom_feature(aes(x = position, y = strand, forward = forward, color = strand),
  #              data = ~subset(., type == 'mRNA')) +
  # geom_gene_arrow(color = NA, data = ~subset(., type != 'mRNA')) +
  geom_gene_arrow(color = NA, data = ~subset(., type == 'mRNA')) +
  scale_y_discrete(breaks = c('-', '+'), labels = c('-', '+')) +
  coord_cartesian(xlim = c(start(reg), end(reg))) +
  scale_x_continuous(seqnames(reg),
                     breaks = xbrks, labels = xlbs,
                     expand = expansion(0)) +
  scale_fill_manual(values = gclrs) +
  scale_color_manual(values = gclrs) +
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

{wrap_plots(p1,p2,p3,nrow=3, heights = c(1,2.5,.3)) &
    theme(plot.background = element_blank()) } %>%
  ggsave('f2_b.pdf', ., width = 5.4, height = 7.5)


i <- fread('../data/aggr/lfc/BT245_K27M_H3K27me3.mat.gz', skip = 1, select = 1:6) %>%
  {which(!grepl('ENS', .$V4))}

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
  }) %>% bind_rows(.id = 'mark') -> agg

readRDS('../data/chip/rx.rds') %>%
  separate(samp, c('line', 'cond', NA), '_') %>%
  mutate(rx = log2(rx)) %>%
  split(.$line) %>%
  lapply(function(x) setNames(x$rx, x$cond)) -> rx

p1 <- agg %>%
  pivot_wider(names_from = 'cond', values_from = 'v') %>%
  mutate(dif = case_when(mark == 'H3K27me3' ~ (K27M + rx$BT245['K27M']) - (KO + rx$BT245['KO']),
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

p2 <- agg %>%
  pivot_wider(names_from = 'cond', values_from = 'v') %>%
  mutate(dif = case_when(mark == 'H3K27me3' ~ (K27M + rx$BT245['K27M']) - (KO + rx$BT245['KO']),
                         T ~ K27M - KO)) %>%
  dplyr::select(-K27M, -KO) %>%
  pivot_wider(names_from = 'mark', values_from = 'dif') %>%
  ggplot(aes(x = H3K27me3, y = RING1B, color = H2Aub)) +
  geom_hline(yintercept = 0, color = 'grey50') +
  geom_vline(xintercept = 0, color = 'grey50') +
  geom_scattermore(alpha = .5, pointsize = 4) +
  #annotate('label', x = Inf, y = Inf, hjust = 1.05, vjust = 1.2,
  #         label = sprintf('Pearson\'s r = %.2f', cor(dd[[m1]] - dd[[m4]], dd[[m2]] - dd[[m5]]))) +
  scale_color_gradientn(expression(Delta*'H2Aub'), colors = rev(brewer.rdylbu(20)), 
                        guide = guide_colorbar(barheight = 3),
                        limits = c(-1,1), oob = scales::squish,
                        breaks = -1:1) +
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

{wrap_plots(p1 + theme(axis.title.x = element_blank()), p2, ncol = 1) } %>%
  ggsave('f2_c.pdf', ., height = 4.93, width = 2.83, bg = 'transparent')


mclrs <- setNames(tableau20()[seq(5,11,2)],c('H2Aub', 'CBX2', 'RING1B', 'H3K27me3'))
names(mclrs) <- paste0('\u0394', names(mclrs))

gcor <- agg %>%
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

ggsave('f2_d.pdf', plt2, height = 2.3, width = 2.35, bg = 'transparent', dev = cairo_pdf)


