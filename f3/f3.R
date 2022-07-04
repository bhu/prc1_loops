library(data.table)
library(tidyverse)
library(rtracklayer)
library(scales)
library(pals)
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

list(NSD1 = list(gn = 'mm10', 
                 nm = 'This study\nmESC\n(NSD1-/-)',
                 reg = GRanges('chr1', IRanges(87134000,87214000)),
                 bws = list('WT' = 'mESC_Par_H3K27me3.bw',
                            'NSD1-KO' = 'mESC_NSD1KO_H3K27me3.bw')),
     NCRM1 = list(gn = 'hg38',
                  nm = 'This study\nhiPSC\n(NSD1+/-)', 
                  reg = GRanges('chr1', IRanges(158150000, 158230000)),
                  bws = list('WT' = 'hiPSC_WT_H3K27me3.bw',
                             'NSD1+/-' = 'hiPSC_NSD1het_H3K27me3.bw'))) %>%
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
      #geom_segment(aes(x = start, xend = end, y = strand, yend = strand),
      #             data = ~subset(., type == 'mRNA')) +
      #geom_gene_arrow(data = ~subset(., type == 'mRNA')) + 
      # geom_feature(aes(x = position, y = strand, forward = forward, color = strand),
      #              data = ~subset(., type == 'mRNA')) +
      # geom_gene_arrow(color = NA, data = ~subset(., type != 'mRNA')) +
      geom_gene_arrow(color = NA, data = ~subset(., type == 'mRNA')) +
      scale_y_discrete(breaks = c('-', '+'), labels = c('-', '+'), position = 'right') +
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
    
    yls <- split(scs, scs$samp) %>%
      lapply(function(s) {
        m <- max(s$score)
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
                           expand = expansion(c(0, .2)),
                           position = 'right')
      }) %>% unname()
    
    p <- ggplot(scs, aes(x = idx, y = score)) +
      geom_vline(aes(xintercept = start), data = cgi, linetype = 'dashed',
                 inherit.aes = F, color = bclr, size = 1) +
      geom_line(aes(color = samp), size = .5) +
      geom_area(aes(fill = samp), alpha = .75) +
      geom_label(aes(label = samp, color = samp), x = -Inf, y = Inf, hjust = 0, vjust = 1,
                 data = distinct(scs, samp, .keep_all = T), label.size = NA, fill = '#ffffff99') +
      scale_x_continuous(sprintf('%s:%s-%skb', seqnames(i$reg),
                                 formatC(as.integer(start(i$reg)/1e3), big.mark = ','), 
                                 formatC(as.integer(end(i$reg)/1e3), big.mark = ',')),
                         expand = expansion(0)) +
      facet_nested(t1 + t2 + t3 + samp ~ ., scales = 'free', switch = 'y',
                   strip = strip_nested(background_y = list(element_rect(fill = 'grey0'),
                                                            element_rect(fill = 'grey20'),
                                                            element_rect(fill = 'grey40'),
                                                            element_blank(),
                                                            element_blank()),
                                        text_y = list(element_text(),
                                                      element_text(),
                                                      element_text(),
                                                      element_blank(),
                                                      element_blank()),
                                        bleed = T))  +
      scale_color_manual(values = tclrs) +
      scale_fill_manual(values = tclrs) +
      facetted_pos_scales(y = yls) +
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
            strip.text = element_text(color = 'white', size = 11))
    
    
    
    
    p
    #if (proj != 'BT245') {
    p <- p + theme(axis.title.y = element_blank())
    #}
    #p
    cowplot::plot_grid(p + theme(axis.title.x = element_blank(),
                                 plot.margin = margin(1, 1, -1, 1)),
                       pg + theme(plot.margin = margin(-1, 1, 1, 1)),
                       nrow = 2, rel_heights = c(3, 1), align = 'v', axis = 'lr')
  }, ., names(.)) -> ps

ggsave('f3_a.m.pdf', ps$NSD1, width = 3, height = 2.1)
ggsave('f3_a.h.pdf', ps$NCRM1, width = 3, height = 2.1)



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
      scale_y_continuous('H3K27me3 confinement', position = 'right') +
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

ggsave('f3_b.m.pdf', ps$`mESC (NSD1)`, height = 2.26, width = 1.6, bg = 'transparent')
ggsave('f3_b.h.pdf', ps$hiPSC, height = 2.26, width = 1.6, bg = 'transparent')



aggr <- list(NSD1KO = c(WT = 'mESC_Par_H3K27me3.mat.gz',
                        NSD1KO = 'mESC_NSD1KO_H3K27me3.mat.gz'),
             NSD1het = c(WT = 'hiPSC_WT_H3K27me3.mat.gz',
                         'NSD1+/-' = 'hiPSC_NSD1het_H3K27me3.mat.gz')) %>%
  lapply(function(x) {
    i <- lapply(x, function(y) {
      p <- file.path('..', 'data', 'aggr', y)
      cgi <- fread(p, skip = 1, select = 1:6)$V4 %>%
        grepl('^EN', .) %>%
        {which(!.)}
      fread(p, skip = 1, drop = 1:6) %>%
        .[,996:1005] %>%
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

names(aggr) %>%
  setNames(., .) %>%
  lapply(function(p) {
    aggr[[p]] %>%
      mutate(proj = p,
             i = as.character(i),
           se = std / sqrt(num)) %>%
      ggplot(aes(x = idx, y = mu, color = i, fill = i, ymin = mu - se, ymax = mu + se)) +
      geom_ribbon(color = NA, alpha = .5) +
      geom_line() +
      facet_wrap(~proj, scales = 'free', nrow = 1) +
      scale_x_continuous(breaks = c(900, 1000, 1100),
                         labels = c('-100kb', 'CGI', '+100kb'))+
      scale_color_manual(values = tclrs) +
      scale_fill_manual(values = tclrs) +
      coord_cartesian(xlim = c(900, 1100)) +
      scale_y_continuous('H3K27me3 CPM', position = 'right') +
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
            axis.line.x = element_line(color = 'black'))
      
  }) -> ps

ggsave('f3_c.m.pdf', ps$NSD1KO, height = 2.26, width = 2.3, bg = 'transparent')
ggsave('f3_c.h.pdf', ps$NSD1het, height = 2.26, width = 2.3, bg = 'transparent')



ps <- readRDS('../data/pup/H3K27me3_CGI.rds') %>%
  separate(samp, c('proj', 'cond'), '_') %>%
  filter(proj == 'hiPSC' | cond %in% c('Par','NSD1KO')) %>%
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

ggsave('f3_d.m.pdf',ps$mESC, height = 2.42,width=2.8,bg='transparent')
ggsave('f3_d.h.pdf',ps$hiPSC, height = 2.42,width=2.8,bg='transparent')



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
  }, ., names(.)) -> ps

{wrap_plots(ps[c('NSD1KO','NSD1het')], nrow = 1) &
    theme(plot.background = element_blank())} %>%
  ggsave('f3_d_p.pdf', ., width = 4, height = 3.4)

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
ggsave('f3_e.csight.pdf',pp$mESC,height = 5.2, width = 2.5, bg = 'transparent', device = cairo_pdf)

pdf('f3_e.euler.pdf', height = 2, width = 2)
u <- GRangesList(pks$mESC) %>% unlist() %>% reduce() %>% subsetByOverlaps(cgis$mm10)
p <- lapply(pks$mESC, function(x) overlapsAny(u, x)) %>% bind_cols()
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

cgis <- list.files('../data/cgi', full.names = T) %>%
  setNames(., sub('.bed', '', basename(.))) %>%
  lapply(import)

reg <- GRanges('chr12', IRanges(56e6, 58e6))
xbrks <- seq(56e6, 58e6, 1e6)
xlbs <- paste0(xbrks/1e6, 'mb')

res <- 10000

mats <- c('WT' = 'mESC_WT', 'NSD1-KO' = 'mESC_NSD1KO') %>%
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
tclrs <- c('#e45756', '#4c78a8')


p1 <-  mats %>%
  bind_rows(.id = 'samp') %>%
  group_by(samp) %>%
  mutate(across(c(x, y), ~.x/sqrt(2))) %>%
  do(., rot(., 45)) %>%
  ungroup() %>%
  filter(between(y, 0, 60)) %>%
  filter(between(x, 101,301)) %>%
  mutate(samp = factor(samp, c('WT','NSD1-KO'))) %>%
  ggplot(aes(x, y, fill = z)) +
  geom_raster() +
  coord_cartesian(expand = 0) +
  scale_fill_gradientn('Interaction',
                       colours = fall, na.value=fall[1],
                       trans = 'log',limits = c(.05, 1),
                       oob = scales::squish, breaks = c(0.1,.5,1),
                       guide = guide_colorbar(barwidth = .5, barheight = 3.8,
                                              title.hjust = .5,
                                              title.vjust = 0.5,
                                              title.position = 'left')) +
  facet_nested('Hi-C' + samp~'mESC (WT)') +
  geom_label(aes(label = samp, color = samp), x = -Inf, y = Inf, hjust = 0, vjust = 1,
             data = ~distinct(., samp), label.size = NA, fill = '#ffffff99',
             show.legend = F) +
  scale_color_manual(values = tclrs) +
  scale_y_continuous(breaks = c(25, 50),
                     labels = c('500kb', '1mb')) +
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
        strip.text = element_text(color = 'white', size = 11))

gn <- 'mm10'
cnds = c('WT','NSD1-KO')
nm = ''

bws = list(H3K36me2 = c(WT = "mESC_Par_H3K36me2.bw", 
                        "NSD1-KO" = "mESC_NSD1KO_H3K36me2.bw"),
           H3K27me3 = c(WT = "mESC_Par_H3K27me3.bw",
                        "NSD1-KO" = "mESC_NSD1KO_H3K27me3.bw"),
           RING1B = c(WT = "mESC_Par_RING1B.bw",
                      "NSD1-KO" = "mESC_NSD1KO_RING1B.bw"),
           CBX2 = c(WT = "mESC_Par_CBX2.bw",
                    "NSD1-KO" = "mESC_NSD1KO_CBX2.bw")) %>%
  lapply(function(m) {
    m %>%
      lapply(function(f) {
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
  mutate(samp = factor(samp, names(bws[[1]])),
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


yl <- scs %>% group_by(mark) %>% summarise(v = max(score)) %>%
  arrange(mark) %>% deframe()


yls <-  scs %>%
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


p2 <- ggplot(scs, aes(x = idx, y = score )) +
  geom_line(aes(color = samp), size = .5) +
  geom_area(aes(fill = samp), alpha = .75) +
  geom_label(aes(label = samp, color = samp), x = -Inf, y = Inf, hjust = 0, vjust = 1,
             data = ~distinct(., samp, mark), label.size = NA, fill = '#ffffff99') +
  scale_x_continuous(expand = expansion(0)) +
  facet_nested(mark + samp ~ proj, scales = 'free')  +
  scale_color_manual(values = tclrs) +
  scale_fill_manual(values = tclrs) +
  ylab('CPM') +
  facetted_pos_scales(y = yls) +
  #facetted_pos_scales(y = list(
  #  mark == 'H3K27me3' ~ scale_y_continuous(expand = expansion(c(0, .2)),
  #                                          breaks = scales::pretty_breaks(2),
  #                                          limits = c(0, yl['H3K27me3'])),
  #  mark == 'H3K36me2' ~ scale_y_continuous(expand = expansion(c(0, .2)),
  #                                          breaks = scales::pretty_breaks(2),
  #                                          limits = c(0, yl['H3K36me2'])),
  #  mark == 'RING1B' ~ scale_y_continuous(expand = expansion(c(0, .2)),
  #                                          breaks = scales::pretty_breaks(2),
  #                                          limits = c(0, yl['RING1B'])),
  #  mark == 'CBX2' ~ scale_y_continuous(expand = expansion(c(0, .2)),
#                                          breaks = scales::pretty_breaks(2),
#                                          limits = c(0, yl['CBX2']))
#)) +
coord_cartesian(clip = 'off') +
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        axis.line.y = element_line(color = 'black'),
        strip.text.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background.x = element_blank(),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(color = 'white', size = 11))



gclrs <- c('#f28e2b', '#b07aa1')



subsetByOverlaps(txs$mm10, reg) %>%
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

{wrap_plots(p1,p2,p3,nrow=3, heights = c(1,2,.3)) &
    theme(plot.background = element_blank()) } %>%
  ggsave('f3_f.pdf', ., width = 5.4, height = 7)


fs <- list.files('../data/aggr', pattern = 'mESC_(NSD1KO|Par)', full.names = T) %>%
  tibble(f = .) %>%
  mutate(s = basename(f)) %>%
  separate(s, c(NA, 'cond', 'mrk'), '_', F) %>%
  mutate(mrk = sub('.mat.gz', '', mrk),
         mark = ifelse(grepl('input', mrk), 'inp', 'chip'),
         mrk = sub('inputFor', '', mrk)) 

odrs <- fs %>%
  split(., .$cond) %>%
  lapply(function(x) {
    split(x, x$mrk) %>%
      lapply(function(y) {
        split(y$f, y$mark) %>%
          lapply(fread, skip = 1, drop = 1:6) %>%
          {log2((.$chip + 1) / (.$inp + 1))} %>%
          .[,990:1010] %>%
          rowMeans(na.rm = T) -> v
        head(order(v, decreasing = T), 1000)
      })
  })

rs <- lapply(odrs, `[[`, 'H3K27me3')
rs$union <- unique(unlist(rs))
rs$intersect <- intersect(rs$Par, rs$NSD1KO)

aggr <- fs %>%
  split(., .$cond) %>%
  lapply(function(x) {
    split(x, x$mrk) %>%
      lapply(function(y) {
        m <- split(y$f, y$mark) %>%
          lapply(fread, skip = 1, drop = 1:6) %>%
          {log2((.$chip + 1) / (.$inp + 1))} %>%
          as.matrix()
        
        lapply(rs, function(z) {
          m[z,] %>%
            {tibble(mu = colMeans(., na.rm = T),
                    std = colSds(., na.rm = T),
                    num = colSums(is.finite(.)))} %>%
            mutate(idx = 1:n())
        }) %>% bind_rows(.id = 'reg')
      }) %>% bind_rows(.id = 'mark')
  }) %>% bind_rows(.id = 'cond')

tclrs <- c('#e45756', '#4c78a8')
aggr %>%
  mutate(se = std / sqrt(num)) %>%
  filter(reg == 'union') %>%
  mutate(cond = factor(cond, c('Par', 'NSD1KO')),
         mark = factor(mark, c('H3K27me3', 'RING1B','CBX2'))) %>%
  ggplot(aes(x = idx, y = mu, color = cond, fill = cond, ymin = mu - se, ymax = mu + se)) +
  geom_ribbon(alpha = .5, color = NA)+
  geom_line() +
  facet_grid(mark ~ 'PRC1 reading of H3K27me3', scales = 'free_y') +
  scale_fill_manual(values = tclrs) +
  scale_color_manual(values = tclrs) +
  coord_cartesian(xlim = c(750, 1250)) +
  scale_x_continuous(breaks = c(750, 1000, 1250),
                     labels = c('-250kb','H3K27me3-\nenriched CGI','+250kb')) +
  scale_y_continuous('ChIP-seq enrichment', breaks = scales::pretty_breaks(2)) +
  theme(axis.text = element_text(size = 11, color = 'black'),
        strip.background = element_rect(fill = 'black'),
        strip.text = element_text(size = 11, color = 'white'),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.position = c(1,1),
        legend.justification = c(1,1),
        axis.title.x = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        plot.background = element_blank()) -> p
ggsave('f3_g.pdf', p, height = 3.4, width = 3.2, bg = 'transparent')

load('../data/pks/mESC.CBX2.rda')
pks <- dba.report(db, th = 1)  

load('../data/chip/mESC.10kb.rda')
bl <- fread('../data/chip/mm10-blacklist.v2.bed.gz', select = 1:3, 
            col.names = c('chr', 'start', 'end')) %>%
  mutate(start = start + 1) %>%
  makeGRangesFromDataFrame()

bins <- mutate(bins, start = start + 1) %>%
  makeGRangesFromDataFrame()
k <- apply(cts, 1, max) & !overlapsAny(bins, bl)
bins <- bins[k]
names(cts) %>%
  grep('inputFor', ., value = T) %>%
  setNames(., sub('inputFor', '', .)) %>%
  lapply(function(x) {
    y <- sub('inputFor', '', x)
    l1 <- lsz[lsz$samp == y,]
    l2 <- lsz[lsz$samp == x,]
    f <- c(l1$endo, l2$endo) %>% 
      {min(.) / .}
    rx <- (l1$endo / l1$exo) / (l2$endo / l2$exo)
    if (grepl('CBX2|RING1B', x)) {
      a <- cts[[y]] * f[1] + 1
    } else {
      a <- cts[[y]] * rx * f[1] +1 
    }
    b <- cts[[x]] * f[2] + 1
    log2(a / b)[k]
  }) %>%
  bind_cols() -> d

hits <- findOverlaps(pks, bins) %>%
  as("List")

mcols(pks) <- lapply(d, function(x) {
  extractList(x, hits) %>%
    mean(na.rm = T)
}) %>% bind_cols() %>%
  cbind(mcols(pks), .)


pd <- as_tibble(pks) %>%
  mutate(grp = case_when(p.value < .05 & Fold < 0 ~ 'Lost',
                         p.value < .05 & Fold > 0 ~ 'Gain',
                         T ~ 'Other'))

clrs <- c('dodgerblue3', 'firebrick2', 'gray50')
p2 <- ggplot(pd, aes(x = Fold, y = -log10(p.value), color = grp)) +
  geom_vline(xintercept = 0) +
  geom_point(pch = 16, size = .5) +
  scale_color_manual(values = clrs,
                     guide = guide_legend(nrow = 3)) +
  facet_wrap(~'\u0394CBX2 binding') +
  labs(x = expression(log[2]~'NSD1-KO'/'WT'),
       y = expression(-log[10]~'p')) +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        legend.background = element_blank(),
        legend.position = c(1,1),
        legend.justification = c(1,1),
        panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.direction = 'horizontal',
        legend.title = element_blank(),
        axis.text = element_text(color = 'black', size = 11),
        panel.grid.major = element_line(color = 'grey80', linetype = 'dashed'),
        strip.text = element_text(color = 'white', size = 11),
        strip.background = element_rect(fill = 'black')) 

compare_means(mESC_Par_H3K27me3 ~ grp, pd) %>%
  mutate(y.position = 3.8 + (1:n()) * .6) -> ann

p3 <- ggplot(pd, aes(x = grp, y = mESC_Par_H3K27me3)) +
  geom_boxplot(aes(color = grp, fill = grp), notch = T, outlier.size = 0.5, outlier.shape = 16) +
  stat_summary(geom = "crossbar", width=0.2, fatten=0, color="white", 
               fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) }) +
  scale_color_manual(values = clrs) +
  scale_fill_manual(values = clrs) +
  stat_pvalue_manual(ann, label = 'p.signif', inherit.aes = F) +
  scale_y_continuous('H3K27me3 (WT)',
                     position = 'right', expand = expansion(c(.01,.1))) +
  theme(plot.background = element_blank(),
        panel.background = element_rect(fill = NA, color = 'black', size = 1),
        panel.grid = element_blank(),
        legend.key = element_blank(),
        legend.direction = 'horizontal',
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(color = 'black', size = 11),
        panel.grid.major = element_line(color = 'grey80', linetype = 'dashed'),
        strip.text = element_text(color = 'white', size = 11),
        strip.background = element_rect(fill = 'black'),
        legend.position = 'none') 

{wrap_plots(p2,p3,nrow=1, widths = c(1.5,1)) &
    theme(plot.background = element_blank()) } %>%
  ggsave('f3_h.pdf',.,height=3.5,width=3.6, device = cairo_pdf, bg = 'transparent')

