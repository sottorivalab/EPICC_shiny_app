turn_off_clipping = function(x) {
  x = ggplot_gtable(ggplot_build(x))
  x$layout$clip = "off"
  
  wh = grep("strip-(t|r)", x$layout$name)
  for (i in wh) {
    x$grobs[[i]]$layout$clip = "off"
  }
  
  return(x)
}

pretty_violin = function(color="gray85", fill="gray10", ...) {
  list(
    ggplot2::geom_violin(fill=color, scale = "width", width=0.8),
    ggplot2::geom_boxplot(fill=fill, color=fill, outlier.shape = NA, width=0.2),
    ggplot2::stat_summary(fun = "mean", aes(shape="Mean"), color=color),
    ggplot2::stat_summary(fun = "median", aes(shape="Median"), color=color, size=2),
    ggplot2::labs(shape=""),
    ggplot2::scale_shape_manual(values=c(20, 95), breaks=c("Mean","Median")),
    ggplot2::guides(shape = ggplot2::guide_legend(override.aes = list(color="gray10", fill=NA, linetype=0))),
    cowplot::theme_cowplot()
  )
}


plot_atac_track = function(click_data, width=1000, n_win=50) {
  
  if (is.null(click_data)) {
    
    title=""
    
  } else {
    
    click_data$peak_alias = click_data$peak #with(recurrence_summary, id[peak == click_data$peak])
    click_data$peak = with(recurrence_summary, peak[id == click_data$peak])
    if (click_data$peak_alias %in% names(.longer_enh_label)) click_data$peak_alias = .longer_enh_label[ click_data$peak_alias]
    
    title =
      sprintf(
        "Case %s (%s)\n%s '%s' at %s",
        click_data$case,
        gsub("\n", " ", click_data$panel_x),
        gsub("s$", "", click_data$panel_y),
        click_data$peak,
        click_data$peak_alias
      )
  }
  
  if (is.null(click_data)) {
    plot = ggplot(NULL) + 
      geom_text(
        label="Please click a rectangle\nin the heatmap to display data", 
        aes(x=1, y=1)
      ) + 
      theme(axis.text = element_blank()) + 
      theme(axis.ticks = element_blank()) + 
      cowplot::theme_cowplot() + 
      ggtitle(title)
    
  } else {
    
    # merge data in groups
    .load_data = function(x) {
      d = readRDS(x)
      if (is.null(d)) d = GRanges()
      n = countOverlaps(wins, d)
      data.frame(pos=pos, n=n)
    }
    
    gr_peak = as(click_data$peak, "GRanges") + width
    wins = GenomicRanges::slidingWindows(gr_peak, width=n_win, step = ifelse(width>5000, 1, 5))[[1]]
    pos = (start(wins) + end(wins))/2
    
    wh = grepl(click_data$case, names(insertion_data)) | names(insertion_data) == "Normals"
    data_group = lapply(insertion_data[wh], .load_data)
    n_per_group = sapply(data_group, NROW)
    
    # calculate number of reads in windows:
    mt = match(names(data_group), group_annot$group)
    cor_factor = with(group_annot[mt,], (n_peak / size_peak - n_bkgr / sizebkgr))
    exp_bkgr_rate = with(group_annot[mt,], n_bkgr / sizebkgr)
    names(cor_factor) = names(exp_bkgr_rate) = names(data_group)
    
    # get window averaged tracks   
    data_per_win = 
      data_group %>% 
      reshape2::melt(measure.vars=c()) %>% 
      dplyr::mutate(L1_c = as.character(L1)) %>%
      dplyr::mutate(cpm = (n/n_win - exp_bkgr_rate[L1_c]) / cor_factor[L1_c])
    
    # get peak data
    if (!exists("peak_data")) peak_data = NULL
    if (any(names(cor_factor) %in% names(peak_data))) {
      
      max_cpm = max(data_per_win$cpm)
      
      peak_data_case = 
        peak_data[names(cor_factor)] %>% 
        (function(x) x[!sapply(x, is.null)]) %>%
        lapply(function(x) x[overlapsAny(x, gr_peak),]) %>% 
        lapply(as.data.frame) %>% 
        reshape2::melt(measure.vars=c()) %>% 
        dplyr::mutate(pos=start+end/2) %>% 
        dplyr::mutate(cpm=(1.05 + 0.025 * as.numeric(factor(L1))) * max_cpm)
      
    } else {
      peak_data_case = NULL
    }
    
    
    # plot peaks and tracks
    plot = ggplot(data_per_win, aes(x=pos, y=cpm, color=L1)) + 
      geom_line(size=1, alpha=0.8) + 
      xlab(paste0("Position on ", seqnames(gr_peak))) +
      ylab("Normalised number of insertions") + 
      scale_color_brewer(palette = "Set1") + 
      labs(color="Region") +
      theme(legend.position = "bottom") + 
      ggtitle(title)
    
    if (!is.null(peak_data_case)) {
      plot = plot  +
        geom_segment(
          data = peak_data_case, 
          aes(x=start, xend=end, y=cpm, yend=cpm), 
          size=2
        )
    }
    
  }
  
  return(plot)
}

plot_sample_stats = function(click_data) {
  
  if (is.null(click_data)) {
    
    plot = 
      ggplot(NULL) + 
      geom_text(
        label="Please click a rectangle\nin the heatmap to display data", 
        aes(x=1, y=1)
      ) + 
      theme(axis.text = element_blank()) + 
      theme(axis.ticks = element_blank()) + 
      cowplot::theme_cowplot()
    
  } else {
    
    # fix peak id
    click_data$peak = with(recurrence_summary, peak[id == click_data$peak])
    
    pl1 = insertion_data_annot %>% 
      dplyr::filter(patient == click_data$case) %>% 
      dplyr::filter(tissue_type != "normal") %>% 
      dplyr::filter(!is.na(purity)) %>% 
      ggplot(aes(
        x=region, 
        y=purity, 
        color=Hmisc::capitalize(as.character(tissue_type))
      )) +
      pretty_violin() + 
      xlab("Region") + 
      ylab("Purity") + 
      ylim(0, 1) + 
      guides(shape='none')+ 
      ggtitle("Sample purity") + 
      labs(color="Tissue") + 
      theme(legend.position="bottom") + 
      scale_color_brewer(palette="Set1")
    
    sids = grep(click_data$case, colnames(.cna_data), value = TRUE)
    
    pl2 = 
      .cna_data[click_data$peak, sids, drop=FALSE] %>% 
      reshape2::melt() %>% 
      cbind(annotation_from_barcode(.$Var2), .) %>% 
      dplyr::filter(tissue_type != "normal") %>% 
      dplyr::filter(!is.na(value)) %>%
      ggplot(aes(
          x=region, 
          y=scales::squish(value, c(0, 6)), 
          color=Hmisc::capitalize(as.character(tissue_type))
      )) +
      pretty_violin() + 
      xlab("Region") + 
      ylab("Copy-number") + 
      ylim(0, 6) + 
      guides(shape='none') + 
      ggtitle("Copy-number") + 
      labs(color="Tissue") + 
      guides(color="none") + 
      scale_color_brewer(palette="Set1")
    
    plot = 
      cowplot::plot_grid(
        pl1, pl2,
        nrow = 1,
        align = "h",
        axis = "tb"
      )
  }
  
  return(plot)
}
