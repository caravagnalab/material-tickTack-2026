is_fittable = function(fit, 
                       min_segment_length,
                       min_purity,
                       min_mutations_number_per_segment,
                       possible_k = c("2:1", "2:2", "2:0"), 
                       alpha = .05, 
                       min_mutations_number = 2) {
  get_clonal_peaks = function(k, purity) {
    multiplicities <- strsplit(k, ":") %>% unlist() %>% as.numeric()
    major <- multiplicities[1]
    n_tot <- sum(multiplicities)
    # get only Major and 1
    multiplicities <- c(1, major)
    peaks <- unique(multiplicities) * purity / (n_tot * purity + 2 * (1 - purity))
    return(sort(peaks))
  }
  
  segments = fit$cna
  mutations = fit$snvs
  purity = fit$purity
  
  segments <- segments %>%
    tidyr::drop_na(Major, minor)
  
  n_segments <- nrow(segments)
  
  accepted_segments_exist <- FALSE
  for (segment_idx in 1:n_segments) {
    print(segment_idx)
    # Segments
    segment <- segments[segment_idx, ]
    chr <- segment$chr
    
    segment_id <- paste(chr, segment$from, segment$to, sep = "_")
    
    # Get karyotype
    Major <- segment$Major
    minor <- segment$minor
    
    k <- paste(Major, minor, sep=':')
    print(k)
    
    peaks <- get_clonal_peaks(k, purity)
    
    
    if (!(k %in% possible_k )){
      print(paste0("seg ", segment_idx," has NOT and acceptable karyotype"))
    }
    if (k %in% possible_k & (segment$to-segment$from) <= min_segment_length){
      print(paste0("seg ", segment_idx," is NOT big enough, length = ", segment$to-segment$from))
    }
    if (k %in% possible_k & (segment$to-segment$from) > min_segment_length & purity >= min_purity & segment$CCF >=1) {
      print(paste0("seg ", segment_idx," is big enough, length = ", segment$to-segment$from))
      # Get info for mutations
      segment_mutations <- mutations %>%
        dplyr::filter(chr == segment$chr,.data$from > segment$from, .data$to < segment$to) %>%
        tidyr::drop_na(DP)
      
      
      accepted_mutations <- data.frame()
      if (nrow(segment_mutations) > 0) {
        # Check if mutation is inside CI
        probs <- c(alpha/2, 1 - alpha/2)
        
        DP <- segment_mutations$DP
        NV <- segment_mutations$NV
        
        accepted_idx <- lapply(1:length(DP), function(i) {
          #fisso i picchi per il segmento e vedo se tutte le mutazioni ricadono in almeno uno dei due intervalli intorno ai picchi
          # che sono diversi a seconda del valore di DP per la specifica mutazione
          for (p in peaks) {
            quantiles <- stats::qbinom(probs, DP[i], p)
            if ((NV[i] >= quantiles[1]) && (NV[i] <= quantiles[2])) {
              return(i)
            }
          }
        }) %>% unlist()
        
        # Get only good mutations
        accepted_mutations <- data.frame(DP = DP[accepted_idx], NV = NV[accepted_idx])
        
      }
      
      if (nrow(accepted_mutations) >= min_mutations_number_per_segment) {
        accepted_segments_exist <- TRUE
        print(paste0("and seg ", segment_idx," has ",nrow(accepted_mutations), " mutations"))
        # return(TRUE) UNCOMMENT to make it faster and comment print + last if else code
        
      } else( print(paste0("but seg ", segment_idx," has ",nrow(accepted_mutations), " mutations")))
    }
  }
  if (  accepted_segments_exist <- TRUE) {
    return(TRUE)
  }else{
    return(FALSE)
  }
  # return(FALSE) UNCOMMENT to make it faster and comment print + last if else code
}




plot_cnaqc_choose_K <- function(x, K, chromosomes = paste0('chr', c(1:22)), add_mobster=FALSE, max_Y_height = 6, cn = 'absolute', highlight = x$most_prevalent_karyotype, highlight_QC = FALSE) {
  
  cnaqc_x = CNAqc::init(mutations = x$mutations, cna = x$cna, purity = x$metadata$purity, ref = x$reference_genome)
  cnaqc_x$results_timing = x$results_timing
  
  plot_CNA = tickTack::plot_segments_tick_tack_CN(cnaqc_x, K = K) +
    ggplot2::theme(legend.position='right', panel.spacing = ggplot2::unit(0, "lines")) +
    ggplot2::labs(caption = NULL) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank())
  
  data_plot <- tickTack::plot_segments_tick_tack_data(x, K = K) +
    ggplot2::theme(legend.position='right',panel.spacing = ggplot2::unit(0, "lines")) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank())
  
  timing_plot <- tickTack::plot_segments_tick_tack(x, colour_by = "clock_mean", K = K) +
    ggplot2::theme(legend.position='right',panel.spacing = ggplot2::unit(0, "lines"))
  
  if(add_mobster){
    vaf_plot <- plot_vaf(x, K) +
      ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", color = "white", size = 20))
    
    # segment_plot <- plot_segments_h(x, chromosomes, max_Y_height, cn, highlight, highlight_QC) +
    #   ggplot2::theme(axis.title.x = element_blank())  # Keep chromosome labels only on this plot
    
    pA = timing_plot + CNAqc:::my_ggplot_theme() +
      ggplot2::theme(
        axis.ticks = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        legend.position = "left"
      )
    pB = plot_CNA + CNAqc:::my_ggplot_theme() +
      ggplot2::theme(
        axis.ticks = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        legend.position = "left"
      )
    pC = data_plot + CNAqc:::my_ggplot_theme() +
      ggplot2::theme(legend.position = "left") +
      ggplot2::labs(y = "VAF")
    pD = vaf_plot + CNAqc:::my_ggplot_theme()
    
    des_left = "
  AAAA#
  AAAAE
  AAAAE
  AAAAE
  AAAAE
  BBBBE
  BBBBE
  BBBBE
  CCCCE
  CCCCE
  DDDDD"
    
    pp = pA + pB + pC + patchwork::guide_area() + pD +
      patchwork::plot_layout(design = des_left, guides = "collect") &
      ggplot2::theme(legend.position = "bottom", legend.direction = "horizontal")
    pp
  } else {
    
    # segment_plot <- plot_segments_h(x, chromosomes, max_Y_height, cn, highlight, highlight_QC) +
    #   ggplot2::theme(axis.title.x = element_blank())  # Keep chromosome labels only on this plot
    
    pA = timing_plot + CNAqc:::my_ggplot_theme() +
      ggplot2::theme(
        axis.ticks = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        legend.position = "left"
      )
    pB = plot_CNA + CNAqc:::my_ggplot_theme() +
      ggplot2::theme(
        axis.ticks = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(),
        legend.position = "left"
      )
    pC = data_plot + CNAqc:::my_ggplot_theme() +
      ggplot2::theme(legend.position = "left") +
      ggplot2::labs(y = "VAF")
    pD = ggplot2::ggplot() + ggplot2::theme_void() + CNAqc:::my_ggplot_theme()
    
    des_left = "
  AAAA#
  AAAAE
  AAAAE
  AAAAE
  AAAAE
  BBBBE
  BBBBE
  BBBBE
  CCCCE
  CCCCE
  DDDDD"
    
    pp = pA + pB + pC + patchwork::guide_area() + pD +
      patchwork::plot_layout(design = des_left, guides = "collect") &
      ggplot2::theme(legend.position = "bottom", legend.direction = "horizontal")
    pp
  }
  
}



plot_segments_tick_tack_CN <- function(cnaqc_x, K = K, max_alleles = 6, chromosomes = paste0("chr", c(1:22))) {
  
  reference_genome <- get_reference(cnaqc_x$reference_genome)
  
  vfrom = reference_genome$from
  names(vfrom) = reference_genome$chr
  
  segments <- cnaqc_x$cna %>%
    dplyr::filter(.data$chr %in% chromosomes)
  
  absolute_segments <- segments %>%
    dplyr::mutate(from = .data$from + vfrom[.data$chr],
                  to = .data$to + vfrom[.data$chr])
  
  absolute_segments$karyotype = paste(absolute_segments$Major, absolute_segments$minor, sep = ":")
  absolute_segments = absolute_segments %>%
    dplyr::filter(.data$Major <= max_alleles & .data$minor <= max_alleles)
  
  k_colors = list(
    '1:1' = '#228B22CC',
    '1:0' = 'steelblue',
    '2:0' = 'turquoise4',
    '2:1' = ggplot2::alpha('orange', .8),
    '2:2' = 'firebrick3'
  )
  p = CNAqc:::blank_genome(cnaqc_x$reference_genome, chromosomes = paste0("chr", c(1:22)))
  
  # Add shadows
  p = p + ggplot2::geom_rect(
    data = absolute_segments %>% dplyr::filter(.data$karyotype %in% c('2:0', '1:0', '1:1', '2:1', '2:2')),
    ggplot2::aes(
      xmin = .data$from,
      xmax = .data$to,
      ymin = -Inf,
      ymax = Inf,
      fill = factor(.data$karyotype, levels = c('2:0', '1:0', '1:1', '2:1', '2:2'))
    ),
    alpha = .3
  ) +
    ggplot2::scale_fill_manual(values = k_colors) +
    ggplot2::guides(fill = ggplot2::guide_legend('', override.aes = list(alpha = 1)))
  
  # Add segments
  p = p +
    ggplot2::geom_segment(
      data = absolute_segments %>%
        dplyr::mutate(Major = as.numeric(.data$Major) + .1, minor = as.numeric(.data$minor) - .1) %>%
        dplyr::select(.data$karyotype, .data$chr, .data$from, .data$to, .data$Major, .data$minor) %>%
        tidyr::pivot_longer(!c(.data$karyotype, .data$chr, .data$from, .data$to)),
      ggplot2::aes(
        x = .data$from,
        xend = .data$to,
        y = .data$value,
        color = as.factor(.data$name)
      ),
      size=1.5
    ) +
    ggplot2::scale_color_manual(values = c("Major" = "red", "minor" ="steelblue")) +
    ggplot2::guides(color = ggplot2::guide_legend('')) +
    ggplot2::labs(y = "Allel count") +
    ggplot2::ylim(c(0, max_alleles))
  
  p
}



#' Plot mutation data from tickTack analysis
#'
#' @param x An object containing mutation data.
#' @param colour_by A character string specifying the variable to color points by (default: 'clock_mean').
#' @param K Numeric, number of clusters.
#' @return A ggplot object showing mutation data.
#' @export


plot_segments_tick_tack_data <- function(x, colour_by = "clock_mean", K = K) {
  
  mutations <- x$mutations
  results <- x$results_timing
  reference_genome <- get_reference(x$reference_genome)
  
  vfrom = reference_genome$from
  names(vfrom) = reference_genome$chr
  
  absolute_mutations <- mutations  %>%
    dplyr::mutate(from = .data$from + vfrom[.data$chr],
                  to = .data$to + vfrom[.data$chr])
  
  segments <- results$data$accepted_cna
  segments$from = lapply(segments$segment_name, function(s) {unlist(strsplit(s, "_"))[2]}) %>% unlist() %>% as.numeric()
  segments$to = lapply(segments$segment_name, function(s) {unlist(strsplit(s, "_"))[3]}) %>% unlist() %>% as.numeric()
  
  absolute_segments <- segments %>%
    dplyr::mutate(from = .data$from + vfrom[.data$chr],
                  to = .data$to + vfrom[.data$chr])
  
  summarized_results <- results$draws_and_summary[[K]]$summarized_results %>%
    dplyr::mutate(from = absolute_segments[.data$segment_id,]$from) %>%
    dplyr::mutate(to = absolute_segments[.data$segment_id,]$to) %>%
    dplyr::mutate(tau_mean = ifelse(.data$clock_mean  < 1, .data$clock_mean , 1)) %>%
    dplyr::mutate(tau_high = ifelse(.data$clock_high < 1, .data$clock_high, 1)) %>%
    dplyr::mutate(tau_low = .data$clock_low)
  
  accepted_mutations = data.frame()
  for (segment_idx in 1:nrow(summarized_results)) {
    segment <- summarized_results[segment_idx, ]
    # print(segment$chr)
    segment_mutations <- absolute_mutations %>%
      dplyr::filter(.data$chr == segment$chr, .data$from > segment$from, .data$to < segment$to) %>%
      tidyr::drop_na(.data$DP)
    segment_mutations <- segment_mutations %>% dplyr::mutate(karyotype = segment$karyotype)
    # print(nrow(segment_mutations))
    # if (nrow(segment_mutations)> 40){
    accepted_mutations <- dplyr::bind_rows(accepted_mutations, segment_mutations)
    # }
  }
  
  matched_mutations <- accepted_mutations
  
  k_colors = list(
    '2:0' = 'turquoise4',
    '2:1' = ggplot2::alpha('orange', .8),
    '2:2' = 'firebrick3'
  )
  
  CNAqc:::blank_genome(ref=x$reference_genome, chromosomes = paste0("chr", c(1:22)))+
    ggplot2::geom_point(
      data = matched_mutations,
      ggplot2::aes(
        x = .data$from,
        y = .data$NV / .data$DP,
        color = as.factor(.data$karyotype)
      ),
      alpha = 0.5,
      size=0.1
    ) +
    ggplot2::scale_color_manual(values = k_colors, name = "CN") +
    ggplot2::guides(
      color = ggplot2::guide_legend(override.aes = list(size = 4))) +
    ggplot2::labs(
      x = "Chromosome",
      y = bquote("Variant Allele Frequency (VAF)")
    )
}