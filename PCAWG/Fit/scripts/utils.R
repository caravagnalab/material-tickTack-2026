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