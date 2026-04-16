
filter_tail_mutations <- function(df, purity, exclude = FALSE){
  variables_names <- colnames(df)
  
  if(exclude == TRUE){
    df <- df %>%
      rowwise() %>%
      mutate(
        peaks = list(get_clonal_peaks(k = karyotype, purity = purity)),
        
        p1 = dbinom(NV, DP, peaks[1]),
        p2 = dbinom(NV, DP, peaks[2]),
        
        mult_estimate = ifelse(p1 > p2, 1, 2),
        
        assigned_peak = peaks[mult_estimate],
        
        q_low  = stats::qbinom(0.025, DP, assigned_peak),
        q_high = stats::qbinom(0.975, DP, assigned_peak),
        
        is_tail = (NV < q_low) | (NV > q_high)
      ) %>%
      ungroup() %>%
      filter(is_tail == FALSE) %>% dplyr::select(c(variables_names,mult_estimate)) %>% 
      mutate( mult_estimate = ifelse(karyotype %in% c("2:2","2:1","2:0"), mult_estimate, NA))
    
  }else{
    
    df <- df %>%
      rowwise() %>%
      mutate(
        peaks = list(get_clonal_peaks(k = karyotype, purity = purity)),
        
        p1 = dbinom(NV, DP, peaks[1]),
        p2 = dbinom(NV, DP, peaks[2]),
        
        mult_estimate = ifelse(p1 > p2, 1, 2),
        
        assigned_peak = peaks[mult_estimate],
        
        q_low  = stats::qbinom(0.025, DP, assigned_peak),
        q_high = stats::qbinom(0.975, DP, assigned_peak),
        
        is_tail = (NV < q_low) | (NV > q_high)
      ) %>%
      ungroup() %>%
      dplyr::select(c(variables_names,mult_estimate)) %>% 
      mutate( mult_estimate = ifelse(karyotype %in% c("2:2","2:1","2:0"), mult_estimate, NA))
  }
  
  
  return(df)
}



get_clonal_peaks = function(k, purity) {
  multiplicities <- strsplit(k, ":") %>% unlist() %>% as.numeric()
  major <- multiplicities[1]
  n_tot <- sum(multiplicities)
  # get only Major and 1
  multiplicities <- c(1, major)
  peaks <- unique(multiplicities) * purity / (n_tot * purity + 2 * (1 - purity))
  return(sort(peaks))
}
