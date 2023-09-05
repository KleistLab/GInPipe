import pandas as pd

def get_seq_info_per_day(seq_info_short_table, smoothing_window=7):
  """Infer sequence statistics per day, including number of sequences, haplotypes and mutations.
  :param seq_info_short_table: The parameter to be minimised
  :type seq_info_short_table: pandas.DataFrame
  :TODO param smoothing_window: Bandwidth for smoothing routine, defaults to 7
  :TODO type smoothing_window: int, optional
  :returns: Table with sequence information per day
  :rtype: pandas.DataFrame
  """
  seq_info_perDay_table = seq_info_short_table.groupby(['t', 'date'], as_index=False).agg( 
                          # count sequences
                          sequences = pd.NamedAgg(column ='index', aggfunc='count'),
                          # count number of haplotypes (unique sequences)
                          n_haplos = pd.NamedAgg(column ='snvs', aggfunc=lambda s: s.nunique()),
                          # count number of unique SNVs in all sequences
                          n_muts = pd.NamedAgg(column ='snvs', aggfunc=(lambda s: s.apply(str.split).explode().nunique())))

  #TODO add current/new haplos and mutation? 

  # count haplotypes and SNVs per day
  #seq_info_perDay.table$new_haplos <- 0
  #seq_info_perDay.table$muts <- 0
  #seq_info_perDay.table$new_muts <- 0
  
  #min_t <- min(seq_info_short.table$t, na.rm = T)
  #curr_haplos <- unique(seq_info_short.table$snvs[seq_info_short.table$t == min_t])
  #curr_muts <- unique(unlist(str_split(seq_info_short.table$snvs[seq_info_short.table$t == min_t], " ")))
  #seq_info_perDay.table$muts[seq_info_perDay.table$t == min_t] <- length(curr_muts)
  
  
  # for(t in seq(min(seq_info_short.table$t)+1, max(seq_info_short.table$t))) {
  #   idx <- which(seq_info_short.table$t == t)
  #   if(length(idx) > 0) {
  #     seq_info_perDay.table$new_haplos[seq_info_perDay.table$t==t] <- sum(!unique(seq_info_short.table$snvs[idx]) %in% curr_haplos)
  #     curr_haplos <- unique(c(curr_haplos, seq_info_short.table$snvs[idx]))
      
  #     muts <- unique(unlist(str_split((seq_info_short.table$snvs[idx]), " ")))
  #     new_muts <- muts[!muts %in% curr_muts]
  #     seq_info_perDay.table$muts[seq_info_perDay.table$t == t] <- length(muts)
  #     seq_info_perDay.table$new_muts[seq_info_perDay.table$t == t] <- length(new_muts)
  #     curr_muts <- unique(c(curr_muts, muts))
  #   }
  # }
  
  #TODO smooth values
  return seq_info_perDay_table

