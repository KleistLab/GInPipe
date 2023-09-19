import pandas as pd

def get_seq_info_per_day(seq_info_short_table):
  """Infer sequence statistics per day, including number of sequences, haplotypes and mutations.
  :param seq_info_short_table: The parameter to be minimised
  :type seq_info_short_table: pandas.DataFrame
  :returns: Table with sequence information per day
  :rtype: pandas.DataFrame
  """
  seq_info_perDay_table = seq_info_short_table.groupby(['t', 'date'], as_index=False).agg( 
                          # count sequences
                          sequences = pd.NamedAgg(column ='index', aggfunc='count'),
                          # count mutant sequences
                          #TODO
                          #mut_sequences = pd.NamedAgg(column ='index', aggfunc=lambda s: s[s not ""].),
                          # count number of haplotypes (unique sequences)
                          n_haplos = pd.NamedAgg(column ='snvs', aggfunc=lambda s: s.nunique()),
                          # count number of unique SNVs in all sequences
                          n_mut_types = pd.NamedAgg(column ='snvs', aggfunc=(lambda s: s.apply(str.split).explode().nunique())))

  #TODO add current/new haplos and mutation? 

  return seq_info_perDay_table

