import pandas as pd
import re
from itertools import compress,chain

#' Filter/mask SNVs.
#'
#' @param snvs_per_seq A pandas.Series with all raw genotypes from the input represented in one string of single nucleotide variants (SNVs) seperated by a blank.
#' @param freq_threshold Threhsold for filtering (default = 2). Positions which don't occur more often than the threshold in the given sequence set are filtered out.
#TODO @param filter_mutation (optional) Boolean indicating if the threshold should be applied to each particular SNV instead of the position. (default = FALSE)
#' @param masked_positions (optional) Vector containing positions which should be masked/filtered out for the genotyping.
#' @param remove_indels (default = True) Boolean indicating if insertions and deletions should be disregarded.
#' @return A vector with the same genotypes from the input, but filterred represented as string of single nucleotide variants (SNVs) seperated by a blank

def filter_snvs(snvs_per_seq, freq_threshold=2, masked_positions=[], remove_indels=True, remove_all_insertions=True):
  #snvs_per_seq = seq_info_short_table['snvs']
  #tic("split snv strings")

  # splitting at blanks
  print("Split string into distinct SNVs")
  snvs_lists = snvs_per_seq.apply(str.split)

  #toc()

  #TODO remove insertions and deletions
  #always remove preceding insertions (pos = 0), optionally remove all insertions
  print("Remove insertions ")
  snvs_lists = snvs_lists.apply(filter_insertions, all=remove_all_insertions)
  if remove_indels:
    print("... and deletions")
    snvs_lists = snvs_lists.apply(filter_indels)

  # extract SNV positions
  print("Extract SNV positions")
  snv_pos_list = snvs_lists.apply(lambda snvs: list(map(extract_position, snvs)))
  
  # merge all lists of positions and count positions
  print("Count positions")
  pos_counts = snv_pos_list.explode().value_counts(sort=False)


  #filter positions which don't exceed the given threshold or which are masked
  # set is faster than list
  print("Filter positions which don't exceed the given threshold or which are masked")
  print("...Positions before filtering " + str(len(pos_counts)))

  pos_filtered = (set(list(pos_counts[pos_counts <= freq_threshold].index) + masked_positions))
  #tic("which pos must be deleted")
  ## 26 sec aufm Mac 
  #pos_counts_filtered = snv_pos_list.apply(lambda snvs: [x for x in snvs if x not in pos_filtered])
  pos_counts_filtered = snv_pos_list.apply(lambda snvs: [x not in pos_filtered for x in snvs])

  print("...Filtering out " + str(len(pos_filtered)) + " positions")

  # filter the actual snv table: check per seuqence if the snv is filtered or not with function compress
  print("Filter and recreate the actual SNV table")
  df = pd.DataFrame({'s': snvs_lists, 'p':pos_counts_filtered})
  snvs_filtered = df.apply(lambda x: ' '.join(list(compress(x.s, x.p))), axis=1)
  
  #snvs_filtered = snvs_lists.reset_index().apply(lambda x: (list(compress(snvs_lists[x.index], pos_counts_filtered[x.index])))).snvs
  # create string again as haplotype
  #snvs_filtered = snvs_filtered.apply(lambda x: ' '.join(x))

  return snvs_filtered 


def extract_position(snv):
  # substitutions and insertions
  pos = re.findall("(?<=[A-Za-z])\\d+(?=[A-Za-z])", snv)
  #deletions
  pos = re.findall("(?<=\\:)\\d+(?=\\:)", snv) if not pos else pos
  #pos =  re.findall("\\d+", snv)
  if pos:
    pos = int(pos[0]) 
  else:
     print("Weird pattern" + snv)
    # ValueError("Error while parsing SNVs: Unknown pattern in " + snv)  
  return pos

def filter_indels(snv_list):
  p = re.compile("^[A-Za-z]\\d+[A-Za-z]$")
  return list(filter(p.match, snv_list))

def filter_insertions(snv_list, all=False):
  if all:
    #all insertions
    p = re.compile("[A-Za-z]?\\d+[A-Za-z][A-Za-z]+$")
    snv_list = list(filter(lambda x: not p.match(x), snv_list))
  #insertions before sequence (also insertion of 1 nucleotide)
  p = re.compile("^0[A-Za-z]+$")
  return list(filter(lambda x: not p.match(x), snv_list))

### SNV mask parsing
def get_positions_from_string(pos_str):
  pos_intervals=pos_str.split(",")
  return list(chain.from_iterable(map(interval_to_range,pos_intervals)))

def interval_to_range(interval_str):
  try:
    i_s = interval_str.split("-")
    if len(i_s) == 2:
      # allow also reverse intervals 
      i_i = [int(i_s[0]), int(i_s[1])]
      return list(range(min(i_i), max(i_i)+1))
    elif len(i_s) == 1:
      # single positions can also be given
      return [int(i_s[0])]
  except ValueError as e:
    print(e)
  # the positions cannot be parsed 
  raise ValueError('The masked positions "' + interval_str + '" could not be parsed.\nPlease provide intervals in the form "pos1-pos2" and separate positions or intervals with a ",".')