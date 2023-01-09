

process BINNING {
  tag "$prefix"

  input:
  path fasta
  path bam
  path index
  path stats

  output:
  path "bins", emit: list_bins

  script:
  prefix = bam[0].getBaseName()
  """
  mkdir bins
  run_binning_count.py $bam $stats $fasta

  """
}

process ERASE_EMPTY {
  tag "$prefix"

  input:
  path list_bins

  output:
  path "bins", emit: list_clean_bins

  script:
  """
  erase_empty_bins.py $list_bins
  """
}

process SORT_AND_BIN_WW {
  tag "$prefix"

  input:
  path list_clean_bins
  val region
  path ref
  val num_repeats

  output:
  path "sorted_binned", emit: sb

  script:
  """
  mkdir sorted_binned
  wastewater_run.py $list_clean_bins $ref $region $num_repeats
  """
}
