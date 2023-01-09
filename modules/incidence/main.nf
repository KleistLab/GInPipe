
process INCIDENCE_EST {
  input:
  path bins
  path ref
  val cutoff
  val min_bin_size
  val min_days
  val max_days
  val group
  path vcf

  output:
  path "*.tsv", emit: inc

  script:
  """
  run_fp.py $bins $ref $cutoff $min_bin_size $min_days $max_days $group $vcf
  """
}

process INTERPOLATE {
  input:
  path inc
  path cases
  val group

  output:
  path "*.csv", emit: traj

  script:
  """

  """
}
