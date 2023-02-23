process RECONSTRUCT_1 {
  publishDir(path: "$prefix", mode: "copy",)
  ignore_selector_warnings = true

  input:
  val prefix
  path sb
  path headers
  path ref
  val region
  val freq_cutoff

  output:
  path "reconstructed_bins/table_merged_phi_estimates_var_from_size.tsv", emit: recon

  script:
  """
  mkdir reconstructed_bins
  run_fp.py $sb $headers $ref $region $freq_cutoff
  """
}

process RECONSTRUCT_2 {
  publishDir(path: "$prefix", mode: "copy",)
  ignore_selector_warnings = true

  input:
  val prefix
  path sb
  path ref
  val region
  val freq_cutoff

  output:
  path "table_merged_phi_estimates_var_from_size.tsv", emit: recon

  script:
  """
  run_fp_2.py $sb $sb $ref $region $freq_cutoff
  """
}

process READS_TO_FP_PAIRED {
  tag "$prefix"
  ignore_selector_warnings = true

  input:
  val prefix
  path bam
  path ref
  val region

  output:
  path "*.tsv", emit: fp_table

  script:
  """
  run_fp_3.py $bam $ref $region
  """
}

process ESTIMATOR {
  tag "$prefix"
  publishDir(path: "$prefix", mode: "copy",)

  input:
  path pos_dir
  
  output:
  path "*.tsv", emit: estimates

  script:
  """
  run_incidence_estimator.py $pos_dir
  """
}
