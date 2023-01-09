process RECONSTRUCT_1 {
  publishDir(path: "ww_$prefix", mode: "copy",)
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
  publishDir(path: "standard_$prefix", mode: "copy",)
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
  run_fp__.py $sb $sb $ref $region $freq_cutoff
  """
}