process PLOT {
  tag "$prefix"

  input:
  path recon
  path trueN

  output:
  path '*.pdf', emit: plot

  script:
  """
  Rscript plot_routines.R $recon $trueN
  """
}
