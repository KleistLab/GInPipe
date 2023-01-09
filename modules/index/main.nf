
process SAMTOOLS_SORT {
  tag "$prefix"
  publishDir params.out, mode: 'copy'

  input:
  path bam

  output:
  path '*.bam', emit: bam

  script:
  prefix = bam[0].getBaseName()
  """
  /Users/mariatrofimova/Documents/samtools-1.14/bin/samtools sort ${bam} > ${prefix}_sorted.bam
  """
}

process SAMTOOLS_INDEX {
  tag "$prefix"
  publishDir params.out, mode: 'copy'

  input:
  path bam from bam_sorted_files

  output:
  path '*.bam.bai' into bam_index_files

  script:
  """
  /Users/mariatrofimova/Documents/samtools-1.14/bin/samtools index ${bam}
  """
}

process SAMTOOLS_STATS {
  tag "$prefix"
  publishDir params.out, mode: 'copy'

  input:
  path bam from bam_sorted_files

  output:
  path "*_stats.txt" into bam_stats_files

  script:
  prefix = bam[0].getBaseName()
  """
  /Users/mariatrofimova/Documents/samtools-1.14/bin/samtools idxstats ${bam} > ${prefix}_stats.txt
  """
}
