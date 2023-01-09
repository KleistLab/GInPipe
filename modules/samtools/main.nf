
process SAMTOOLS_INDEX {
  tag "$prefix"

  input:
  path bam

  output:
  path '*.bam.bai', emit: samtools_idx

  script:
  """
  samtools index $bam
  """
}

process SAMTOOLS_STATS {
  tag "$prefix"

  input:
  path bam
  path index

  output:
  path "*_stats.txt", emit: samtools_stats

  script:
  prefix = bam[0].getBaseName()
  """
  samtools idxstats $bam > "${prefix}_stats.txt"
  """
}

process SAMTOOLS_FAIDX {
  tag "$prefix"

  input:
  path ref

  output:
  path "*.fai", emit: faidx

  script:
  """
  samtools faidx $ref
  """
}

process SAMTOOLS_DICT {
  tag "$prefix"

  input:
  path ref

  output:
  path "*.dict", emit: dict

  script:
  """
  samtools dict $ref -o "${prefix}.dict"
  """
}

process SAMTOOLS_FLAGSTAT {
  tag "$prefix"

  input:
  path bam

  output:
  path "*_flagstat.txt", emit: flagstat

  script:
  """
  samtools flagstat $bam > "${prefix}_flagstat.txt"
  """
}

process SAMTOOLS_SORT {
  tag "$prefix"

  input:
  path bam

  output:
  path "*.bam", emit: bam_sort

  script:
  """
  samtools sort $bam > "${prefix}_sort.bam"
  """
}
