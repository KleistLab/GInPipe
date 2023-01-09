
process SAMTOOLS_INDEX {
  tag "$prefix"

  input:
  path bam

  output:
  path '*.bam.bai', emit: samtools_idx

  script:
  """
  /Users/mariatrofimova/Documents/samtools-1.14/bin/samtools index $bam
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
  /Users/mariatrofimova/Documents/samtools-1.14/bin/samtools idxstats $bam > "${prefix}_stats.txt"
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
  /Users/mariatrofimova/Documents/samtools-1.14/bin/samtools faidx $ref
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
  /Users/mariatrofimova/Documents/samtools-1.14/bin/samtools dict $ref -o "${prefix}.dict"
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
  /Users/mariatrofimova/Documents/samtools-1.14/bin/samtools flagstat $bam > "${prefix}_flagstat.txt"
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
  /Users/mariatrofimova/Documents/samtools-1.14/bin/samtools sort $bam > "${prefix}_sort.bam"
  """
}
