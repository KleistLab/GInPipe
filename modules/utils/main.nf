
process REMOVE_WS {
  tag "$prefix"

  input:
  path fasta

  output:
  path "fasta", emit: fasta_ws

  script:
  """
  reformat.sh in=$fasta out="${prefix}_ws.fasta" underscore ignorejunk overwrite=true
  """
}

process REMOVE_DASHES {
  tag "$prefix"

  input:
  path fasta

  output:
  path "fasta", emit: fasta_rd

  script:
  """
  seqkit replace -s -p '-' -r 'N' $fasta -o "${prefix}_rd.fasta"
  """
}

process FASTP {
  tag "$prefix"

  input:
  path fastq1
  path fastq2

  output:
  tuple file("${fastq1.baseName}.trimmed.fq.gz"),
        file("${fastq2.baseName}.trimmed.fq.gz"), emit: trim_fastq
  
  script:
  """
  fastp --in1 $fastq1 --in2 $fastq2 --out1 ${fastq1.baseName}.trimmed.fq.gz \
        --out2 ${fastq2.baseName}.trimmed.fq.gz
  """
}
