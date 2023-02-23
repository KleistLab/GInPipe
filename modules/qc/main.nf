process FASTQC{
    tag "$prefix"

    input:
    file fastq

    output:
    file('fastqc_report.html'), emit: fastqc_out

    """
    fastqc $fastq
    """
}

process MULTIQC{
    tag "$prefix"
    publishDir "${prefix}/qc", mode: 'copy', overwrite: false

    input:
    file fastqc

    output:
    file('multiqc_report.html')

    """
    multiqc .
    """
}

process TRIM_LOW_QUAL {
    tag "$prefix"
    publishDir "${prefix}/qc", mode: 'copy', overwrite: false

    input:
    path fastqc1
    path fastqc2

    output:
    path "*1_trimmed.fastqc", emit: trim_fastqc1
    path "*2_trimmed.fastqc", emit: trim_fastqc2

    """
    cutadapt -q 20,20 -o "${fastqc1.baseName}_trimmed.fastqc" -p "${fastqc2.baseName}_trimmed.fastqc" $fastqc1 $fastqc2
    """
}