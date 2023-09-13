


import bam.bam_to_fingerprints as bam
import utils.io_routines as io

bam_file = snakemake.input.bam
ref_file = snakemake.input.ref
result_file =  snakemake.output.snv

bam_converter =  bam.SAMtoFP(bam_file,ref_file)

snv_table = bam_converter.get_snv_table()


io.write_table(snv_table, result_file)

