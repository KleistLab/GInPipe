#!/usr/bin/env nextflow

nextflow.enable.dsl=2
/*
 * Process import
 */
include { FASTQC; MULTIQC } from './modules/qc/main.nf'
include { MINIMAP_PAIRED; MINIMAP_INDEX; BWA_ALIGN_PAIRED; BWA_INDEX }  from './modules/align/main'
include { SAMTOOLS_INDEX; SAMTOOLS_STATS; SAMTOOLS_SORT_BY_NAME } from './modules/samtools/main'
include { BINNING_TEMP; BINNING_POS; ERASE_EMPTY; SORT_AND_BIN_WW }  from './modules/binning/main'
include { READS_TO_FP_PAIRED; ESTIMATOR } from './modules/reconstruction/main'

params.reads1 = '/Users/mariatrofimova/Documents/GitHub/SWAMPy/simulation_output/example_R1.fastq'
params.reads2 = '/Users/mariatrofimova/Documents/GitHub/SWAMPy/simulation_output/example_R2.fastq'
params.ref = '/Users/mariatrofimova/Documents/GitHub/SWAMPy/simulation_output/genomes_ref.fasta'
params.region = 'initialSequence'
//params.cons = '/Users/mariatrofimova/Documents/GitHub/GInSim-main/results_pmut_0.0001/fasta/sim_1_NS.fasta'
//params.true_table = '/Users/mariatrofimova/Documents/GitHub/GInSim-main/add_fitness_10_pos_exp_growth/table/sim_1_NS_shredded_0.005.tsv'
//params.masking_vcf = ''
// Number of reconstruction runs
//sed -n 's/\/1/|2020-01-01\/1/p' example_R1.fastq > example_R1_dates.fastq
params.num_repeats = 1
//params.mode = 'size','size'
//params.size = [5000,10000]
prefix = '1'

workflow {
  ref_ch = Channel.of(params.ref)
  reads1_ch = Channel.of(params.reads1)
  reads2_ch = Channel.of(params.reads2)
  //reads_ch = Channel.of(params.reads1, params.reads2)
  region_ch = Channel.value(params.region)
  // Only works with one mode at a time! 
  // TODO: Find a workaround
  modes = ['days']
  sizes = ['1']

  num_repeats_ch = Channel.value(params.num_repeats)
  //FASTQC(reads_ch)
  //MULTIQC(FASTQC.out.fastqc_out)
  MINIMAP_INDEX(ref_ch)
  MINIMAP_PAIRED(ref_ch, reads1_ch, reads2_ch, MINIMAP_INDEX.out.index)
  READS_TO_FP_PAIRED(prefix, MINIMAP_PAIRED.out.bam, ref_ch, region_ch)
  BINNING_TEMP(READS_TO_FP_PAIRED.out.fp_table)
  BINNING_POS(BINNING_TEMP.out.temp_bins, ref_ch)
  ESTIMATOR(BINNING_POS.out.pos_bins).view()
  //SAMTOOLS_INDEX(MINIMAP_PAIRED.out.bam)
  //SAMTOOLS_STATS(MINIMAP_PAIRED.out.bam, SAMTOOLS_INDEX.out.samtools_idx)
  //BINNING(ref_ch, MINIMAP.out.bam, SAMTOOLS_INDEX.out.idx ,SAMTOOLS_STATS.out.samtools_stats)
  //ERASE_EMPTY(BINNING.out.list_bins)
  //SORT_AND_BIN_WW(ERASE_EMPTY.out.list_clean_bins, region_ch, ref_ch, num_repeats_ch).view()
  //RECONSTRUCT_3(prefix, SORT_AND_BIN_WW.out.sb, ERASE_EMPTY.out.list_clean_bins, ref_ch, region_ch, freq_cutoff).view()
}


