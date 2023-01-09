#!/usr/bin/env nextflow

nextflow.enable.dsl=2
/*
 * Process import
 */
include { MINIMAP_INDEX; MINIMAP }  from './modules/align/main'
include { SAMTOOLS_INDEX; SAMTOOLS_STATS } from './modules/samtools/main'
include { BINNING; ERASE_EMPTY; SORT_AND_BIN_WW }  from './modules/binning/main'
include { RECONSTRUCT_2 } from './modules/reconstruction/main'

//params.reads = '/Users/mariatrofimova/Documents/GitHub/GInSim-main/results_pmut_0.0001/fasta/sim_1_NS_shredded_0.05.fasta'
params.ref = '/Users/mariatrofimova/Documents/GitHub/GInSim-main/results_pmut_0.0001/reference/sim_3.fasta'
params.region = 'initialSequence'
params.cons = '/Users/mariatrofimova/Documents/GitHub/GInSim-main/results_pmut_0.0001/fasta/sim_3_NS.fasta'
//params.true_table = '/Users/mariatrofimova/Documents/GitHub/GInSim-main/add_fitness_10_pos_exp_growth/table/sim_1_NS_shredded_0.005.tsv'
//params.masking_vcf = ''
// Number of reconstruction runs
params.num_repeats = 1
freq_cutoff = 0
prefix = '1'

workflow {
  ref_ch = Channel.fromPath(params.ref)
  reads_ch = Channel.fromPath(params.cons)
  region_ch = Channel.value(params.region)

  num_repeats_ch = Channel.value(params.num_repeats)
  MINIMAP_INDEX(ref_ch)
  MINIMAP(ref_ch, reads_ch, MINIMAP_INDEX.out.index)
  SAMTOOLS_INDEX(MINIMAP.out.bam)
  SAMTOOLS_STATS(MINIMAP.out.bam, SAMTOOLS_INDEX.out.samtools_idx)
  BINNING(ref_ch, MINIMAP.out.bam, SAMTOOLS_INDEX.out.samtools_idx ,SAMTOOLS_STATS.out.samtools_stats)
  ERASE_EMPTY(BINNING.out.list_bins)
  RECONSTRUCT_2(prefix, ERASE_EMPTY.out.list_clean_bins, ref_ch, region_ch, freq_cutoff).view()
}

