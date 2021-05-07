#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 18:426:09 2021

@authors: Maureen Smith, Maria Trofimova
"""

import sys
import argparse
import os
import random
import math

from sequence_evolution import sequenceEvol
from output_writer import writer



# initialize argument parser
parser = argparse.ArgumentParser(description='Simulate sequence evolution.')

# Required arguments
parser.add_argument('-file_prefix', required=True,
                    help='file_prefix for fasta file to write the simulated sequences to')

parser.add_argument('-output', required=True,
                    help='output path were the simulated fasta and table files are saved')

parser.add_argument('-L', type=int, required=True,
                    help='sequence length')

parser.add_argument('-N_init', type=int, required=True,
                    help='number of initial sequences')

parser.add_argument('-p_mut', type=float, required=True,
                    help='mutation rate')

parser.add_argument('-t_final', type=int, required=True,
                    help='number of time steps')

parser.add_argument('-p_death', type=float, nargs=1, required=False,
                    help='death rate')


# Optional arguments
parser.add_argument('-p_rep', type=float,
                    help='(optional) replication rate')

parser.add_argument('-p_rep2', type=float,
                    help='(optional) second replication rate')

parser.add_argument('-subsampling', type=float,
                    help='(optional) subsampling of all samples. if zero, take random subsample')

# Switch
parser.add_argument('--switch_orig', action='store_true', required=False,
					help='(optional) In the case of introductions and two replication rates: '
                         'If true: switch to secon replication rate at the same time point as the original outbreak. '
                         'If false: switch is at the half of the respective time frame.')

parser.add_argument('-init_seq', nargs='+',
                    help='(optional) Initial sequence. If not given, a random sequence is created.')

parser.add_argument('-sub_rel', type=float, nargs='+', required=False,
                    help='(optional) relative subsampling of the complete sequence set')

parser.add_argument('-sub_abs', type=int, nargs='+', required=False,
                    help='(optional) absolute amount of subsampled sequences per day')


args = parser.parse_args()

print("*"*100)
print("Running simulation of sequence evolution with argruments\n")
for a in vars(args):
    print(a, getattr(args, a))
    #print(' {} {}'.format(a, getattr(args, a) or ''))
print("*"*100)

evol_run = sequenceEvol(length=args.L,
                            p_repl=args.p_rep,
                            p_repl2=args.p_rep2,
                            p_death=args.p_death,
                            p_mut=args.p_mut,
                            N_init=args.N_init,
                            t_start=0,
                            t_final=args.t_final,
                            t_switch=math.floor(args.t_final / 2)
                        )

# run simulation
time_trajectory = evol_run.evolve_poi()

wr = writer(outputpath=args.output, file_prefix=args.file_prefix)

# write initial sequence as reference
wr.write_reference(evol_run.init_seq)

# header=hCoV-19/Italy/LAZ-INMI1-isl/2020|EPI_ISL_410545|2020-01-29
header_prefix=">NS|"
file_suffix = "NS"


# write fasta with all sequences
df_NS = wr.write_fasta(file_suffix=file_suffix, header_prefix=header_prefix, species_dict=time_trajectory)
df_NS["true_N"] = df_NS["sampled_N"]
wr.write_table(table=df_NS, file_suffix=file_suffix)
wr.write_config_yaml(file_suffix=file_suffix)

# write fasta with all absolute subsample
if args.sub_abs is not None:
    ts = range(args.t_final)
    for s_abs in args.sub_abs:
        print("---  Subsample sequence set taking " + str(s_abs) + " ---")
        header_prefix = header_prefix=">WS|" + str(s_abs) + "|"
        file_suffix = "WSABS_"+ str(s_abs)

        subsampled_time_trajectory = []
        for t in ts:
            num_seq = sum(time_trajectory[t].values())
            time_trajectory_sub = {}
            # take abolsute subsample or N(t) if less
            if(num_seq > s_abs):
                # sampling the particular sequences for each time point without replacemenet
                seq_subset = random.sample(list(time_trajectory[t].keys()), counts=time_trajectory[t].values(), k=s_abs)
                # counts sampled sequences
                for s in seq_subset:
                    time_trajectory_sub[s] = time_trajectory_sub.get(s, 0) + 1
            else:
                time_trajectory_sub = time_trajectory[t]
            subsampled_time_trajectory.append(time_trajectory_sub)

        df_WS_abs = wr.write_fasta(file_suffix=file_suffix, header_prefix=header_prefix, species_dict=subsampled_time_trajectory)
        df_WS_abs["true_N"] = df_NS["true_N"]
        wr.write_table(table=df_WS_abs, file_suffix=file_suffix)
        wr.write_config_yaml(file_suffix=file_suffix)


# write fasta with all relative subsample
if args.sub_rel is not None:
    ts = range(args.t_final)
    for s_rel in args.sub_rel:
        print("--- Subsample sequence set taking " + str(s_rel) + " ---")
        header_prefix = header_prefix=">WS|" + str(s_rel) + "|"
        file_suffix = "WSREL_"+ str(s_rel)
        # total sample set size
        subsample_size = round(sum(df_NS["true_N"]) * s_rel)

        # sampling without replacement from which time point the sequences are coming (weighted by the number seqs)
        time_subset = random.sample(ts, k=subsample_size, counts=df_NS["true_N"])
        subsampled_time_trajectory = []
        for t in ts:
            #sampling the particular sequences for each time point without replacemenet
            seq_subset = random.sample(list(time_trajectory[t].keys()), counts=time_trajectory[t].values(), k=time_subset.count(t))

            # counts sampled sequences
            time_trajectory_sub = {}
            for s in seq_subset:
                time_trajectory_sub[s] = time_trajectory_sub.get(s, 0) + 1
            subsampled_time_trajectory.append(time_trajectory_sub)

        df_WS_rel = wr.write_fasta(file_suffix=file_suffix, header_prefix=header_prefix, species_dict=subsampled_time_trajectory)
        df_WS_rel["true_N"] = df_NS["true_N"]
        wr.write_table(table=df_WS_rel, file_suffix=file_suffix)
        wr.write_config_yaml(file_suffix=file_suffix)