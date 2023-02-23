#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import numpy as np
import math
import re
import datetime
import os
from datetime import date
strptime = datetime.datetime.strptime



class BinningTemporal:
    def __init__(self, table, size):
        self.table = pd.read_table(table, delimiter="\t")
        self.size = size

    def _getReferenceInfo(self):
        """Get reference name, sequence and length
        """
        name = ''
        ref_sequence = ''
        for fasta in self.ref:
            name, ref_sequence = fasta.id, str(fasta.seq)
        list_sequence = list(ref_sequence.upper())
        self.ref_length = len(list_sequence)
        self.ref_seq_list = list_sequence
        # TODO: name without ">"
        self.name = name

    def _makeDateColumn(self):
        names = self.table["name"]
        dates = [n.split("|")[-1] for n in names]

        return dates

    def _dateToDaysSince(self, start, dates):
        days_since = [strptime(d, "%Y-%d-%m").date() - start for d in dates]

        return days_since

    def _writeFile(self, folder, filename, table):
        outfile = '%s/%s' % (folder, filename)
        table.to_csv(outfile, sep='\t', index=False, mode='w')
    
    def binByDays(self):
        # Make date column
        dates = self._makeDateColumn()
        print(self.table)
        # Table sorted by date
        self.table["date"] = dates
        # Sort table by date
        self.table = self.table.sort_values(by=["date"])
        print(self.table["date"])
        date_col = self.table["date"]
        start = strptime(date_col[0], "%Y-%d-%m").date()
        end = strptime(date_col[len(date_col)-1], "%Y-%d-%m").date() 
        n_days = (end - start).days + 1
        days_per_bin = int(self.size)
        n_bins = math.ceil(n_days / days_per_bin)

        # Binning borders
        bins_ranges = np.arange(0,days_per_bin*(n_bins+1),step=days_per_bin)

        days_since = self._dateToDaysSince(start, dates)
        self.table['days_since'] = days_since

        # Get mean read length
        pair_lengths = [e-s for e,s in zip(self.table['end_pair'],self.table['start_pair'])]
        mean_pair_length = np.repeat(math.ceil(np.mean(np.array(pair_lengths))),len(np.array(pair_lengths)))
        self.table['mean_pair_length'] = mean_pair_length

        line_ind = 0
        bin_index = []
        for index, line in self.table.iterrows():
            line_ind += 1
            day = line['days_since'].days
            # Check if current day bigger than bin borders
            bins_bool = [day >= b for b in bins_ranges]
            # Sum  bool list to get bin index
            bin_ind = sum(bins_bool)-1
            bin_index.append(bin_ind)
        self.table["bin_index"] = bin_index

        print('Bins')
        print(n_bins)
        table_grouped = self.table.groupby(['bin_index'])
        print(table_grouped.head(5))
        folder = f"bins/eq_days_{self.size}"
        os.mkdir(folder)
        for i in range(0,n_bins):
            #ind_str = "%s" % i
            bin_table = table_grouped.get_group(i)
            filename = f"eq_days_{self.size}_bin_{i+1}.tsv"
            self._writeFile(folder, filename, bin_table)

        return 0
    
    def binBySize(self):
        # Write lines into bins of equal size line by line
        # Make date column
        dates = self._makeDateColumn()
        # Table sorted by date
        self.table["date"] = dates
        # Sort table by date
        self.table = self.table.sort_values(by=["date"])
        folder = f"bins/eq_size_{self.size}"
        os.mkdir(folder)
        for i, df in enumerate(np.array_split(self.table, math.ceil(len(self.table.index)/int(self.size)))):
            filename = f"eq_size_{self.size}_bin_{i+1}.tsv"
            self._writeFile(folder, filename, df)

        return 0

    
