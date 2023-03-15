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
        dates = [strptime(re.search("\d{4}-\d{2}-\d{2}",n).group(0),"%Y-%m-%d").date() for n in names]

        return dates

    def _dateToDaysSince(self, start, dates):
        days_since = [d - start for d in dates]

        return days_since

    def _writeFile(self, folder, filename, table):
        outfile = '%s/%s' % (folder, filename)
        table.to_csv(outfile, sep='\t', index=False, mode='w')
    
    def binByDays(self):
        # Make date column
        # TODO - incorrect date sorting? 
        dates = self._makeDateColumn()
        # Table sorted by date
        self.table["date"] = dates
        
        # Sort dates
        dates_sort = np.array(sorted(dates))
        start = dates_sort[0]
        end = dates_sort[-1]
        n_days = (end - start).days + 1
        days_per_bin = int(self.size)
        n_bins = math.ceil(n_days / days_per_bin)
        
        # Binning borders
        bins_ranges = np.arange(0,days_per_bin*(n_bins+1),step=days_per_bin)
        # days_siince on unsorted dates array
        days_since = self._dateToDaysSince(start, dates)
        self.table["days_since"] = days_since

        # Get mean read length
        pair_lengths = [e-s for e,s in zip(self.table['end_pair'],self.table['start_pair'])]
        mean_pair_length = np.repeat(math.ceil(np.mean(np.array(pair_lengths))),len(np.array(pair_lengths)))
        self.table["mean_pair_length"] = mean_pair_length

        line_ind = 0
        bin_index = []
        for index, line in self.table.iterrows():
            line_ind += 1
            day = line['days_since'].days
            # Check if current day bigger than bin borders
            bins_bool = [day > b for b in bins_ranges]
            # Sum  bool list to get bin index
            bin_ind = sum(bins_bool)
            bin_index.append(bin_ind)

        self.table["bin_index"] = bin_index
        table_grouped = self.table.groupby("bin_index")
        folder = f"bins/eq_days_{self.size}"
        os.mkdir(folder)
        for ind, group in table_grouped:
            filename = f"eq_days_{self.size}_bin_{int(ind)+1}.tsv"
            self._writeFile(folder, filename, group)

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

    
