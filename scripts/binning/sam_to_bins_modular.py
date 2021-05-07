#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  3 17:19:40 2020

@authors: Yannick Duport, Maria Trofimova

Getting base counts with pysam and binning routine 
"""
import csv
import datetime
from datetime import date
import math
import numpy as np
import pandas as pd
import os
import pysam
import re
from warnings import warn
strptime = datetime.datetime.strptime


class SAM:
    def __init__(self, samfile_path, bins_dir, meta_dir, num_per_bin, days_per_bin, seq_name):
        self.samfile_path = samfile_path
        self.bins_dir = bins_dir
        self.meta_dir = meta_dir
        self.seq_name = seq_name
        self.samfile = pysam.AlignmentFile(self.samfile_path, "rb")
        self.header = self.samfile.header
        self.sam_dict = self._index_dates()
        self.eq_num_param = num_per_bin
        self.eq_days_param = days_per_bin
        #workaround to make it faster: create for unique sequences a map to its index
        self.seq_idx_dict = self._create_sequence_map()

    def _to_dict(self, l):
        sam_dict = {}
        for i in range(len(l)):
            idx, date, name = l[i]
            sam_dict[i] = {
                'index': idx,
                'date': date,
                'header': name
            }
        return sam_dict

    def _create_sequence_map(self):
        '''
        WORKAROUND: make binning faster, by mapping for each unique sequence its index
        :return:
        '''
        # crate map of each sequence where to find it
        seq_dict = {}
        for idx, attr in self.sam_dict.items():
            if not seq_dict.get(attr['header']):
                seq_dict[attr['header']] = idx
        return seq_dict
    
    def _write_dates_meta(self, idx_dates):
        filename = "meta_dates.tsv"
        os.makedirs(self.meta_dir, exist_ok=True)
        filepath = str(self.meta_dir)+'/'+filename
        with open(filepath, 'w+', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t',
                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
            writer.writerow(["Sequence_ID","Collection_date"])
            for i in range(len(idx_dates)):
                writer.writerow([idx_dates[i][2],idx_dates[i][1]])
    
    def _index_dates(self):
        """Extracts headers, dates and index from a samfile
        Extracts headers, dates, index from samfile, sorts them by date and adds them to a dictionary
        Dictionary keys contain the index from the sorted dates, values are header, date, index
        :return: dict[sorted_index]={header:string, date:string, index:int}
        """
        idx_dates = []
        index = 0
        for read in self.samfile.fetch(self.seq_name):
            name = read.query_name
            # date_format = "[|]20[\d]{2}\-[\d]{2}\-[\d]{2}"
            #date_format = r"[|]20[\d]{2}\-*\d*-*\d*"
            # adding negative look behind to allow only for exact 4 digits for year and not more
            date_format = r"[|]20[\d]{2}(?!\d)\-*\d*-*\d*"
            dt = re.search(date_format, name)
            if dt:
                dt = dt[0][1:]
                if len(dt) == 7:    # dates without days (format YYYY-mm)
                    # FIXME: include date in case for certain binnings (e.g. by month)
                    # dt = dt+"-15"
                    continue
                elif len(dt) == 4:  # dates without months and years (format YYYY)
                    # dt = f"{dt}-02-15"
                    continue
                try:
                    Date = strptime(dt, "%Y-%m-%d").date()
                except ValueError:
                    Date = strptime(dt, "%Y-%d-%m").date()
                    dt = f"{dt[0:5]}{dt[8:]}{dt[4:7]}"
                    warn(f"Date-format different from '%Y-%m-%d' in: {name}\n"
                         f"'%Y-%d-%m' used instead")
                Today = date.today()
                if Date > Today:
                    dt = f"{dt[0:5]}{dt[8:]}{dt[4:7]}"
                idx_dates.append((index, dt, name))
                index += 1
            else:
                warn(f"No date found in the following query name: {name}")
        idx_dates.sort(key=lambda x: x[1])
        self._write_dates_meta(idx_dates)
        sam_dict = self._to_dict(idx_dates)
        return sam_dict

    def _create_filenames(self, folder, n_bins, infix=None):
        """Creates folder, files and returns filenames

        :param method: string
        :param n_bins: int
        :return:
        """
        filepath = self.bins_dir / folder
        basename = os.path.basename(self.samfile_path)
        basename = os.path.splitext(basename)[0]

        os.makedirs(filepath, exist_ok=True)
        filenames_bin = sorted([filepath / f"bin_{basename}_{idx:04d}.bam" for idx in range(n_bins)])
        filenames_header = sorted([filepath / f"header_{basename}_{idx:04d}.tsv" for idx in range(n_bins)])
        filenames_range = sorted([filepath / f"range_{basename}_{idx:04d}.tsv" for idx in range(n_bins)])
        
        return filenames_bin, filenames_header, filenames_range

    def _write_bins_old(self, filenames, indices):
        bin_idx = 0
        for filename in filenames:
            print(f"Writing Bin {bin_idx}")
            names = [self.sam_dict[i]['header'] for i in indices[bin_idx]]
            with pysam.AlignmentFile(filename, "wb", header=self.header) as outfile:
                for read in self.samfile.fetch(self.seq_name):
                    # https://github.com/pysam-developers/pysam/issues/509!!!
                    read.tid = 0
                    if read.query_name in names:
                        outfile.write(read)
            bin_idx = bin_idx + 1

    def _write_bins(self, filenames, seq_bin_dict):
        print(f"Writing to bins")
        for read in self.samfile.fetch(self.seq_name):
            # https://github.com/pysam-developers/pysam/issues/509!!!
            read.tid = 0
            if seq_bin_dict is not None:
                if read.query_name in seq_bin_dict:
                    bin_idx = seq_bin_dict[read.query_name]
                    tmpfile = str(filenames[bin_idx])+".tmp"
                    if not os.path.exists(tmpfile):
                        print("Writing bin " + str(bin_idx))
                        with pysam.AlignmentFile(str(tmpfile), "w", header=self.header) as outfile:
                            outfile.write(read)
                    else:
                        with open(str(tmpfile), 'a') as outfile:
                            outfile.write((read.to_string()+"\n"))
        #Workaround: How to append in a BAM file, or how to print binary into a BAM file?
        print("Convert SAM to BAM")
        for filename in filenames:
            tmpfile = str(filename)+".tmp"
            if os.path.exists(tmpfile):
                infile = pysam.AlignmentFile(tmpfile, "r")
                with  pysam.AlignmentFile(filename, "wb", template=infile) as bamfile:
                    for s in infile:
                        bamfile.write(s)
                os.remove(tmpfile)
            #also create empty bam for downstream stuff
            #else:
             #  open(filename, 'w+')





    def _write_header(self, filenames, indices):
        # uniquify header names
        for i in range(len(indices)):
            headers = []
            idx_list = []
            for idx in indices[i]:
                if self.sam_dict[idx]['header'] not in headers:
                   headers.append(self.sam_dict[idx]['header'])
                   idx_list.append(idx)
            indices[i] = idx_list
        fieldnames = ['header', 'date']
        for i in range(len(filenames)):
            with open(filenames[i], 'w+') as csvfile:
                writer = csv.DictWriter(csvfile,
                                        delimiter='\t',
                                        fieldnames=fieldnames,
                                        extrasaction='ignore')
                writer.writeheader()
                for j in indices[i]:
                    writer.writerow(self.sam_dict[j])
                    
    def _write_days_ranges(self, filenames, days_ranges):
        for i, filename in enumerate(filenames):
            with open(filename, 'w+') as csvfile:
                writer = csv.writer(csvfile,
                                    delimiter='\t')
                writer.writerow(['start_day','end_day'])
                writer.writerow([days_ranges[i][0],days_ranges[i][1]])
    

    def bin_eq_size_names(self):
        """Creates bins of equal size (n=10), but files with same name are treated as one

        :return:
        """
        print(f"Reads per bin: {self.eq_num_param}")

        binsize = self.eq_num_param
        n_reads = len(self.sam_dict)
        names = [self.sam_dict[i]['header'] for i in range(n_reads)]
        names_unique = pd.unique(names)
        n_bins = math.ceil(len(names_unique)/binsize)

        filenames_bin, filenames_header, filenames_range  = self._create_filenames(f"eq_size_names_{binsize}", n_bins)
        #for filename in filenames_bin:
        #    open(filename, 'w+')

        bins_n = [[] for _ in range(n_bins)]
        indices = [[] for _ in range(n_bins)]
        #workaround to make it faster
        seq_bin_dict={}
        for i in range(len(names_unique)):
            bin_id = math.floor(i / binsize)
            bins_n[bin_id].append(names_unique[i])
            seq_bin_dict[names_unique[i]] = bin_id
        
        for i in range(len(bins_n)):
            #for ii in range(len(bins_n[i])):
            #    for j in range(n_reads):
            #        if self.sam_dict[j]['header'] == bins_n[i][ii]:
            #            indices[i].append(j)
            for seq_name in bins_n[i]:
                if seq_name in self.seq_idx_dict:
                    indices[i].append(self.seq_idx_dict[seq_name])

        self._write_bins(filenames_bin, seq_bin_dict)
        self._write_header(filenames_header, indices)

    def bin_eq_days(self, fuzzy=False):
        """Create bins containing equal number of days

        :param fuzzy:
        :return:
        """
        if fuzzy:
            folder = f"fuzzy_days"
        else:
            folder = f"eq_days"

        n_reads = len(self.sam_dict)
        date_format = "%Y-%m-%d"
        start_day = strptime(self.sam_dict[0]['date'], date_format).date()
        end_day = strptime(self.sam_dict[n_reads - 1]['date'], date_format).date()
        n_days = (end_day - start_day).days + 1
        days_per_bin = self.eq_days_param
        n_bins = math.ceil(n_days / days_per_bin)
        
        days_ranges = []
        
        for i in range(n_bins):
            days_ranges.append((start_day+datetime.timedelta(days=(days_per_bin*i)),start_day+datetime.timedelta(days=(days_per_bin*(i+1)-1))))
            

        print(f"total number of days: {n_days}")
        print(f"days per bin: {days_per_bin}")

        # Create empty files
        filenames_bin, filenames_header, filenames_range = self._create_filenames(f"{folder}_{days_per_bin}", n_bins)
        for filename in filenames_bin:
            open(filename, 'w+')

        # fill list of indices and bins
        indices = [[] for _ in range(n_bins)]
        curr_bin = 0
        seq_bin_dict = {}
        for i in range(n_reads):
            curr_day = strptime(self.sam_dict[i]['date'], date_format).date()
            curr_n_days = (curr_day - start_day).days + 1
            if curr_n_days > days_per_bin:
                # r1 = np.random.uniform(low=0.0, high=1.0, size=1)
                # r2 = np.random.uniform(low=0.0, high=1.0, size=1)
                # if fuzzy and (r1 > r2):
                if fuzzy and np.random.randint(2):
                    pass
                else:
                    curr_bin += 1
                    start_day = strptime(self.sam_dict[i]['date'], date_format).date()
            indices[curr_bin].append(i)
            #workaround to make it faster
            seq_bin_dict[self.sam_dict[i]['header']] = curr_bin

        # write bam-files (1 file per bin)
        self._write_bins(filenames_bin, seq_bin_dict)
        self._write_header(filenames_header, indices)
        self._write_days_ranges(filenames_range, days_ranges)
        
 
    def bin_cal_week(self):
        """Binning by calendar week

        :return:
        """
        n_reads = len(self.sam_dict)
        date_format = "%Y-%m-%d"
        start_day = strptime(self.sam_dict[0]['date'], date_format).date()
        end_day = strptime(self.sam_dict[n_reads - 1]['date'], date_format).date()
        start_day_offset = start_day.isocalendar()[2]
        n_days = (end_day - start_day).days
        n_bins = math.ceil((n_days + start_day_offset) / 7)
        
        days_ranges = []
        
        for i in range(n_bins):
            days_ranges.append((start_day+datetime.timedelta(days=(7*i)),start_day+datetime.timedelta(days=(7*(i+1)-1))))
            

        filenames_bin, filenames_header, filenames_range = self._create_filenames('cal_week', n_bins)
        for filename in filenames_bin:
            open(filename, 'w+')

        indices = [[] for _ in range(n_bins)]
        counter = 0
        curr_year, curr_week = strptime(self.sam_dict[0]['date'], date_format).isocalendar()[:2]
        seq_bin_dict = {}
        for i in range(n_reads):
            next_year, next_week = strptime(self.sam_dict[i]['date'], date_format).isocalendar()[:2]
            if next_week > curr_week or next_year > curr_year:
                counter += 1
                curr_year, curr_week = next_year, next_week
            indices[counter].append(i)
            # workaround to make it faster
            seq_bin_dict[self.sam_dict[i]['header']] = counter

        # write bam-files (1 file per bin)
        self._write_bins(filenames_bin, seq_bin_dict)
        self._write_header(filenames_header, indices)
        self._write_days_ranges(filenames_range, days_ranges)
        
