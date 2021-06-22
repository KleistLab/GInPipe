#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 11:36:24 2020

@authors: Maria Trofimova, Yannick Duport
"""
import pysam
import datetime
import re
import os

'Testing the binning'

class testDIR:
    
    def __init__(self, samfile_dir, outfile_dir):
        self.samfile_dir = samfile_dir
        self.outfile_dir = outfile_dir

    def test_binning(self):
        # List bam files
        filename_bam = os.listdir(self.samfile_dir)
        filenames = [k for k in filename_bam if '.bai' not in k and '.bam' in k]
        filenames.sort()
        dates = []
        
        for i in range(len(filenames)):
            dates.append([])
        
        count = 0
        for filename in filenames:
            file_path = str(self.samfile_dir / filename)
            pysam.index(file_path)
            samfile = pysam.AlignmentFile(file_path, "rb")
            for read in samfile.fetch("EMBOSS_001"):
                name = read.query_name
                date_ = re.search('[|]20[\d]{2}\-[\d]{2}\-[\d]{2}', name)
                if date_:
                    date = re.search('[|]20[\d]{2}\-[\d]{2}\-[\d]{2}', name)[0][1:]
                    if len(date)==7:
                        date = date+"-15"
                    dates[count].append((date, name))
            count += 1
        
        for date_bin in dates:
            date_bin.sort(key=lambda x: x[0])
        
        date_format = "%Y-%m-%d"
        
        strptime = datetime.datetime.strptime
        
        sample = open(self.outfile_dir, 'w') 
        for i in range(len(dates)):
            if dates[i]!=[]:
                start_day = strptime(dates[i][0][0], date_format).date()
                end_day = strptime(dates[i][-1][0], date_format).date()
                print ('Start bin on %s and end on %s' % (start_day, end_day), file = sample)
                print('\tTotal number of days: ',end_day-start_day, file = sample)
        sample.close()