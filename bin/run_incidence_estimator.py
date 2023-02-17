#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
import os
from incidence_estimator_paired import IncidenceEstimator


parser = argparse.ArgumentParser()
parser.add_argument("dir", help="Directory of positional window binnings")

args = parser.parse_args()

# List binning directory
bin_modes = args.dir

estimator = IncidenceEstimator(bin_modes)
estimator.run()
