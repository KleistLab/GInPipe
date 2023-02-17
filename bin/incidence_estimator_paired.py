#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import math
import csv
import os
import datetime
from scipy import optimize
from scipy.stats import gamma
from scipy.stats import t
from scipy.stats import sem
from scipy.optimize import minimize
from numdifftools import Jacobian, Hessian
from datetime import date
strptime = datetime.datetime.strptime

class IncidenceEstimator:
    def __init__(self, traj_dir):
        self.traj_dir = traj_dir

    def _fmle_short(self,x,nu,ns):
        """
        Estimate function for root finding.

        :param nu: number of origins
        :type nu: int
        :param ns: number of mutant sequences
        :type ns: int
        :param x: param to optimise
        :type x: float
        :returns: population size estimate
        :rtype: float
        """
        return (nu-(x*np.log(1+ns/x)))**2
    
    def _hess(self, x, nu, ns):
        return Hessian(lambda x: self._fmle_short(x,nu,ns))(x)
    
    def _jac(self, x, nu, ns):
        return Jacobian(lambda x: self._fmle_short(x,nu,ns))(x).ravel()

    def _optimize(self, nu, ns):
        """
        Optimize function.

        :param nu: number of origins
        :type nu: int
        :param ns: number of mutant sequences
        :type ns: int
        :returns: optimal population size estimate
        :rtype: float
        """
        # Constraint
        con1 = {'type': 'ineq', 'fun': lambda x: x}
        hess = lambda x: np.zeros((1, 1))
        # If no data available -- set estimate to zero
        if nu==0 or ns==0:
            return 0
        # Else optimize with trust-constr option to keep evaluation in feasible region
        else:
            x0 = 1
            #sol = optimize.minimize(self._fmle, x0=x0, method='trust-constr', args=(nu,ns), constraints=con1)
            # use least squares without MLE
            #sol = optimize.minimize(self._fmle_short, x0=x0, method='trust-constr', args=(nu,ns), constraints=con1)
            sol_f = minimize(self._fmle_short, x0, args=(nu,ns,), method='dogleg', jac=self._jac, hess=self._hess)
            return sol_f.x[0]

    def _countHaplotypesAndMutants(self, table):
        fps = set(zip(table['fp1'],table['fp2']))
        n_mutants = sum([str(a)!='nan' and str(a)!='nan' for a, b in zip(table['fp1'],table['fp2'])])

        return fps, n_mutants


    def _makePositionalEstimate(self, dir, filenames):
        """Estimate incidence for a single positional window table
        """
        estimates = []
        start_dates = []
        end_dates = []
        for filename in filenames:
            table_path = '%s/%s' % (dir, filename)
            table = pd.read_table(table_path, delimiter="\t")
            # Check if empty
            if not table.empty:
                # Get first and last date in the posit window
                table = table.sort_values(by=["date"])
                dates = np.array(table['date'])
                start = strptime(dates[0], "%Y-%d-%m").date()
                start_dates.append(start)
                end = strptime(dates[-1], "%Y-%d-%m").date()
                end_dates.append(end)

                #fp_table = pd.DataFrame(list(zip(table['fp1'], table['fp2'])),columns =['fp1', 'fp2'])
                haplos, n_mutants = self._countHaplotypesAndMutants(table)
                print(len(haplos),n_mutants)
                sol = self._optimize(len(haplos), n_mutants)
                estimates.append(sol)

        return estimates, start_dates, end_dates

    def _confidenceInterval(self, vals):
        interval = t.interval(alpha=0.95, df=len(vals)-1,
              loc=np.mean(vals),
              scale=sem(vals))
        return interval

    def run(self):
        """Combine positional window estimate for each time point into 
        a trajectory
        """
        # Start with binning mode level 
        modes = os.listdir(self.traj_dir)
        for mode in modes:
            filename = '%s_%s.tsv' % (mode,'incidence')
            with open(filename, 'w+') as outfile:
                outwriter = csv.writer(outfile, delimiter='\t')
                # posit_binning_date_start/end are the first and last dates encountered
                # in all positional bins at time point of estimate
                header = ['posit_binning_date_start','posit_binning_date_start','estimate','lower_05','upper_05','min','max']
                outwriter.writerow(header)
                times_dir = '%s/%s' % (self.traj_dir,mode)
                list_times_dir = os.listdir(times_dir)
                
                for time_point in list_times_dir:
                    pos_dir = '%s/%s' % (times_dir,time_point)
                    pos_filenames = os.listdir(pos_dir)

                    # Collect estimates
                    estimates, start_dates, end_dates = self._makePositionalEstimate(pos_dir,pos_filenames)
                    vals = np.array(estimates)
                    interval = self._confidenceInterval(vals)
                    print(min(start_dates))
                    mean = np.mean(vals)
                    minimum = np.min(vals)
                    maximum = np.max(vals)
                    line = [min(start_dates),max(end_dates),mean,interval[0],interval[1],minimum,maximum]
                    outwriter.writerow(line)

        return 0