#import gpflow as gp
import numpy as np
import pandas as pd
import utils.date_routines as dr
#from sklearn.neighbors import KernelDensity
#from scipy.ndimage import gaussian_filter1d

import skfda
from skfda import FDataGrid
import skfda.preprocessing.smoothing.kernel_smoothers as ks

#from skfda.misc.hat_matrix import NadarayaWatsonHatMatrix
#from skfda.preprocessing.smoothing.validation import SmoothingParameterSearch
#from skfda.preprocessing.smoothing import KernelSmoother


def smooth_values(x, y, bandwidth):   
    # TODO gpflow macht nur ärger
    # Try from https://fda.readthedocs.io/en/latest/auto_examples/plot_kernel_smoothing.html
    # Exakt meine Vorhaben: https://stackoverflow.com/questions/64700958/python-kernel-smoothing
    # => to get the same result as in ksmooth multiply the smoothing parameter with 0.3706506.... (necessary?)
    # Although ksmooth in R is based on Nadaraya--Watson... (also in this file)
    # Local linear regression kernel smoothing.
    #TODO wieder raus nehmen
    bandwidth = bandwidth*0.3706506
    
#######

    fd = FDataGrid(sample_points=[x],
           data_matrix=[y])
    output_points = np.arange(min(x), max(x)+1)#[:, None] #np.linspace(min(x), max(x)+1)
    smoother = ks.NadarayaWatsonSmoother(smoothing_parameter=bandwidth, 
                                         output_points=output_points)
    smoothed = smoother.fit_transform(fd)

    # choose best smoothing parameter from a range of bandwiths. 
    # TODO: although done like in examples, it throws error that 'SmoothingParameterSearch' object has no attribute 'param_name'???
    # bandwidths = np.arange(1, 31)
    # nw = SmoothingParameterSearch(
    #       KernelSmoother(kernel_estimator=NadarayaWatsonHatMatrix(),
    #                      output_points=output_points),
    #       param_values=bandwidths,
    #       param_name='kernel_estimator__bandwidth'
    #       )   
    # nw.fit(fd)
    # smoothed_opt = nw.transform(fd)

    return  pd.DataFrame({
          'x': output_points,
          'y_smoothed': smoothed.data_matrix.ravel(),
        })

def smooth_phi_estimates(phiPerBin_table, smoothingBandwidth=7):

    #sort by date
    phiPerBin_table = phiPerBin_table.sort_values(by=['t'])

    date_zero = dr.get_date_zero(phiPerBin_table['t'].iloc[0], phiPerBin_table['date'].iloc[0])
    
    smoothedPhi_table = smooth_values(x=phiPerBin_table['t'], y=phiPerBin_table['phi'], bandwidth=smoothingBandwidth)
    
    return pd.DataFrame({'t': smoothedPhi_table.x.astype(int),
		'date': dr.add_days_to_date(smoothedPhi_table.x.to_list(), date_zero),
		'phi_smoothed': smoothedPhi_table.y_smoothed})
		#'phi_lower': smoothedPhi_table.f_lower,
		#'phi_upper': smoothedPhi_table.f_upper}) 